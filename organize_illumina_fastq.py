#! /usr/bin/env python

## last updated 4/4/2022
## update makes it so it splits file name on "_L" for nexseq runs this will be less restrictive than "_L00"

import os
import glob
import re
import shutil
import pandas as pd 
from datetime import date
import time
import sys

import argparse
from google.cloud import storage
import subprocess


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--seq_run", help="sequencing run name; should match your current working directory")
    parser.add_argument("--bucket_path", help="google cloud bucket path; e.g. gs://covid_terra")
    parser.add_argument("--read_type",  help="either 'paired' or 'single'")
    parser.add_argument("--sample_sheet", help='the sequencing excel (.xlsx) workbook from the j drive', default = 'not provided')
    parser.add_argument('--terra_output_dir', help='(optional) default = gs://covid_terra/{seq_run}/terra_outputs/', default = 'not provided')
    parser.add_argument('--remove_bs_dir', help = '(optional) either TRUE or FALSE; default = TRUE and the bs directories will be removed', default = 'TRUE')
    options = parser.parse_args(args)
    return options


def read_in_sample_sheet(sample_sheet_path, seq_run, output_dir, terra_output_dir, tech_platform, download_date, read_type, bucket_path):
    print('  *** reading in sample sheet')
    
    # open excel sequencing workbook and idnetify the sample sheet tab
    xl = pd.ExcelFile(sample_sheet_path)
    sheet_names = xl.sheet_names
    
    if 'Sample Sheet' in sheet_names: # NEXSEQ runs
        target_sheet = 'Sample Sheet'
    elif 'SampleSheet' in sheet_names: # MISEQ runs
        target_sheet = 'SampleSheet'
        
    # save teh sheet as its own tsv
    xl = pd.read_excel(sample_sheet_path, sheet_name = target_sheet)
    target_sheet_tsv = os.path.join(output_dir, '%s_sample_sheet.tsv' % seq_run)
    xl.to_csv(target_sheet_tsv, sep = '\t', index = False)
    
    # figure out what row the header begins
    count = -1
    with open(target_sheet_tsv, 'r') as f:
        lines = f.readlines()
        for line in lines:
            count += 1
            if re.search('Sample_ID', line.strip()):
                header_line = count
                break
    
    df = pd.read_csv(target_sheet_tsv, sep = '\t', dtype = {'Sample_ID' : object}, header = header_line)
    df = df.rename(columns = {'Sample_ID': 'accession_id'})
    df = df.dropna(subset = ['accession_id'])
    df = df.reset_index(drop = True)
    columns = df.columns
    
#     if tech_platform == 'Illumina NovaSeq':
#         for row in range(df.shape[0]):
#             sample_name = df.accession_id[row]
#              new_sample_name = sample_name.replace('_', '')
#              df.at[row, 'accession_id'] = new_sample_name
    
    # add some additional columns to df
    df['out_dir'] = terra_output_dir
    df['seq_run'] = seq_run
    df['download_date'] = download_date
    df['tech_platform'] = tech_platform
    df['read_type'] = read_type
    
    # check for primer_set
    if not 'primer_set' in df.columns:
        print('  .... primer_set not included in sample sheet; primer_set will be recorded as "not specified"')
        df['primer_set'] = 'not specified'
    
    # check for plate and sample well location in columns (different between NEXSEQ and MISEQ...urgh...)
    if 'COVSEQ_Plate' in columns: #NEXSEQruns
        df = df.rename(columns = {'COVSEQ_Plate' : 'plate_name', 'Well_Location' : 'plate_sample_well'})
    
    elif 'Sample_Plate' in columns: # MISEQruns
        df = df.rename(columns = {'Sample_Plate' : 'plate_name', 'Sample_Well' : 'plate_sample_well'})
    
    else:
        print('  .... sample_plate and plate_sample_well are not included in the sample sheet; will be recorded as "not provided"')
        df['plate_name'] = 'not provided'
        df['plate_sample_well'] = 'not provided'
        
   
    # reorder the columns
    if tech_platform == 'Illumina NovaSeq':
        col_order = ['accession_id', 'Lane','primer_set', 'seq_run', 'download_date','tech_platform', 'read_type', 
                  'plate_name', 'plate_sample_well', 'out_dir']
    
    else:
        col_order = ['accession_id', 'primer_set', 'seq_run', 'download_date','tech_platform', 'read_type', 
                  'plate_name', 'plate_sample_well', 'out_dir']
    
    df = df[col_order]
    df = df.reset_index(drop = True)
    
    # remove the temp tsv
    os.remove(target_sheet_tsv)
    
    return df

def get_sample_list(df, tech_platform):
    if tech_platform == 'Illumina NovaSeq':
        sample_list = []
        # combine the lane info with the sample id....
        for row in range(df.shape[0]):
            accession_id = df.accession_id[row]
            lane = df.lane[row]
            
            sample_name = '%s_00%d' % (accession_id, lane)
            sample_list.append(lane)
        print('  .... there are %d samples in this run' % len(sample_list))
    
    else:
        sample_list = df.accession_id.unique().tolist()
        print('  ..... there are %d samples in this run' % len(sample_list))
    
    return sample_list


def create_fastq_files_directory(seq_run, current_dir):
    print('')
    print('  *** creating fastq files directory')
    
    fastq_files_dir = os.path.join(current_dir, 'fastq_files')

    if os.path.exists(fastq_files_dir):
        shutil.rmtree(fastq_files_dir)
    
    os.mkdir(fastq_files_dir)
                             
    return fastq_files_dir    
    

def concat_fastq_gz_files_single(current_dir, fastq_files_dir, sample_list, tech_platform):
    print('')
    print('  *** organizing fastq files')
    print('  ..... read_type == single')
    print('  ..... tech_platform == %s' % tech_platform)
    
    if tech_platform == 'Illumina NextSeq':
        print('  ..... each sample should have one fastq file per 4 lanes for a total of 4 fastq files')
        print('  ..... will concatenate 4 fastq files into one')
        print('')
        time.sleep(2)
             
        # concatenate using teh subprocess command
        #### track which samples do not have 4 fastq files
        n =0
        samples_without_4_fastq_files = {}
        for sample in sample_list:
            # track progress
            n = n + 1
            mult_50 = n % 50
            if mult_50 == 0 or n == len(sample_list) or n == 1:
                print('  ..... concatenating %d of %d' % (n, len(sample_list)))
                      
            # create empty list to hold fastq files associated with sample
            #### track which samples don't have 4 fastq files
            file_list = []
            for file in glob.glob('*/*.fastq.gz'):
                sample_name = file.split('_L')[0]
                if sample_name == sample:
                    file_list.append(file)
            
            if len(file_list) !=4:
                samples_without_4_fastq_files[sample] = len(file_list)
                
            if len(file_list) != 0:    
                # concatenated the fastq files in file list
                combined_fastq_file_name = '%s.fastq.gz' % sample
                outpath = os.path.join(fastq_files_dir, combined_fastq_file_name)

                file_list_string = ' '.join(file_list)
                shell_command = 'cat %s > %s' % (file_list_string, outpath)
                subprocess.run(args = shell_command, shell = True, check = True)
            
            # delete the files from the file list
            for file in file_list:
                if os.path.exists(file):
                    os.remove(file)
            
        # print warning to screen if samples have more or less than 4 fastq files
        if len(samples_without_4_fastq_files) > 0:
            print('')
            print('  WARNING!')
            print('  ..... the following fastq files do not have 4 fastq files; contact Alex; sequencing run may need to be requeued in bs')
            for sample in samples_without_4_fastq_files:
                print('  ..... ..... %s: %d fastq files' % (sample, samples_without_4_fastq_files[sample]))
            time.sleep(5)   

    elif tech_platform == 'Illumina NovaSeq':
        print('  ..... each sample sample will only have one fastq file')
        print('  ..... will move and rename each fastq file')
        time.sleep(2)
           
            
        # move and rename file...
        # rename fastq file so name matches sample name
        
        # get a list of all fastq files so can do a quick count
        fastq_file_list = []
        for file in glob.glob('*/*.fastq.gz'):
            fastq_file_list.append(file)
        print('  .... there are %d fastq files listed in the directory' % len(fastq_file_list))
        
        # copy to fastq_files_dir
        n = 0    
        for file in glob.glob('*/*.fastq.gz'):
            n = n +1
            mult_100 = n % 100
            if mult_100 ==0 or n == len(fastq_file_list) or n == 1:
                print('  ..... moving and renaming %d of %d fastq files' % (n, len(fastq_file_list)))

            sample_name = file.split('_ds')[0]   
            new_file_name = '%s.fastq.gz' % sample_name
            old_file_name = file.split('/')[-1]
            

            # move the file to save room
            shutil.move(file, fastq_files_dir)

            # rename the file
            old= os.path.join(fastq_files_dir, old_file_name)
            new = os.path.join(fastq_files_dir, new_file_name)
            os.rename(old, new)

        # check that each sample has forward and reverse fastq file
        file_list = os.listdir(fastq_files_dir)
        samples_without_fastq_files = {}
        for sample in sample_list:
            num_fastq_files = 0
            for file in file_list:
                if re.search(sample, file):
                    num_fastq_files = num_fastq_files + 1
            if num_fastq_files != 1:
                samples_without_fastq_files[sample] = num_fastq_files
        
        # print warning to screen if samples have more or less than 4 fastq file
        if len(samples_without_fastq_files) > 0:
            print('')
            print('  WARNING!')
            print('  ..... the following samples do not have a fastq file')
            for sample in samples_without_fastq_files:
                print(' ..... ..... %s: %d fastq files' % (sample, samples_without_fastq_files[sample]))
            time.sleep(5)              
        
    
def copy_and_rename_fastq_gz_files_paired(seq_run, current_dir, fastq_files_dir, sample_list, tech_platform):
    print('')
    print('  *** organizing fastq files')
    print('  ..... read_type == paired')
    print('  ..... tech_platform == %s' % tech_platform)
    print('  ..... each sample will have a fastq file for the forward and reverse reads')
    time.sleep(2)
     
    # get a list of all fastq files so can do a quick count
    fastq_file_list = []
    for file in glob.glob('*/*.fastq.gz'):
        fastq_file_list.append(file)

    # rename fastq file so name matches sample name
    # copy to fastq_files_dir
    n = 0    
    for file in glob.glob('*/*.fastq.gz'):
        n = n +1
        mult_25 = n % 25
        if mult_25 ==0 or n == len(fastq_file_list) or n == 1:
            print('  ..... moving and renaming %d of %d fastq files' % (n, len(fastq_file_list)))
        
        sample_name = file.split('_L00')[0]
        sequencing_info = re.findall('_R\d_\d\d\d', file)[0]   
        new_file_name = '%s%s.fastq.gz' % (sample_name, sequencing_info)
        old_file_name = file.split('/')[-1]
        
        # copy the file
        shutil.move(file, fastq_files_dir)
        
        # rename the file
        old= os.path.join(fastq_files_dir, old_file_name)
        new = os.path.join(fastq_files_dir, new_file_name)
        os.rename(old, new)
        
    # check that each sample has forward and reverse fastq file
    file_list = os.listdir(fastq_files_dir)
    samples_without_2_fastq_files = {}
    for sample in sample_list:
        num_fastq_files = 0
        for file in file_list:
            if re.search(sample, file):
                num_fastq_files = num_fastq_files + 1
        if num_fastq_files != 2:
            samples_without_2_fastq_files[sample] = num_fastq_files
            
    # print warning to screen if samples have more or less than 4 fastq file
    if len(samples_without_2_fastq_files) > 0:
        print('')
        print('  WARNING!')
        print('  ..... the following samples do not have 2 fastq files')
        for sample in samples_without_2_fastq_files:
            print(' ..... ..... %s: %d fastq files' % (sample, samples_without_2_fastq_files[sample]))
        time.sleep(5)              


def push_fastq_files_to_bucket(seq_run, current_dir, fastq_files_dir, bucket_path):
    print('')
    print('  *** pushing fastq_files directory to google bucket')
    print('  ..... this might take a minute or two or three.... or four.... or 5 thousand')
    
    #append the seq run to the bucket path
    full_bucket_path = os.path.join(bucket_path, seq_run, 'fastq_files')
    
    # use gsutil command and supress output to screen
    shell_command = 'gsutil -m cp -r %s/* %s' % (fastq_files_dir, full_bucket_path)
    subprocess.run(args = shell_command, shell = True, check = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    
def create_terra_data_table(df, bucket_path, read_type, seq_run):
    print('')
    print('  *** creating terra datatable and pushing to bucket')
    print('  ..... read_type = %s' % read_type)
    time.sleep(2)
    
    df = df.reset_index(drop = True)
    
    if read_type == 'single':
        for row in range(df.shape[0]):
            accession_id = df.accession_id[row]
            fastq_file_name = '%s.fastq.gz' % accession_id
            fastq_bucket_path = os.path.join(bucket_path, seq_run, 'fastq_files', fastq_file_name)
            df.at[row, 'fastq'] = fastq_bucket_path 
            
        # rename the accession_id column with entity:sampleXXXX_id
        ### find the base of seq_run name
        seq_run_base = seq_run.split('_')[0]
        seq_run_suffix = seq_run.split('_')[1]
        entity_col_header = 'entity:sample%s%s_id' % (seq_run_base, seq_run_suffix)
        df = df.rename(columns = {'accession_id' : entity_col_header})
        
        # save the data table
        print('  ..... saving terra datatable')
        outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
        df.to_csv(outfile, index = False, sep = '\t')
        
        # push to bucket
        print('  ..... pushing terra datatable to bucket')
        #append the seq run to the bucket path
        full_bucket_path = os.path.join(bucket_path, seq_run)

        # use gsutil command and supress output to screen
        shell_command = 'gsutil -m cp -r %s %s' % (outfile, full_bucket_path)
        subprocess.run(args = shell_command, shell = True, check = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    

    elif read_type == 'paired':
        
        for row in range(df.shape[0]):
            accession_id = df.accession_id[row]
            fastq_file_name_R1 = '%s_R1_001.fastq.gz' % accession_id
            fastq_file_name_R2 = '%s_R2_001.fastq.gz' % accession_id
            fastq_bucket_path_R1 = os.path.join(bucket_path, seq_run, 'fastq_files', fastq_file_name_R1)
            fastq_bucket_path_R2 = os.path.join(bucket_path, seq_run, 'fastq_files', fastq_file_name_R2)
            df.at[row, 'fastq_1'] = fastq_bucket_path_R1
            df.at[row, 'fastq_2'] = fastq_bucket_path_R2
              
        # rename the accession_id column with entity:sampleXXXX_id
        ### find the base of seq_run name
        seq_run_base = seq_run.split('_')[0]
        seq_run_suffix = seq_run.split('_')[1]
        entity_col_header = 'entity:sample%s%s_id' % (seq_run_base, seq_run_suffix)
        df = df.rename(columns = {'accession_id' : entity_col_header})

        # save the data table
        print('  ..... saving terra datatable')
        outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
        df.to_csv(outfile, index = False, sep = '\t')
        
        # push to bucket
        print('  ..... pushing terra datatable to bucket')
        #append the seq run to the bucket path
        full_bucket_path = os.path.join(bucket_path, seq_run)

        # use gsutil command and supress output to screen
        shell_command = 'gsutil -m cp -r %s %s' % (outfile, full_bucket_path)
        subprocess.run(args = shell_command, shell = True, check = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        

def remove_fastq_gz_directories(current_dir):
    print('')
    print('  ***removing fastq.gz directories')
    os.chdir(current_dir)
    for directory in os.listdir():
        if re.search('_ds.', directory):
            shutil.rmtree(directory)    
              
              
if __name__ == '__main__':
    
    print('')
    print('  *************************************************************************')
    print('  *** starting ORGANIZE_ILLUMNA_FASTQ')
    print('      .... last updated (major) 2021-10-28')
    print('      .... lastest update revamps script, now sample sheet dependent')
    print('')
    print('')
    
    time.sleep(2) # delay 2 seconds so can read output to screen
    
    
    options = getOptions()
    
    # create variables from user input
    current_dir = os.getcwd()
    seq_run = options.seq_run
    seq_run_prefix = seq_run.split('_')[0]
    if seq_run_prefix == 'COVSEQ':
        tech_platform = 'Illumina MiSeq'
    elif seq_run_prefix == 'NEXSEQ':
        tech_platform = 'Illumina NextSeq'
    elif seq_run_prefix == 'NOVASEQ':
        tech_platform = 'Illumina NovaSeq'
    
    read_type = options.read_type
    sample_sheet_path = options.sample_sheet
    if not re.search('.xlsx', sample_sheet_path):
        print('  ERROR!: sample sheet not formatted correctly, must be excel workbook')
        print('  .... be sure to use the excel sequencing workbook, not the csv import sheet')
        print('  .... exiting')
        print('')
        print('')
        sys.exit()
        
    bucket_path = options.bucket_path
    download_date = str(date.today())
    
    terra_output_dir = options.terra_output_dir
    if terra_output_dir == 'not provided':
        terra_output_dir = os.path.join('gs://covid_terra', seq_run, 'terra_outputs')
    else:
        terra_output_dir = os.path.join(terra_output_dir, seq_run, 'terra_outputs')
        
    remove_bs_dir = options.remove_bs_dir
     
    
    print('  User Inputs and Parameters:')
    print('  current_working_directory: %s' % current_dir)
    print('  download_date: %s' % download_date)
    print('  seq_run: %s' % seq_run)
    print('  sequencing_tech_platform: %s' % tech_platform)
    print('  bucket_path: %s' % bucket_path)
    print('  read_type: %s' % read_type)
    print('  sample_sheet: %s' % sample_sheet_path)
    print('  terra_output_bucket_directory: %s' % terra_output_dir)
    print('  remove_bs_directories: %s' % remove_bs_dir)
    print('')
    print('')
    time.sleep(6)
    
    # okay now run the functions
    df = read_in_sample_sheet(sample_sheet_path = sample_sheet_path, 
                              seq_run = seq_run, 
                              output_dir = current_dir, 
                              terra_output_dir = terra_output_dir, 
                              tech_platform = tech_platform, 
                              download_date = download_date, 
                              read_type = read_type,
                              bucket_path = bucket_path)
    sample_list = get_sample_list(df = df, tech_platform = tech_platform)
    fastq_files_dir = create_fastq_files_directory(seq_run = seq_run,  
                                                   current_dir = current_dir)
    
    # split organizing fastq files depending on read type
    if read_type == 'single':
        concat_fastq_gz_files_single(current_dir = current_dir, 
                                     fastq_files_dir = fastq_files_dir, 
                                     sample_list = sample_list, 
                                     tech_platform = tech_platform)

    elif read_type == 'paired':
        copy_and_rename_fastq_gz_files_paired(seq_run = seq_run, 
                                              current_dir = current_dir, 
                                              fastq_files_dir = fastq_files_dir, 
                                              sample_list = sample_list, 
                                              tech_platform = tech_platform)
        
    
    push_fastq_files_to_bucket(seq_run = seq_run, 
                               current_dir = current_dir, 
                               fastq_files_dir = fastq_files_dir, 
                               bucket_path = bucket_path)
    create_terra_data_table(df =df, 
                            seq_run = seq_run,
                            bucket_path = bucket_path, 
                            read_type = read_type)
    
    if remove_bs_dir == 'TRUE':
        remove_fastq_gz_directories(current_dir)
        
    print('')
    print('  *** DONE!')
    print('')
    print('')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    