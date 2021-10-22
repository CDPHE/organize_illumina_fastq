#! /usr/bin/env python

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


# before running this script....

# 1. create a directory for the sequence run

# 2. download the fastq files from bs into a directory named the same as the seq_run
# for example if you are downloading COVSEQ_0043, the directory must be named 'COVSEQ_0043'

# 3. cd to the directory with the downloaded fastq files (right this directory will 
# have a series of directories with each sample's fastq.gz files)

# 4. run this script from inside the directory
# note this script takes 4 arguments

## Example input ##
# organiz_fastq_files.py --seq_run <seq_run_name> --bucket_path <gs://path_to_bucket> --run_type <paired or single> --sample_sheet {path to sample sheet.xlxs}
# organize_fastq_files.py --seq_run COVSEQ_0050rr --bucket_path gs://covid_terra --run_type paired --sample_sheet some_name.xlxs



#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--seq_run", help="sequencing run name; must match directory")
    parser.add_argument("--bucket_path", help="google cloud bucket path; gs://covid_terra")
    parser.add_argument("--run_type",  help="either 'paired' or 'single'")
    parser.add_argument("--sample_sheet", help='the import .xlsx sheet from the j drive with the plate well location', default = 'not provided')
    parser.add_argument('--terra_output_dir', help='optional, default = gs://covid_terra/{seq_run}/terra_outputs/', default = 'not provided')
    options = parser.parse_args(args)
    return options


def read_in_sample_sheet(sample_sheet, seq_run, output_dir, terra_output_dir):
    print('')
    print('')
    print('  ******* reading in sample sheet and getting plate location info *******')
    xl = pd.ExcelFile(sample_sheet)
    sheet_names = xl.sheet_names

    if 'Sample Sheet' in sheet_names: # NEXSEQ runs
        target_sheet = 'Sample Sheet'
    elif 'SampleSheet' in sheet_names: # MISEQ (WWT/COVSEQ)
        target_sheet = 'SampleSheet'

    # save the sheet as is own tsv
    xl = pd.read_excel(sample_sheet, sheet_name = target_sheet)
    target_sheet_tsv = os.path.join(output_dir, '%s_sample_sheet.tsv' % seq_run)
    xl.to_csv(target_sheet_tsv, sep = '\t', index = False)


    count = -1
    with open(target_sheet_tsv, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            count += 1
            if re.search('Sample_ID', line.strip()):
                header_line = count
                break
    
    df = pd.read_csv(target_sheet_tsv, sep = '\t', dtype = {'Sample_ID' : object}, header = header_line)
    columns = df.columns
    
    # check if primer_set in columns:
    if not 'primer_set' in df.columns:
        print('  ... primer set not specificed; primer_set will be recorded as "not specified"')
        df['primer_set'] = 'not specified'
    
    if "COVSEQ_Plate" in columns: # NEXSEQ runs
        col_order = ['Sample_ID', 'COVSEQ_Plate', 'Well_Location', 'primer_set']
        col_rename = {'Sample_ID' : 'accession_id', 'COVSEQ_Plate' : 'plate_name', 'Well_Location': 'plate_sample_well'}
        df = df[col_order]
        df = df.rename(columns = col_rename)
        sample_list = df.accession_id.unique().tolist()

    elif "Sample_Plate" in columns: # MISEQ (WWT/COVSEQ)
        col_order = ['Sample_ID', 'Sample_Plate', 'Sample_Well', 'primer_set']
        col_rename = {'Sample_ID' : 'accession_id', 'Sample_Plate' : 'plate_name', 'Sample_Well': 'plate_sample_well'}
        df = df[col_order]
        df = df.rename(columns = col_rename)
        sample_list = df.accession_id.unique().tolist()
        
    else:
        print('  ... sample sheet does not contain columns with plate name and sample location')
        print('  ... will continue by specifying plate name and sample location as "not provided"')
        time.sleep(10)
        col_order = ['Sample_ID']
        col_rename= {'Sample_ID' : 'accession_id'}
        df = df[col_order]
        df = df.rename(columns = col_rename)
        df['plate_name'] = 'not provided'
        df['plate_sample_well'] = 'not provided'
        
        col_order = ['Sample_ID', 'Sample_Plate', 'Sample_Well', 'primer_set']
        df = df[col_order]
        
        sample_list = df.accession_id.unique().tolist()        
        
    
    # add out_dir and seq_run columns to df
    df['out_dir'] = terra_output_dir
    df['seq_run'] = seq_run
    
    
    # remove temporary target sheet tsv
    os.remove(target_sheet_tsv)
    
    return {'df': df, 'sample_list' : sample_list}
    
def create_fastq_files_directory(seq_run, current_dir):
    print('')
    print('')    
    print('  ******* creating fastq_files directory *******')
    time.sleep(5)
    
    fastq_files_dir = os.path.join(current_dir, 'fastq_files')

    if os.path.exists(fastq_files_dir):
        shutil.rmtree(fastq_files_dir)
    
    os.mkdir(fastq_files_dir)
                             
    return fastq_files_dir



def concate_fastq_gz_files_single(current_dir, fastq_files_dir, sample_sheet_sample_list):
    print('')
    print('')
    print('  ******* concatenating 4 fastq files per sample id *******')
    time.sleep(5)
    
    os.chdir(current_dir)
    
    # get a list of the sample ids
    sample_list = []
    for file in glob.glob('*_L00*/*.fastq.gz'):
        sample_name = file.split('_L00')[0]
        if sample_name not in sample_list:
            sample_list.append(sample_name)
    
    print('  .... there are %d sample ids in %s' % (len(sample_list), seq_run))
    print('')
    time.sleep(5)
    
    # now concatenate using subprocesss command
    n = 0
    samples_without_4_fastq_files = {}
    for sample in sample_list:
        n = n + 1
        mult_50 = n % 50
        if mult_50 == 0 or n == len(sample_list):
            print('  .... concatenated %d of %d' % (n, len(sample_list)))
        
#         print('  ............%d: concatenating files for sample %s' % (n, sample))

        # create empty list to hold fastq files assoicated with sample
        file_list = []
        for file in glob.glob('*_L00*/*.fastq.gz'):

            sample_name = file.split('_L00')[0]
            if sample_name == sample:
                file_list.append(file)

        combined_fastq_file_name = '%s.fastq.gz' % sample
        out_path = os.path.join(fastq_files_dir, combined_fastq_file_name)

#         print('  concatenating the following files into single file named: %s' % (combined_fastq_file_name))
        i = 0
        for file_name in file_list:
            i = i + 1
            fastq_file_name = file_name.split('/')[1]
#             print('  %d: %s' % (i,fastq_file_name))
        file_list_string = ' '.join(file_list)

        shell_command = 'cat %s > %s' % (file_list_string, out_path)

        subprocess.run(args= shell_command, shell= True, check = True)
#         print('')
        if i != 4:
            samples_without_4_fastq_files[sample] = i
    
    if len(samples_without_4_fastq_files) > 0:
        print('')
        print('  WARNING!!! check the following samples')
        print('  these samples had either less than or more than 4 fastq files; sequencing run may need to be requeued in basespace after removing current fastq files')
        for sample in samples_without_4_fastq_files:
            print('  .... %s: %d fastq files' % (sample, samples_without_4_fastq_files[sample]))
        print('')
        
    if len(sample_sheet_sample_list) > 0: 
        # check that all samples from run have a fastq file:
        missing = []
        for sample in sample_sheet_sample_list:
            if sample not in sample_list:
                missing.append(sample)
        if len(missing) > 0:
            print('  WARNING!! found %d samples in sample sheet without fastq files' % len(missing))
            print('  the missing samples are:')
            for item in missing:
                print('  .... %s' % item)
            print('')
            print('')
            time.sleep(5)

    
    return samples_without_4_fastq_files

def copy_and_rename_fastq_gz_files_paired(seq_run, current_dir, fastq_files_dir, sample_sheet_sample_list):
    print('')
    print('')
    print('  ******* copying and renaming fastq.gz files *******')
    time.sleep(5)
    
    os.chdir(current_dir)
    file_list = os.listdir(current_dir)
    sample_list = []
    num_fastq_files = 0
    for file in glob.glob('*_L001_ds.*/*.fastq.gz'):
        num_fastq_files = num_fastq_files + 1
    n = 0
    for file in glob.glob('*_L001_ds.*/*.fastq.gz'):
        n = n+1
        mult_25 = n % 25
        if mult_25 == 0 or n == num_fastq_files:
            print('  .... renaming %d of %d' % (n, num_fastq_files))
        
        sample_name = file.split('_L001')[0] # the base directory name
        sample_list.append(sample_name)

        sequencing_info = re.findall('_R\d_\d\d\d', file)


        new_file_name = '%s%s.fastq.gz' % (sample_name, sequencing_info[0])
        old_file_name = file.split('/')[1]

        # check that new name makes sense:
#         print('  %d: %s is renamed to %s' % (n, old_file_name, new_file_name))

        # copy file into new directory
        shutil.copy(file, fastq_files_dir)

        # rename the file
        old= os.path.join(fastq_files_dir, old_file_name)
        new = os.path.join(fastq_files_dir, new_file_name)
        os.rename(old, new)
    
    if len(sample_sheet_sample_list) > 0: 
        # check that all samples from run have a fastq file:
        missing = []
        for sample in sample_sheet_sample_list:
            if sample not in sample_list:
                missing.append(sample)
        if len(missing) > 0:
            print('  WARNING!! found %d samples in sample sheet without fastq files' % len(missing))
            print('  the missing samples are:')
            for item in missing:
                print('  ... %s' % item)
            print('')
            print('')
            time.sleep(10)


def upload_fastq_files_to_bucket(seq_run, current_dir, fastq_files_dir, bucket_name, bucket_prefix):
    print('')
    print('')
    print('  ******* uploading fastq_file_directory to google bucket *******')
    print('  .... this might take a minute or two')
    time.sleep(5)
    
#     local_path = fastq_files_dir
    os.chdir(fastq_files_dir)  
    client_storage = storage.Client()
    bucket = client_storage.get_bucket(bucket_name)
    
    n = 0
    for file in glob.glob('*.fastq.gz'):
        n = n + 1
        localFile = os.path.join(fastq_files_dir, file)
        bucket_path = os.path.join(bucket_prefix, file)
        blob = bucket.blob(bucket_path)
        blob.upload_from_filename(localFile)
#         print('  {}: uploaded {} to {} bucket'.format(n, file, os.path.join(bucket_name, bucket_prefix)))
    

          
def create_terra_data_table_single(bucket_path, bucket_name, bucket_prefix, current_dir, seq_run, plate_loc_df, terra_output_dir):
    print('')
    print('')    
    print('  ******* creating datatable for terra *******')
    time.sleep(5)

    seq_run_mod = re.sub('_', '', seq_run)

    df = pd.DataFrame()
    
    sample_name_list = []
    fastq_path_list = []
    
    client_storage = storage.Client()
    blobs = client_storage.list_blobs(bucket_name, prefix = bucket_prefix)
    
    for blob in blobs:
#         print(blob.name)
        if re.search('.fastq.gz', blob.name):
            fastq_file_name = blob.name.split('/')[-1]
            sample_name = fastq_file_name.split('.fastq.gz')[0]
            if re.search('_$', sample_name):
                sample_name = re.sub('_$', '', sample_name)
            sample_name_list.append(sample_name)
          
            file_name = '%s/%s' % (bucket_name, blob.name)
            fastq_path_list.append('gs://%s' % file_name)
        

    df['accession_id'] = sample_name_list
    df['fastq'] = fastq_path_list
    
    if plate_loc_df.shape[0] != 0:
        # add the plate location info to the terra datatable
        df = df.set_index('accession_id')
        plate_loc_df = plate_loc_df.set_index('accession_id')

        df = df.join(plate_loc_df, how = 'left')
        df = df.reset_index()
    else:
        # if no sample sheet is provided
        df['plate_name'] = 'not provided'
        df['plate_sample_well'] = 'not provided'
        df['primer_set'] = 'not specified'
        df['seq_run'] = seq_run
        df['out_dir'] = terra_output_dir
    
    col_order = ['accession_id', 'fastq', 'plate_name', 'plate_sample_well', 'primer_set',  'out_dir', 'seq_run']
    df = df[col_order]
    df = df.rename(columns = {'accession_id' : 'entity:sample%s_id' % seq_run_mod})
   
    outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
    df.to_csv(outfile, index = False, sep = '\t')


    os.chdir(current_dir)
    # upload terra datatable to bucket
    print('')
    print('')
    print('  ******* uploading terra datatable to bucket *******')
    time.sleep(5)
    client_storage = storage.Client()
    bucket = client_storage.get_bucket(bucket_name)
    
#     localFile = os.path.join(fastq_files_dir, file)
    bucket_path = os.path.join(bucket_outerdir_prefix, '%s_terra_data_table.tsv' % seq_run)
    blob = bucket.blob(bucket_path)
    blob.upload_from_filename(outfile)
#     print('  uploaded {} to {} bucket'.format( outfile, os.path.join(bucket_name, bucket_outerdir_prefix)))
    
      
    
def create_terra_data_table_paired(bucket_path, bucket_name, bucket_prefix, current_dir, seq_run, plate_loc_df, terra_output_dir ):
    print('')
    print('')
    print('  ******* creating datatable for terra *******')
    time.sleep(5)

    seq_run_mod = re.sub('_', '', seq_run)
    if re.search('WWT', seq_run_mod):
        seq_run_mod = re.sub('COVSEQ', '', seq_run_mod)
        
    R1_df = pd.DataFrame()
    R2_df = pd.DataFrame()
    
    R1_sample_name_list = []
    R2_sample_name_list = []
    R1_fastq_list = []
    R2_fastq_list = []
    
    client_storage = storage.Client()
    blobs = client_storage.list_blobs(bucket_name, prefix = bucket_prefix)
    
    for blob in blobs:
#         print(blob.name)
        if re.search('R1_001.fastq.gz', blob.name):
            fastq_file_name = blob.name.split('/')[-1]
            sample_name = fastq_file_name.split('_R1_001.fastq.gz')[0]
            if re.search('_$', sample_name):
                sample_name = re.sub('_$', '', sample_name)
            R1_sample_name_list.append(sample_name)
          
            file_name = '%s/%s' % (bucket_name, blob.name)
            R1_fastq_list.append('gs://%s' % file_name)
        
        elif re.search('R2_001.fastq.gz', blob.name):
            fastq_file_name = blob.name.split('/')[-1]
            sample_name = fastq_file_name.split('_R2_001.fastq.gz')[0]
            if re.search('_$', sample_name):
                sample_name = re.sub('_$', '', sample_name)
            R2_sample_name_list.append(sample_name)
          
            file_name = '%s/%s' % (bucket_name, blob.name)
            R2_fastq_list.append('gs://%s' % file_name)

    R1_df['accession_id'] = R1_sample_name_list
    R1_df['fastq_1'] = R1_fastq_list
    R1_df = R1_df.set_index('accession_id')      
    
    
    R2_df['accession_id'] = R2_sample_name_list
    R2_df['fastq_2'] = R2_fastq_list
    R2_df = R2_df.set_index('accession_id')
          
    df = R1_df.join(R2_df, how='left')
    df = df.reset_index()
    
    if plate_loc_df.shape[0] != 0:
        # add the plate location info to the terra datatable
        df = df.set_index('accession_id')
        plate_loc_df = plate_loc_df.set_index('accession_id')

        df = df.join(plate_loc_df, how = 'left')
        df = df.reset_index()
    else:
        df['plate_name'] = 'not provided'
        df['plate_sample_well'] = 'not provided'
        df['seq_run'] = seq_run
        df['out_dir'] = terra_output_dir
        df['primer_set'] = 'not specified'
    
    
    col_order = ['accession_id', 'fastq_1', 'fastq_2', 'plate_name', 'plate_sample_well', 'primer_set', 'out_dir', 'seq_run']
    df = df[col_order]
    df = df.rename(columns = {'accession_id' : 'entity:sample%s_id' % seq_run_mod})
    
    outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
    df.to_csv(outfile, index = False, sep = '\t')
    
    os.chdir(current_dir)
    
    # upload terra datatable to bucket
    print('')
    print('')
    print('  ******* uploading terra datatable to bucket *******')
    time.sleep(5)
    client_storage = storage.Client()
    bucket = client_storage.get_bucket(bucket_name)
    
#     localFile = os.path.join(fastq_files_dir, file)
    bucket_path = os.path.join(bucket_outerdir_prefix, '%s_terra_data_table.tsv' % seq_run)
    blob = bucket.blob(bucket_path)
    blob.upload_from_filename(outfile)
#     print('  uploaded {} to {} bucket'.format( outfile, os.path.join(bucket_name, bucket_outerdir_prefix)))    
    
def remove_fastq_gz_directories(current_dir):
    print('')
    print('')
    print('  ******* removing fastq.gz directories *******')
    os.chdir(current_dir)

    for directory in os.listdir():
        if re.search('_ds.', directory):
            shutil.rmtree(directory)
            
            
            
if __name__ == '__main__':
    
    print('')
    print('  *************************************************************************')
    print('  *** starting ORGANIZE_ILLUMNA_FASTQ ***')
    print('      .... last updated 2021-10-21')
    print('      .... lastest update includes adding plate location to terra datatable')
    print('      .... lastest update includes adding seq_run, out_dir, and primer_set to terra datatable')
    print('')
    print('')
    
    time.sleep(2) # delay 3 seconds so can read output to screen)
    
    
    options = getOptions()
    
    # get teh seq_run name
    seq_run = options.seq_run # this is the sequence run name
    run_type = options.run_type # this is the run type either single or paired
    sample_sheet = options.sample_sheet # path to the excel sample sheet (optional, default is 'not provided')
    terra_output_dir = options.terra_output_dir # path to google bucket for terra output to go (optiona, default is 'not provided' 
                                                # and will be assigned gs://covid_terra/{seq_run}/terra_outputs/
    print('  USER INPUTS')
    print('  sequence run: %s' % seq_run)
    print('  run type: %s' % run_type)
    print('  sample_sheet: %s' % sample_sheet)
    
    # do some magic to create the terra output dir
    if terra_output_dir == 'not provided':
        print('  terra_output_dir: %s' % terra_output_dir)
        print('  terra_output_dir will default to gs://covid_terra/%s/terra_outputs/' % seq_run)
        terra_output_dir = 'gs://covid_terra/%s/terra_outputs/' % (seq_run)
    else:
        print('  terra_output_dir: %s' % terra_output_dir)
        
        if re.search('gs://', terra_output_dir):
            terra_output_dir = terra_output_dir.replace('gs://', '')
        if re.search('/$', terra_output_dir):
            terra_output_dir = re.sub('/$', '', terra_output_dir)
            
        print('  terra_output_dir will become: gs://%s/%s/terra_outputs/' % (terra_output_dir, seq_run))
        terra_output_dir = 'gs://%s/%s/terra_outputs/' % (terra_output_dir, seq_run)
    print('')
    print('')
    time.sleep(3)
    
    if sample_sheet == 'not provided':
        print('  .... no sample sheet provided')
        print('  .... plate locations will be recorded as "not provided"')
        print('  .... primer set will be recorded as "not provided"')
        print("  .... if you'd like the plate location and/or primer_set to be included in the sequencing results output, specify --sample_sheet flag")
        print('  .... to get sample import sheet see j drive')
        print('  .... be sure to use the excel sequencing workbook, not the csv import sheet')
        print('  .... continuing without sample sheet')
        print('')
        print('')
        time.sleep(10)
    
    elif not re.search('.xlsx', sample_sheet):
        print('  ERROR')
        print('  ... sample sheet not formatted correctly; must be excel workbook')
        print('  ... to get sample import sheet see j drive')
        print('  ... be sure to use the excel sequencing workbook, not the csv import sheet')
        print('  ... exiting ...')
        print('')
        print('')
        sys.exit()

        
    
    # get bucket details
    bucket_path = options.bucket_path # this is the path inside the bucket
    if re.search('gs://', bucket_path): 
        bucket_path = re.sub('gs://', '', bucket_path) 
    bucket_path_components = bucket_path.split('/')
    bucket_name = bucket_path_components[0] # this is the bucket name 
    bucket_outerdir_prefix = os.path.join('/'.join(bucket_path_components[1:]), seq_run)
    bucket_prefix = os.path.join('/'.join(bucket_path_components[1:]), '%s/fastq_files' % seq_run) # this is the bucket prefix for the files
    
    # get current directory full path
    current_dir = os.getcwd()

    # check that the current directory has the same name as the seq_run name
    if current_dir.split('/')[-1] != seq_run:
        print('  ERROR: you must name your current directory the same as your seq_run name')
        print('  ....exiting....')
        print('')
        print('')
        sys.exit()
    
    # now run the list of functions
    if sample_sheet != 'not provided':
        sample_sheet_output = read_in_sample_sheet(sample_sheet = sample_sheet, seq_run = seq_run, 
                                                   output_dir = current_dir, terra_output_dir = terra_output_dir)
        sample_sheet_sample_list = sample_sheet_output['sample_list']
        plate_loc_df = sample_sheet_output['df']
    else:
        sample_sheet_sample_list = []
        plate_loc_df = pd.DataFrame()
        
    
    if run_type == 'single':

        fastq_files_dir = create_fastq_files_directory(seq_run = seq_run, 
                                                       current_dir = current_dir)

        samples_without_4_fastq_files = concate_fastq_gz_files_single(current_dir = current_dir, 
                               fastq_files_dir = fastq_files_dir,
                                     sample_sheet_sample_list = sample_sheet_sample_list)

        upload_fastq_files_to_bucket(seq_run = seq_run, 
                                     current_dir = current_dir, 
                                     fastq_files_dir = fastq_files_dir, 
                                     bucket_name = bucket_name, 
                                     bucket_prefix = bucket_prefix)

        create_terra_data_table_single(bucket_path = bucket_path,
                                bucket_name = bucket_name, 
                                bucket_prefix = bucket_prefix, 
                                current_dir = current_dir, 
                                seq_run= seq_run,
                                      plate_loc_df = plate_loc_df,
                                      terra_output_dir = terra_output_dir)                   

#         remove_fastq_gz_directories(current_dir = current_dir)


    elif run_type == 'paired':

        fastq_files_dir = create_fastq_files_directory(seq_run = seq_run, 
                                                       current_dir = current_dir)

        copy_and_rename_fastq_gz_files_paired(seq_run = seq_run, 
                                       current_dir = current_dir,
                                       fastq_files_dir = fastq_files_dir,
                                             sample_sheet_sample_list = sample_sheet_sample_list)

        upload_fastq_files_to_bucket(seq_run = seq_run, 
                                     current_dir = current_dir, 
                                     fastq_files_dir = fastq_files_dir, 
                                     bucket_name = bucket_name, 
                                     bucket_prefix = bucket_prefix)

        create_terra_data_table_paired(bucket_path = bucket_path,
                                bucket_name = bucket_name, 
                                bucket_prefix = bucket_prefix, 
                                current_dir = current_dir, 
                                seq_run= seq_run,
                                      plate_loc_df = plate_loc_df,
                                      terra_output_dir = terra_output_dir)                   

#         remove_fastq_gz_directories(current_dir = current_dir)
        print('')
        print('')
        print('  DONE!')
        print('')
        print('')


