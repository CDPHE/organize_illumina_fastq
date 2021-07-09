#! /usr/bin/env python

import os
import glob
import re
import shutil
import pandas as pd 
from datetime import date
import sys

import argparse
from google.cloud import storage
import subprocess


# before running this script....

# 1. create a directory for the sequence run

# 2. download the fastq files from bs into a directory named the same as the seq_run
# for example if you are downloading COVSEQ_0043, the directory must be named 'COVSEQ_0043'

# 3. cd to the directory with the downloaded fastq files (right now this directory will 
# have a series of directories with each sample's fastq.gz files)

# 4. run this script from inside the directory
# note this script takes 3 arguments

## Example input ##
# organiz_fastq_files.py --seq_run <seq_run_name> --bucket_path <gs://path_to_bucket> --run_type <paired or single>
# organize_fastq_files.py --seq_run COVSEQ_0050rr --bucket_path gs://covid_terra --run_type paired



#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--seq_run", help="sequencing run name; must match directory")
    parser.add_argument("--bucket_path", help="google cloud bucket path; gs://covid_terra")
    parser.add_argument("--run_type",  help="either 'paired' or 'single'")
    options = parser.parse_args(args)
    return options


def create_fastq_files_directory(seq_run, current_dir):
    print('******* creating fastq_files directory *******')
    
    fastq_files_dir = os.path.join(current_dir, 'fastq_files')

    if os.path.exists(fastq_files_dir):
        shutil.rmtree(fastq_files_dir)
    
    os.mkdir(fastq_files_dir)
                             
    return fastq_files_dir



def concate_fastq_gz_files_single(current_dir, fastq_files_dir):
    print('******* concatenating 4 fastq files per sample id *******')
    
    os.chdir(current_dir)
    
    # get a list of the sample ids
    sample_list = []
    for file in glob.glob('*_L001*/*.fastq.gz'):
        sample_name = file.split('_L00')[0]
        if sample_name not in sample_list:
            sample_list.append(sample_name)
    
    print('there are %d sample ids in %s' % (len(sample_list), seq_run))
    print('')
    
    # now concatenate using subprocesss command
    n = 0
    for sample in sample_list:
        n = n + 1
        print('............%d: concatenating files for sample %s' % (n, sample))
        
        # create empty list to hold fastq files assoicated with sample
        file_list = []
        for file in glob.glob('*_L00*/*.fastq.gz'):
        
            sample_name = file.split('_L00')[0]
            if sample_name == sample:
                file_list.append(file)

        combined_fastq_file_name = '%s.fastq.gz' % sample
        out_path = os.path.join(fastq_files_dir, combined_fastq_file_name)
        
        print('concatenating the following files into single file named: %s' % (combined_fastq_file_name))
        i = 1
        for file_name in file_list:
            fastq_file_name = file_name.split('/')[1]
            print('%d: %s' % (i,fastq_file_name))
            i = i + 1
        file_list_string = ' '.join(file_list)
        
        shell_command = 'cat %s > %s' % (file_list_string, out_path)

        subprocess.run(args= shell_command, shell= True, check = True)
        print('')

def copy_and_rename_fastq_gz_files_paired(seq_run, current_dir, fastq_files_dir):
    print('')
    print('******* copying and renaming fastq.gz files *******')
    
    os.chdir(current_dir)

    n = 0
    for file in glob.glob('*_L001_ds.*/*.fastq.gz'):
        n = n+1
        sample_name = file.split('_L001')[0] # the base directory name


        sequencing_info = re.findall('_R\d_\d\d\d', file)


        new_file_name = '%s%s.fastq.gz' % (sample_name, sequencing_info[0])
        old_file_name = file.split('/')[1]

        # check that new name makes sense:
        print('%d: %s is renamed to %s' % (n, old_file_name, new_file_name))

        # copy file into new directory
        shutil.copy(file, fastq_files_dir)

        # rename the file
        old= os.path.join(fastq_files_dir, old_file_name)
        new = os.path.join(fastq_files_dir, new_file_name)
        os.rename(old, new)

def upload_fastq_files_to_bucket(seq_run, current_dir, fastq_files_dir, bucket_name, bucket_prefix):
    print('******* uploading fastq_file_directory to google bucket *******')
    
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
        print('{}: uploaded {} to {} bucket'.format(n, file, os.path.join(bucket_name, bucket_prefix)))
    

          
def create_terra_data_table_single(bucket_path, bucket_name, bucket_prefix, current_dir, seq_run ):
    print('******* creating datatable for terra *******')

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
        

    df['sample_name'] = sample_name_list
    df['fastq'] = fastq_path_list 

    df = df.rename(columns = {'sample_name' : 'entity:sample%s_id' % seq_run_mod})
          
    outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
    df.to_csv(outfile, index = False, sep = '\t')


    os.chdir(current_dir)
    # upload terra datatable to bucket
    print('******* uploading terra datatable to bucket *******')
    client_storage = storage.Client()
    bucket = client_storage.get_bucket(bucket_name)
    
#     localFile = os.path.join(fastq_files_dir, file)
    bucket_path = os.path.join(bucket_outerdir_prefix, '%s_terra_data_table.tsv' % seq_run)
    blob = bucket.blob(bucket_path)
    blob.upload_from_filename(outfile)
    print('uploaded {} to {} bucket'.format( outfile, os.path.join(bucket_name, bucket_outerdir_prefix)))
    
def create_terra_data_table_paired(bucket_path, bucket_name, bucket_prefix, current_dir, seq_run ):
    print('******* creating datatable for terra *******')

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

    R1_df['sample_name'] = R1_sample_name_list
    R1_df['fastq_1'] = R1_fastq_list
    R1_df = R1_df.set_index('sample_name')      
    
    
    R2_df['sample_name'] = R2_sample_name_list
    R2_df['fastq_2'] = R2_fastq_list
    R2_df = R2_df.set_index('sample_name')
          
    df = R1_df.join(R2_df, how='left')
    df = df.reset_index()
    df = df.rename(columns = {'sample_name' : 'entity:sample%s_id' % seq_run_mod})
          
    outfile = os.path.join(current_dir, '%s_terra_data_table.tsv' % seq_run)
    df.to_csv(outfile, index = False, sep = '\t')
    
    os.chdir(current_dir)
    
    # upload terra datatable to bucket
    print('******* uploading terra datatable to bucket *******')
    client_storage = storage.Client()
    bucket = client_storage.get_bucket(bucket_name)
    
#     localFile = os.path.join(fastq_files_dir, file)
    bucket_path = os.path.join(bucket_outerdir_prefix, '%s_terra_data_table.tsv' % seq_run)
    blob = bucket.blob(bucket_path)
    blob.upload_from_filename(outfile)
    print('uploaded {} to {} bucket'.format( outfile, os.path.join(bucket_name, bucket_outerdir_prefix)))    
    
def remove_fastq_gz_directories(current_dir):
    print('******* removing fastq.gz directories *******')
    os.chdir(current_dir)

    for directory in os.listdir():
        if re.search('_ds.', directory):
            shutil.rmtree(directory)
            
            
            
if __name__ == '__main__':
    
    options = getOptions()
    
    # get teh seq_run name
    seq_run = options.seq_run # this is the sequence run name
    run_type = options.run_type # this is the run type either single or paired
    
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
        print('you must name your directory the same as your seq_run name')
        go_ahead = 'NO'
    else:
        go_ahead = 'YES'
    
    # now run the list of functions
    if go_ahead == 'YES':
        
        if run_type == 'single':

            print('Going to organize fastq.gz file for sequence run: %s' % seq_run)
            print('run type = single')

            fastq_files_dir = create_fastq_files_directory(seq_run = seq_run, 
                                                           current_dir = current_dir)

            concate_fastq_gz_files_single(current_dir = current_dir, 
                                   fastq_files_dir = fastq_files_dir)

            upload_fastq_files_to_bucket(seq_run = seq_run, 
                                         current_dir = current_dir, 
                                         fastq_files_dir = fastq_files_dir, 
                                         bucket_name = bucket_name, 
                                         bucket_prefix = bucket_prefix)

            create_terra_data_table_single(bucket_path = bucket_path,
                                    bucket_name = bucket_name, 
                                    bucket_prefix = bucket_prefix, 
                                    current_dir = current_dir, 
                                    seq_run= seq_run)                   

            remove_fastq_gz_directories(current_dir = current_dir)

            print('DONE!')
            
            
        elif run_type == 'paired':
            
            print('Going to organize fastq.gz file for sequence run: %s' % seq_run)
            print('run type = paired')

            fastq_files_dir = create_fastq_files_directory(seq_run = seq_run, 
                                                           current_dir = current_dir)

            copy_and_rename_fastq_gz_files_paired(seq_run = seq_run, 
                                           current_dir = current_dir,
                                           fastq_files_dir = fastq_files_dir)

            upload_fastq_files_to_bucket(seq_run = seq_run, 
                                         current_dir = current_dir, 
                                         fastq_files_dir = fastq_files_dir, 
                                         bucket_name = bucket_name, 
                                         bucket_prefix = bucket_prefix)

            create_terra_data_table_paired(bucket_path = bucket_path,
                                    bucket_name = bucket_name, 
                                    bucket_prefix = bucket_prefix, 
                                    current_dir = current_dir, 
                                    seq_run= seq_run)                   

            remove_fastq_gz_directories(current_dir = current_dir)
            print('DONE!')


