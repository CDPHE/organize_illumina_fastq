# seq_assembly_prep_for_terra

## Purpose:
This python script can be used for illumina single end or paired end reads, by specifying ```--run_type``` as "single" or "paired", respectively.
This python script completes the following tasks:
- removes fastq files downloaded from illuminia basespace from their individual directories
- for single end reads, the script concatenates all reads into a single .fastq.gz file named according to the sample name.
- for paired end reads the script renames the fastq files based on the sample name (e.g. 2100000000_R1_001.fastq.gz)
- moves all fastq files into a single directory called ```fastq_files```
- uploads the fastq file directory to a designated google bucket that can be accessed as input for the illumina-preprocess-assembly workflow on terra
- creates a datatable for input into the illumin-preprocess-assembly workflow on terra and pushes this datatable to the same designated google bucket.

## Requirements:
- The python modules neccessary to run this script are contained in a conda environment; therefore so you must have Anaconda or miniconda installed.

## Preparing your environment:
This only needs to be performed the first time you run the script.
1. Clone the repository to your machine and change directories to the repository on your machine.

```git clone https://github.com/CDPHE/seq_assembly_prep_for_terra.git```

``` cd seq_assembly_prep_for_terra```

2. Create the conda environment using the ```environment.yml``` file. The environment's name should be ```terra_seq_prep```

```conda env create -f environment.yml```

3. Check that the environment exists. The name of the enivornment should be ```terra_seq_prep``

```conda env list```

4. Activate the conda environment

```conda activate terra_seq_prep```

## Running the script:
### Part 1: Preparing your data
1. create a directory that has the same name as the sequencing run name and move into that directory

```mkdir COVSEQ_0000```

```cd COVSEQ_0000```

2. Download the fastq files from basespace.

```bs list projects```

```bs download project -n COVSEQ_0000 -o .```

### Part 2: Running the script
1. activate the conda environment

```conda activate terra_seq_prep```

2. run the script (be sure to still be in the ```COVSEQ_0000``` directory you initally made above). The order of the arguments should be the path to the python script, followed by the name of the sequencing run (which should be the same name as the directory you created above) and finally the path to the google storage bucket where you want the fastq files directory to be stored.

```<path_to_seq_assembly_for_terra>/prep_illumina_fastq_files_for_terra.py --seq_run <COVSEQ_0000> --bucket_path gs://<path_to_bucket> --run_type <paired or single>```

## Outputs
1. ```fastq_files``` directory

2. ```COVSEQ_0000_terra_data_table.tsv```
