# organize_illumna_fastq

## latest updates:
- Major update as of 11/02/2021.
  - sample sheet is now a required input!
  - script works with NEXSEQ single end reads, NOVASEQ single end reads, and MISEQ paired end reads
  - instead of copying fastq files into the fastq file directory, the script now moves them, inorder to speed up the process and save memory
  - datatable columns include: sample, primer_set, seq_run, download_date, tech_platform, read_type, plate_name, plate_sample_well_and out_dir
  - provides warning when a sample in teh sample sheet does not have the expected number of fastq files.
      - NEXSEQ = 4 fastq files per sample (one for each of 4 lanes)
      - NOVASEQ = 1 fastq file per sample_sheet
      - MiSEQ = 2 fastq files per sample (one forward and one reverse reads)


## Purpose/general jist:
This python script can be used for illumina single end or paired end reads, by specifying ```--read_type``` as "single" or "paired", respectively.
This python script completes the following tasks:
- moves fastq files downloaded from illuminia basespace from their individual directories and into a single directory named ``fastq_files``
- for single end reads (NEXSEQ) the script concatenates all reads into a single .fastq.gz file named according to the sample name
- for paired end reads the script renames the fastq files based on the sample name (e.g. 2100000000_R1_001.fastq.gz)
- moves all fastq files into a single directory called ```fastq_files```
- uploads the fastq file directory to a designated google bucket that can be accessed as input for the illumina-preprocess-assembly workflow on terra
- creates a datatable for input into the illumin-preprocess-assembly workflow on terra and pushes this datatable to the same designated google bucket.

## Requirements:
- The python modules neccessary to run this script are contained in a conda environment; therefore so you must have Anaconda or miniconda installed.

## Preparing your environment:
This only needs to be performed the first time you run the script.
1. Clone the repository to your machine and change directories to the repository on your machine.

``git clone https://github.com/CDPHE/organzie_illumina_fastq.git``

``cd organize_illumna_fastq``

2. Create the conda environment using the ```environment.yml``` file. The environment's name should be ```terra_seq_prep``

```conda env create -f environment.yml```

2. If the environment already exists, then to update the environment
``conda env update -f environment.yml``

3. Check that the environment exists. The name of the enivornment should be ```terra_seq_prep``

```conda env list```

4. Activate the conda environment

```conda activate terra_seq_prep```

## Preparing to run the script:
1. Clone the repository to your machine and change directories ot the repository on your machine:
``git clone https://github.com/CDPHE/organzie_illumina_fastq.git``

2. Make the script executable:
``chmod 755 create_COVMIN_terra_data_table.py``

3. Option 1: add the script to a scripts directory that is listed in your machine's ``$PATH`` variable. In this case you will only need to specify the name of the script each time you run the script.  

4. OPtion 2: specify the full path to the script each time you run the script.

## Running the script:
### Part 1: Preparing your data
1. create a directory that has the same name as the sequencing run name and move into that directory

``mkdir COVSEQ_0000``

``cd COVSEQ_0000``

2. Download the fastq files from basespace.

``bs list projects``

``bs download project -n COVSEQ_0000 -o .``

### Part 2: Running the script
1. activate the conda environment

``conda activate terra_seq_prep``

2. run the script (be sure to still be in the ```COVSEQ_0000``` directory you initally made above). Specify the following flags:
  - ``--seq_run`` : sequecing run name; must match the current output_directory
  - ``--bucket_path`` : google coloud bucket i.e. gs://covid_terra
  - ``--read_type`` : either "paired" or "single". Used paired for COVSEQ and WWT runs; use single for NEXSEQ and NOVASEQ runs
  - ``--sample_sheet`` : path to the excel sample sheet workbook; downloaded from the j drive
  - ``--terra_output_dir`` : (optional) if supplied changes the terra_output directory; defaut is gs://covid_terra/{seq_run}/terra_outputs/
  - ``--remove_bs_dir`` : (optional) default = TRUE, which will remove the bs outer directories which will just be shells (empty directories) by the end of the script anyways


``organize_illumna_fastq.py --seq_run <COVSEQ_0000> --bucket_path gs://<path_to_bucket> --read_type <paired or single>``` --sample_sheet <path to sample sheet>``

## Outputs
1. ```fastq_files``` directory

2. ```COVSEQ_0000_terra_data_table.tsv```
