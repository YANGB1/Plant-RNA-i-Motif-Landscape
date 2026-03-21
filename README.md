# Plant-RNA-i-Motif-Landscape

This repository contains scripts used for the statistical analysis of RNA i-motifs (iMs) in plants and serves as a companion resource to the manuscript “RNA i-motif landscapes in the plant kingdom and their potential functional roles” enabling full reproducibility of the analyses.

# dependency packages

Download and install package pbccs from https://github.com/pacificbiosciences/unanimity/. Program ccs in pbccs package is used to generate highly accurate single-molecule consensus reads (HiFi reads) for PacBio data analysis part.

Download and install package lima from https://github.com/pacificbiosciences/barcoding/. Program lima is used to demultiplex the reads if reads have been multiplexed before sequencing for PacBio data analysis part.

Download and install package pbtk from https://github.com/pacificbiosciences/pbtk/. Program bam2fasta in pbtk package is used to convert the reads from BAM to FASTA format for PacBio data analysis part.

Download and install package Blasr from https://github.com/pacificbiosciences/blasr. Program blasr is used to map the clean reads to the transcriptome reference for PacBio data analysis part.

Download and install package Dorado from https://github.com/nanoporetech/dorado. Program Dorado is used to process Nanopore raw data including base calling and demultiplexing.

Download and install package IsoQuant from https://github.com/ablab/IsoQuant. Program IsoQuant is used to assign clean reads to the transcriptome reference for Nanopore data analysis part.

Download and install package SeqKit from https://github.com/shenwei356/seqkit. Program SeqKit is used to convert reads in FASTQ format into FASTA format for Nanopore data analysis part.

Download and install package LAST from https://github.com/mcfrith/last-rna. Program LAST is used to align long reads with specific RNA reference for Nanopore data analysis part.

Download and install package Perbase from https://github.com/sstadick/perbase. Program Perbase is used to calculate mismatch landscape for Nanopore data analysis part.

Download and install package SAMtools from https://www.htslib.org/. Program SAMtools is used to convert file format for Nanopore data analysis part.

Download and install package CONTRAfold from http://contra.stanford.edu/contrafold/. Program CONTRAfold is used to fold RNA structure by context-free grammar. CONTRAfold is needed for both Nanopore data and PacBio data.

Download and install Python (v3.7 or above) and following packages: 

    DotMap (https://pypi.org/project/dotmap/)
  
    dottree (https://pypi.org/project/dottree/)
  
    Biopython (https://biopython.org/)
  
    Matplotlib (https://matplotlib.org/)
  
    pandas (https://pandas.pydata.org/)
  
    NumPy (https://numpy.org/)
  
    scikit-learn (https://scikit-learn.org/)
  
    Forgi (https://github.com/ViennaRNA/forgi)
  

The dependency packages can be installed by:
``` 
pip3 install -r requirements.txt
``` 
The pipeline can be freely access on GitHub using following command:
``` 
git clone https://github.com/YANGB1/DaVinci-pipeline.git
``` 
The DaVinci pipeline directory can be added to the ‘PATH’ environmental variable or the scripts with full path can be run alternatively.

Alternatively, the python script 'Putative-iM-Searcher.py' can be downloaded directly from Github. The stored directory can be added to the ‘PATH’ environmental variable or the scripts with full path can be run alternatively.

# Step by step workflow for PacBio analysis

Step 1 Obtain the raw PacBio sequencing data. The raw subread data is usually stored in a BAM file. Generate the single-molecule consensus reads (HiFi reads) using program ccs by the following command:
``` 
ccs -j 10 --minPasses=10 raw_subreads.bam HiFi_reads_ccs.bam
``` 
Step 2 Demultiplex the HiFi reads generated in last step by the barcodes using lima by the following command: 
``` 
lima --ccs HiFi_reads_ccs.bam primers.fa lima_out.bam --same --split-bam-named --bam-handles 11 --bam-handles-verbose
``` 
Step 3 Convert the HiFi reads in BAM files from last step to FASTA format by the following command: 
``` 
bam2fasta -u --output barcode1 lima_out.barcode1.consensusreadset.xml
``` 
Step 4 Map the HiFi reads to the transcripts reference by the following command: 
``` 
blasr --hitPolicy leftmost --nproc 8 barcode1.fasta reference.fasta --minMatch 10 -m 5 --out barcode1.m5
``` 
Step 5 Merge the mapped files in different multiplexed primers from same experimental condition from last step using the following linux command: 
``` 
cat barcode1.m5 barcode2.m5 | grep -v '^$' > barcode1_2.m5
``` 
Step 6 Add a prefix for each mapped item of mapping profile from last step by the following command: 
``` 
python3 label_name.py --input_file barcode1.m5 --output_file barcode1_new.m5
```
Step 7 Generate single-molecule resolved RNA structural conformation ensembles for all target RNAs, cluster the ensembles, and select representative conformations. Run the pipeline by the following command: 
``` 
python3 DaVinci_PacBio.py --bit_folder bit_folder --reference_file reference.fasta --number_of_clusters 3 --target_list all --min_reads_thre 20
```

# Step by step workflow for Nanopore analysis
Step 1 Obtain the raw Nanopore sequencing data. The raw sequencing data with ‘.pod5’ suffix is usually stored in a folder ‘pod5_folder’. Basecall the raw data by the following command: 
``` 
dorado basecaller --estimate-poly-a dna_r10.4.1_e8.2_400bps_sup@v4.3.0 pod5_folder --kit-name SQK-NBD114-24 > basecalled_data.bam
``` 
Step 2 Demultiplex the result in BAM format from last step into per-barcode BAMs using the following command: 
``` 
dorado demux --output-dir barcode_divided_folder --no-classify basecalled_data.bam
``` 
Step 3A Convert the per-barcode BAM files (e.g. barcode1.bam) from last step into FASTQ format. The command for single BAM file is as following: 
``` 
samtools fastq barcode1.bam | gzip > barcode1.fastq.gz
``` 
Step 3B The command for BAM files in same folder is as following:
``` 
python3 convert_BAM_to_fastq.py --input_folder divided_barcode_folder --output_folder
``` 
Step 4 Merge the read files in different multiplexed barcodes from same experiment condition from last step using the following linux command:: 
``` 
cat barcode1.fastq.gz barcode2.fastq.gz > barcode1_2.fastq.gz
``` 
Step 5 Assign the reads to transcriptome reference using IsoQuant by the following command: 
``` 
isoquant.py -d nanopore -g gtf_file -r genome_reference_file --fastq barcode1.fastq.gz --stranded none -t 10 -o isoquant_folder --clean_start --fl_data
``` 
Step 6 Split the reads into FASTA files according to their mapped RNA reference by the following command: 
``` 
python3 sort_fasta_by_isoform.py -fi barcode1.fastq.gz -o divided_by_RNA_folder -ri isoquant_OUT.read_assignments.tsv
``` 
Step 7 Generate single-molecule resolved RNA structural conformation ensembles for all target RNAs, cluster the ensembles, and select representative conformations. Run the pipeline by the following command: 
``` 
python3 DaVinci_Nanopore.py --input_folder divided_by_RNA_folder --reference_file reference.fasta --output_folder DaVinci_Nanopore_folder --number_of_clusters 3 --target_list all --min_reads_thre 20
``` 
