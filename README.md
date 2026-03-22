# RNA i-Motif Landscape in Plant Kingdom

This repository contains scripts used for the statistical analysis of RNA i-motifs (iMs) in plants and serves as a companion resource to the manuscript “RNA i-motif landscapes in the plant kingdom and their potential functional roles” enabling full reproducibility of the analyses.

# Key dependency packages

Download and install Python and following packages: 

    scipy 
  
    pandas 
  
    NumPy 
  
    scikit-learn 

Download and install R and following packages: 

    ape 
  
    nlme 
  

# Association between iM and environmental variables: Pearson Correlation Coefficient (PCC) and Phylogenetic Generalized Least Squares (PGLS)

The phylogenetic generalized least squares (PGLS) calculation: 

``` 
library(ape)
library(nlme)

#tree file
nwk = "plant_tree.tre"
tree= read.tree(nwk)



df <- read.csv("iMdesnity_enrichment.txt")

#data format is as following (in CSV format):
#species,iM_density,BIO5
#Plant1,120.5,28.3
#Plant2,98.2,32.1
#Plant3,110.7,30.5

model <- gls(iM_density~BIO5, correlation = corBrownian(1, phy = tree), data = df)
summary(model)

``` 


# TE-related iM features selection: Spearman Correlation Coefficient (SPCC), Mutual information (MI), u-test, and permutation-based test

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
