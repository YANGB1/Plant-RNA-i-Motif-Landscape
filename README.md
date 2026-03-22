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
  

# Association between iM and environmental variables

The Pearson Correlation Coefficient (PCC) calculation: 


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
python3 DaVinci_Nanopore.py --input_folder divided_by_RNA_folder --reference_file reference.fasta --output_folder DaVinci_Nanopore_folder --number_of_clusters 3 --target_list all --min_reads_thre 20
``` 
