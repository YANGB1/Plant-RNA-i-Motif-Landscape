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

Example datasets are included in "example" folder.

# Association between iM and environmental variables

The Pearson Correlation Coefficient (PCC): 

``` 
python3 PCC_density_environmental.py

``` 


The phylogenetic generalized least squares (PGLS): 

``` 
library(ape)
library(nlme)

#tree file
nwk = "plant_tree.tre"
tree= read.tree(nwk)



df <- read.csv("iMdesnity_bioclimateVaria.txt")

#data format is as following (in CSV format):
#species,iM_density,BIO5
#Plant1,120.5,28.3
#Plant2,98.2,32.1
#Plant3,110.7,30.5

model <- gls(iM_density~BIO5, correlation = corBrownian(1, phy = tree), data = df)
summary(model)

``` 


# TE-related iM features selection


Spearman Correlation Coefficient (SPCC), Mutual information (MI), u-test, and permutation-based test

``` 
python3 feature_importance_rice_rebuttal.py
``` 
