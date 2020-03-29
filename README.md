This repository contains data and scripts used in the paper

[***Dissecting vancomycin intermediate resistance in Staphylococcus aureus using genome-wide association 
Md Tauqeer Alam; Robert A. Petit III; Emily K. Crispelll III; Timothy A. Thornton; Karen N. Conneely; Yunxuan Jiang; Sarah W. Satola; Timothy D. Read
Genome Biology and Evolution 2014;
doi: 10.1093/gbe/evu092***](http://gbe.oxfordjournals.org/content/early/2014/04/30/gbe.evu092.abstract?keytype=ref&ijkey=zDyKdOE4oaYODaI)

### WEBSITE
  We have published an [interactive webserver](https://tread.shinyapps.io/VISA-shiny/) to help exploration of the data.  The site was built using the RStudio Shiny package. 

### SCRIPTS

###### visa-table-manipulation.R
    Script to automate rare mutation counting.

###### visa_gwas.R
    Script to automate ROADTRIPS and QROADTRIPS analysis using E-Test values.
    
###### visa_gwas_bmd.R
    Script to automate ROADTRIPS and QROADTRIPS analysis using BMD values.
    
###### visa_gwas_pap-auc.R
    Script to automate ROADTRIPS and QROADTRIPS analysis using PAP-AUC values.
    
###### visa_gwas_functions.R
    Script containing functions for process data to be used with ROADTRIPS and 
    QROADTRIPS.

##### Interactive_VISA_model
    This directory contains scripts, data and tests for the shiny interactive website.  

### DATA

###### superAlignProt124.phy
    Alignment of core genes 

###### superAlignDNA124.phy
    Alignment of core proteins 

###### superAlignDNA124.fna
    Alignment of core proteins

###### VISA_feature_table.xls
    Created by subsetting the SNP matrix by the 16 candidate genes.  Also added 
    in the vector of SNP frequency for in the Staph genome database.

###### VISA-feature-table-v2-PAPAUC
    Transpose VISA_feature_table.xls and manually removed all SNPs with a count 
    > 81 (ie common SNPs).  Add PAP-AUC data. Used to make prelim version of 
    candidate SNP frequency table.

###### VISA-feature-table-v2-etest
    From VISA-feature-table-v2-PAPAUC but used for Etest counts.

###### VISA-feature-table-v3.csv
    Adapted from VISA-feature-table-v2-etest.  Stripped out unnecessary fields.  This is the input file for the Shiny server and the random forest model. (Note - this was moved into the interactive_VISA_model directory in order to publish the webserver.)  
    
###### genotype.txt
    A matrix of SNPs for the total 75 samples.

###### phenotype.txt
    Phenotype values from E-Test for each sample.

###### phenotype_bmd.txt
    Phenotype values from BMD for each sample.

###### phenotype_pap-auc.txt
    Phenotype values from PAP-AUC for each sample.