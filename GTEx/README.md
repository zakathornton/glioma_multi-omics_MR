# Instructions for Analysis (GTEx)

### Data is available from https://gtexportal.org/home/datasets - GTEx sQTL v8

## Step 1 - Clean exposure data \\ Data/data_cleaning_exposure.R
Clean the raw dataset to ensure all data is available for MR analysis, then apply a genome-wide associtation p-value filter (5*10^-8)

## Step 2 - Clean outcome data \\ Data/data_cleaning_outcome
Extract the SNPs identifed in step 1 in the glioma (outcome) data.

## Step 3 - MR analysis \\ MR/MR_analysis.R
Perform two-sample MR on the exposure and outcome data.

## Step 4 - Manhattan Plot \\ Figures/Manhattan_plot.R
Create a Manhattan plot for the MR results generated from all tissues.

## Step 5 - Retrieve SNPs for PWCoCo (Exposure) \\ PWCoCo/retrieve_SNP_exposure
Extract the cis-SNPs within a 1Mb window of the gene start and end site of the exposure.

## Step 6 - Clean SNPs for PWCoCo (Exposure) \\ PWCoCo/data_cleaning_exposure.R
Clean the SNPs for PWCoCo

## Step 7 - Retrieve SNPs for PWCoCo (Outcome) \\ PWCoCo/retrieve_SNP_outcome

## Step 8 -

## Step 9 - 

## Step 10 - 

## Step 11 - 

## Step 12 - 


