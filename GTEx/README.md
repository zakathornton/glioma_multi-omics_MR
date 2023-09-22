# Instructions for Analysis (GTEx)
Data is available from https://gtexportal.org/home/datasets - GTEx sQTL v8

## Step 1 - Clean exposure data - [Data/Data Cleaning.R]
Clean the raw dataset to ensure all data is available for MR analysis and in correct format, then apply a genome-wide associtation p-value filter (P > 5E-08).

## Step 2 - MR analysis - [MR/MR_analysis.R]
Perform two-sample MR on the exposure and outcome data.

## Step 3 - Manhattan Plot - [Figures/Manhattan_plot.R]
Create a Manhattan plot for the MR results generated from all tissues.

## Step 4 - Retrieve SNPs for PWCoCo (Exposure) - [PWCoCo/retrieve_SNP_exposure]
Extract the cis-SNPs within a 1Mb window of the gene start and end site of the exposure.

## Step 6 - Clean SNPs for PWCoCo (Exposure) - [PWCoCo/data_cleaning_exposure.R]
Clean the exposure SNPs for PWCoCo

## Step 7 - Retrieve SNPs for PWCoCo (Outcome) - [PWCoCo/retrieve_SNP_outcome]


## Step 8 - Clean SNPs for PWCoCo (Outcome) - [PWCoCo/data_cleaning_outcome.R]


## Step 9 - PWCoCo - [PWCoCo/pwcoco]
Perform pair-wise conditional analysis and colocalisation.

## Step 10 - Steiger Filtering - [Steiger Filtering/


## Step 11 - Perform adhoc analysis for risk factors in OpenGWAS [Adhoc Analysis/OpenGWAS.R]


## Step 12 - Forest Plot [Figures/Forest Plot.R]
Create a forest plot for the robust results generated from all tissues.

