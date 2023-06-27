# Instructions for Analysis (GTEx)

### Data is available from https://gtexportal.org/home/datasets - GTEx sQTL v8

## Step 1 - Clean data // data_cleaning.R
Clean the raw dataset to ensure all data is available for MR analysis, then restrict data to a genome-wide level p-value (5*10^-8)

**Data necessary for analysis** - Chromosome, Base pair location, SNP (rsID), Exposure (gene), Effect allele (A1), Other/reference allele (A2), Beta, Standard Error, P Value

Other data - Effect allele frequency (EAF), Sample size (N)

## Step 2 - Retrieve outcome data // retrieve_outcome_data
Find and extract the SNPs identifed in step 1, in the glioma (outcome) data.


