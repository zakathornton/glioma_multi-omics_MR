## Load modules

library(ieugwasr)
library(TwoSampleMR)

#set up the api

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

## Load exposure data

exp_data<-read.delim("hippocampus_data.txt", sep = " ")

## Filter by gene

exp_data<-exp_data[exp_data$GENE%in%c("RTEL1","PHLDB1","FAIM"),]

## Format exposure data

exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1", other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="FREQA1")

## Clump SNPs

exp_data<-clump_data(exp_data,clump_r2=0.001)

## Load outcome data for trait of interest from OpenGWAS API

out_data<-extract_outcome_data(outcomes = c('ebi-a-GCST90002412','ieu-a-1283','ukb-b-20208','ukb-b-110','ukb-b-16702','ukb-b-17219','ukb-b-4601'), snps = exp_data$SNP, proxies = T, rsq=0.8)

## Format outcome data

out_data<-format_data(out_data, type="outcome", header=T, snp_col="SNP", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome", pval_col="pval.outcome", phenotype_col="outcome", eaf_col="eaf.outcome")

## Harmonise strands

dat<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

## Run MR analysis

res<-mr(dat)

## Generate odds ratios and 95% confidence intervals

res<-generate_odds_ratios(res)

## Create a flag for if the p-value passes the Bonferroni-correction.

res$flag<-res$pval<(0.05/res$nsnp)

## Save the results dataframe

write.table(res,"mr_res_all.txt", row.names=F, quote=F, sep = "\t")
