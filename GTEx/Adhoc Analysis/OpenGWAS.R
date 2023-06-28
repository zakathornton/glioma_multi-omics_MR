
#load modules

library(ieugwasr)
library(TwoSampleMR)

#load exposure data

exp_data<-read.delim("~/Documents/University/PhD Population Health Sciences/QTL/Exposure/GTEx/Data/Cerebellum/cerebellum_data.txt", sep = " ")

#filter by gene

exp_data<-exp_data[exp_data$GENE%in%c("IFT46","RTEL1"),]

#format exposure data

exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1",
                      other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="EAF")


#clump SNPs

exp_data<-clump_data(exp_data,clump_r2=0.001)

#load outcome data

out_data<-extract_outcome_data(outcomes = 'ebi-a-GCST90001390', snps = exp_data$SNP)

#format outcome data

out_data<-format_data(out_data, type="outcome", header=T, snp_col="SNP", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome", pval_col="pval.outcome", phenotype_col="outcome", eaf_col="eaf.outcome")

#harmonise strands

dat<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

#run MR analysis

res<-mr(dat)

#generate odds ratios and 95% confidence intervals

res<-generate_odds_ratios(res)

#create a flag for if the p-value passes the Bonferroni-correction.

res$flag<-res$pval<(0.05/res$nsnp)

#save the results dataframe

write.table(res,"~/Documents/University/PhD Population Health Sciences/QTL/Exposure/GTEx/Adhoc Analysis/mr_res_LBD.txt", row.names=F, quote=F)
