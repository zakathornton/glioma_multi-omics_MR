
#load modules

library(ieugwasr)
library(TwoSampleMR)

#load exposure data

exp_data<-read.delim("amygdala_data.txt", sep = " ")

#format exposure data

exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1",
                      other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="FREQA1")

#clump SNPs

exp_data<-clump_data(exp_data,clump_r2=0.001)

#load outcome data

out_data<-read.delim("am_SNP_outcome.txt")

#format outcome data

out_data=format_data(out_data, type="outcome", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1",
                     other_allele_col="A2", pval_col="P", phenotype_col="DISEASE", eaf_col="FREQA1")

#harmonise data into new dataframe

dat<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

#run MR

res<-mr(dat)

#generate odds ratios and 95% confidence intervals

res<-generate_odds_ratios(res)

#save results

write.table(res, "am_mr_res.txt", sep = "\t" row.names=F, quote=F)
