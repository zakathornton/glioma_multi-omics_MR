
#load modules

library(ieugwasr)
library(TwoSampleMR)

#load exposure data

exp_data<-read.delim("amygdala_data.txt", sep=" ")

#Filter by gene(s) of interest

exp_data<-exp_data[exp_data$GENE%in%c("RTEL1"),]

#format exposure data

exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1",
                      other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="FREQA1")

#clump SNPs

exp_data<-clump_data(exp_data,clump_r2=0.001)

#load outcome data

out_data<-read.delim("am_SNP_outcome.txt")

#format outcome data

out_data<-format_data(out_data, type="outcome", header=T, snp_col="SNP", beta_col="BETA", se_col="SE",effect_allele_col="A1",
                      other_allele_col="A2", pval_col="P", phenotype_col="DISEASE", eaf_col="FREQA1")

#harmonise data

data<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

#add units of exposure as SD (continuous) and outcome as log odds (binary)

data$units.exposure<-"SD"
data$units.outcome<-"log odds"

#add sample size, ncase and ncontrol, depending on phenotype

data$samplesize.outcome<-ifelse(data$outcome=="all",30686,ifelse(data$outcome=="gbm",24381,24009))
data$ncase.outcome<-ifelse(data$outcome=="all",12496,ifelse(data$outcome=="gbm",6191,5819))
data$ncontrol.outcome<-18190

#perform steiger filtering

res<-steiger_filtering(data)

#save the results data frame

write.table(res, "am_steig.txt", row.names = F, quote = F)
