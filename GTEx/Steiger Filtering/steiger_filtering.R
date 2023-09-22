
#load modules.

library(ieugwasr)
library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(dplyr)
library(genetics.binaRies)

#set the location for plink and bcftools.

gwasvcf::set_plink("/mnt/storage/software/apps/Plink-1.90/plink")
gwasvcf::set_bcftools()
ldfile<-"EUR"

#load exposure data

exp_data<-read.delim("amygdala_data.txt", sep=" ")

#Filter by gene(s) of interest

exp_data<-exp_data[exp_data$GENE%in%c("RTEL1"),]

#format exposure data

exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1",
                      other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="FREQA1")

#load outcome VCF data for all glioma data.

all_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/all_glioma.vcf")

#run a query to find proxy SNPs for SNPs that cannot be found in the outcome data.

all_data<-gwasvcf::query_gwas(all_vcf, rsid=unique(exp_data), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "all_index.rsidx")

#format VCF for Steiger filtering.

all_data<-gwasglue::gwasvcf_to_TwoSampleMR(all_data, type="outcome")

#change the name of outcome column to correct phenotype.

all_data$outcome<-"all"

#repeat the same steps for the gbm and non-gbm data.

gbm_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/gbm_glioma.vcf")
gbm_data<-gwasvcf::query_gwas(gbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "gbm_index.rsidx")
gbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(gbm_data, type="outcome")
gbm_data$outcome<-"gbm"

nongbm_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/nongbm_glioma.vcf")
nongbm_data<-gwasvcf::query_gwas(nongbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "nongbm_index.rsidx")
nongbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(nongbm_data, type="outcome")
nongbm_data$outcome<-"nongbm"

#combine all the new data into a new dataframe.

out_data<-rbind(all_data,gbm_data,nongbm_data)

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
