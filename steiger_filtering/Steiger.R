## Load modules

library(ieugwasr)
library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(dplyr)
library(genetics.binaRies)

## Load plink, bcftools and the reference panel

gwasvcf::set_plink("/mnt/storage/software/apps/Plink-1.90/plink")
gwasvcf::set_bcftools()
ldfile<-"EUR"

## Load the exposure data, and subset the genes to perform analysis on.

exp_data<-read.delim("substantia_nigra_data.txt", header=T, sep=" ")
exp_data<-exp_data[exp_data$GENE%in%c("SLC8A1"),]
exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1", other_allele_col="A2", pval_col="P", samplesize_col="N", phenotype_col="GENE$

## 

all_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/all_glioma.vcf")
all_data<-gwasvcf::query_gwas(all_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "all_index.rsidx")
all_data<-gwasglue::gwasvcf_to_TwoSampleMR(all_data, type="outcome")
all_data$outcome<-"all"

##

gbm_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/gbm_glioma.vcf")
gbm_data<-gwasvcf::query_gwas(gbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "gbm_index.rsidx")
gbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(gbm_data, type="outcome")
gbm_data$outcome<-"gbm"

##

nongbm_vcf<-VariantAnnotation::readVcf("~/../../work/bp20214/Glioma/nongbm_glioma.vcf")
nongbm_data<-gwasvcf::query_gwas(nongbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "nongbm_index.rsidx")
nongbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(nongbm_data, type="outcome")
nongbm_data$outcome<-"nongbm"

##

out_data<-rbind(all_data,gbm_data,nongbm_data)

##

dat<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

##

dat$units.exposure<-"SD"
dat$units.outcome<-"log odds"
dat$samplesize.outcome<-ifelse(dat$outcome=="all",30686,ifelse(dat$outcome=="gbm",24381,24009))
dat$ncase.outcome<-ifelse(dat$outcome=="all",12496,ifelse(dat$outcome=="gbm",6191,5819))
dat$ncontrol.outcome<-18190

## Perform steiger filtering

res<-steiger_filtering(dat)

## Save the table

write.table(res, "~/../../work/bp20214/QTL/exposure/GTEx/steiger/substantia_nigra/sn_steig.txt", quote = F, row.names = F, sep = "\t")
