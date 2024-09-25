## Load modules

library(ieugwasr)
library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(dplyr)
library(genetics.binaRies)

## Locate PLINK and the reference genome

gwasvcf::set_plink("/mnt/storage/software/apps/Plink-1.90/plink")
gwasvcf::set_bcftools()
ldfile<-"GRCh37_EUR_dedup"

## Load the list of SNPs to be exctracted in the outcome GWAS

exp_data<-read.delim("brainQTL_D2HGDH_SNPs.txt", sep = " ", header=F)

## Find the SNPs or proxy SNPs in the outcome data, and save to a dataframe

all_vcf<-VariantAnnotation::readVcf("all_glioma.vcf")
all_data<-gwasvcf::query_gwas(all_vcf, rsid=unique(exp_data$V1), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "all_index.rsidx")
all_data<-gwasglue::gwasvcf_to_TwoSampleMR(all_data, type="outcome")
all_data$outcome<-"all"

## Find the SNPs or proxy SNPs in the outcome data, and save to a dataframe

gbm_vcf<-VariantAnnotation::readVcf("gbm_glioma.vcf")
gbm_data<-gwasvcf::query_gwas(gbm_vcf, rsid=unique(exp_data$V1), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "gbm_index.rsidx")
gbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(gbm_data, type="outcome")
gbm_data$outcome<-"gbm"

## Find the SNPs or proxy SNPs in the outcome data, and save to a dataframe

nongbm_vcf<-VariantAnnotation::readVcf("nongbm_glioma.vcf")
nongbm_data<-gwasvcf::query_gwas(nongbm_vcf, rsid=unique(exp_data$V1), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "nongbm_index.rsidx")
nongbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(nongbm_data, type="outcome")
nongbm_data$outcome<-"nongbm"

## Combine all three datasets together

out_data<-rbind(all_data,gbm_data,nongbm_data)

## Re-order and rename the columns

out_data<-out_data[c(10,1,2,4,3,8,5,6,7,13)]
colnames(out_data)[1]<-"SNP"
colnames(out_data)[2]<-"CHR"
colnames(out_data)[3]<-"POS"
colnames(out_data)[4]<-"A1"
colnames(out_data)[5]<-"A2"
colnames(out_data)[6]<-"FREQA1"
colnames(out_data)[7]<-"BETA"
colnames(out_data)[8]<-"SE"
colnames(out_data)[9]<-"P"
colnames(out_data)[10]<-"DISEASE"

## Save the data

write.table(out_data, "~/../../work/bp20214/QTL/outcome/BrainQTL/pwcoco/D2HGDH/brainQTL_D2HGDH_SNP_outcome.txt", quote=F, row.names=F, sep = "\t")
