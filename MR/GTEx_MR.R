## Load modules

library(ieugwasr)
library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(dplyr)
library(genetics.binaRies)

## Set the wd for plink and reference genome

gwasvcf::set_plink("/mnt/storage/software/apps/Plink-1.90/plink")
gwasvcf::set_bcftools()
ldfile<-"EUR"

## Load the exposure data and format for MR, then perform clumping

exp_data<-read.delim("nucleus_accumbens_data.txt", sep = " ")
exp_data<-format_data(exp_data,type="exposure", header=T, snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1", other_allele_col="A2",
                      pval_col="P", samplesize_col="N", phenotype_col="GENE", eaf_col="FREQA1")
exp_data<-clump_data(exp_data,clump_r2=0.001)

## Create the rsidx files from vcf files

create_rsidx_index_from_vcf("all_glioma.vcf", "all_index.rsidx")
create_rsidx_index_from_vcf("gbm_glioma.vcf", "gbm_index.rsidx")
create_rsidx_index_from_vcf("nongbm_glioma.vcf", "nongbm_index.rsidx")

## Look for the exact SNPs, or proxies if not available.

all_vcf<-VariantAnnotation::readVcf("all_glioma.vcf")
all_data<-gwasvcf::query_gwas(all_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "all_index.rsidx")
all_data<-gwasglue::gwasvcf_to_TwoSampleMR(all_data, type="outcome")
all_data$outcome<-"all"

## Look for the exact SNPs, or proxies if not available.

gbm_vcf<-VariantAnnotation::readVcf("gbm_glioma.vcf")
gbm_data<-gwasvcf::query_gwas(gbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "gbm_index.rsidx")
gbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(gbm_data, type="outcome")
gbm_data$outcome<-"gbm"

## Look for the exact SNPs, or proxies if not available.

nongbm_vcf<-VariantAnnotation::readVcf("nongbm_glioma.vcf")
nongbm_data<-gwasvcf::query_gwas(nongbm_vcf, rsid=unique(exp_data$SNP), bfile = ldfile, proxies = "yes", tag_r2=0.8, rsidx = "nongbm_index.rsidx")
nongbm_data<-gwasglue::gwasvcf_to_TwoSampleMR(nongbm_data, type="outcome")
nongbm_data$outcome<-"nongbm"

## Combine the datasets back together

out_data<-rbind(all_data,gbm_data,nongbm_data)

## Harmonise the data

dat<-harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)

## Perform MR analysis

res<-mr(dat)
res<-generate_odds_ratios(res)

## Save the table

write.table(res, "na_mr_res.txt", sep="\t", col.names=T, row.names=F, quote=F)
