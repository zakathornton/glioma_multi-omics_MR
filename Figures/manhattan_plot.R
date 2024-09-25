## Load modules

library(biomaRt)
library(ggrepel)
library(dplyr)
library(qqman)

## Load MR results

data<-read.delim("comp_mr_res.txt", sep="\t")

## Set p-value as a numeric class.

data$pval<-as.numeric(data$pval)

## Create variable for the results which meet bonferroni correction

signif<-data[data$pval<5E-06,]

## Create a vector of all the unique genes in the dataset

hgnc<-unique(data$exposure)

## Load a biomart for finding missing data

grch37<-useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               host = "grch37.ensembl.org",
               path = "/biomart/martservice",
               dataset = "hsapiens_gene_ensembl")

## Retrieve the chromosome number, and start and end position of each gene

table<-getBM(attributes = c("hgnc_symbol","start_position", "end_position","chromosome_name"),
             filters = "hgnc_symbol",
             values = hgnc,
             mart = grch37)

## Change the name of of the variable to exposure (to merge with the main dataset)

colnames(table)[1]<-"exposure"

## Merge the new data with the MR results

data<-merge(data,table,by="exposure")

## Make chromosome and start position numeric class, and remove any SNPs with missing info

data$start_position<-as.numeric(data$start_position)
data$chromosome_name<-as.numeric(data$chromosome_name)
data<-data[!is.na(data$chromosome_name),]

## Remove the major histocompatability complex in chromosome 6

data<-data[!(data$chromosome_name==6 & data$start_position > 29000000 & data$end_position < 33000000),]

## Create Manhattan plot and save

png("MB_manhattan_plot.png", width = 2000, height = 1000, pointsize = 18)
manhattan(data, chr = "chromosome_name", bp = "start_position", snp = "exposure", p = "pval", logp = T, suggestiveline = -log10(5E-06), genomewideline = FALSE, annotateTop = TRUE, main="Manhattan plot of Mendelian randomisation results in five brain tissues (eQTL)", cex=0.6, cex.axis=0.9, ylim=c(0,50), col = c("blue","orange"))
dev.off()
