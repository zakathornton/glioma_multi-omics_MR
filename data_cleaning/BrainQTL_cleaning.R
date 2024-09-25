
## Load modules

library(tidyr)
library(biomaRt)
library(dplyr)

## Load exposure data

data<-read.csv("brainQTL_raw_filter.txt", sep = " ")

## Add sample sizes

data$N<-(1277-data$NMISS)

## Set up BioMaRt

ensembl<-useMart("ensembl")
ensembl<-useDataset("hsapiens_gene_ensembl", mart = ensembl)

## Get external gene name from BioMaRt

biomart<-getBM(attributes = c("ensembl_gene_id","external_gene_name"),
      filters = "ensembl_gene_id",
      values = data$GENE,
      mart = ensembl)
biomart<-distinct(biomart)

## Rename the columns

colnames(biomart)[1]<-"ENSG"
colnames(biomart)[2]<-"GENE"
colnames(data)[13]<-"ENSG"

## Merge the two datasets

data<-merge(data, biomart, by="ENSG")

## Reorder the dataset

data<-data[c(2,4,3,17,5,12,13,8,11,10,16)]

## Change column names

colnames(data)[2]<-"POS"
colnames(data)[7]<-"FREQA1"

## Save the table

write.table(data, "brainQTL_data.txt", quote = F, row.names = F, sep = " ")
