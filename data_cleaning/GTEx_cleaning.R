## Load modules

library(ieugwasr)
library(TwoSampleMR)

## Load exposure data

data<-read.delim("Brain_Substantia_nigra.v8.sgenes.txt", sep="\t")

## Replace chromosome with just chromosome number

data$chr<-gsub("chr","",data$chr)
data$chr<-as.numeric(data$chr)

## Make variant position numeric

data$variant_pos<-as.numeric(data$variant_pos)

## If ref_factor is -1, swap beta and maf

data$maf<-ifelse(data$ref_factor==-1, 1-data$maf, data$maf)
data$slope<-ifelse(data$ref_factor==-1, data$slope*-1, data$slope)

## Filter for data needed then re-order

data<-data[c(14,15,19,2,17,16,25,26,24,22)]

## Appropriately rename columns

colnames(data)[1]<-"CHR"
colnames(data)[2]<-"BP"
colnames(data)[3]<-"SNP"
colnames(data)[4]<-"GENE"
colnames(data)[5]<-"A1"
colnames(data)[6]<-"A2"
colnames(data)[7]<-"BETA"
colnames(data)[8]<-"SE"
colnames(data)[9]<-"P"
colnames(data)[10]<-"FREQA1"

## Manually add the sample size.

data$N<-114

## Save cleaned raw data

write.table(data,"substantia_nigra_raw.txt", quote = F, row.names = F)

## Apply a genome-wide association filter

data<-data[data$P<5E-08,]

## Save the new dataframe

write.table(data,"substantia_nigra_data.txt", quote = F, row.names = F)
