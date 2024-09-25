## Load modules

library(tidyr)

## Load exposure data

data<-read.delim("basal_ganglia_raw.txt", sep = " ")

## Separate chr, bp, rsid and allele

data<-separate(data=data, col = SNP, into = c("Chr","BP","rsID","SNPAllele"), sep = ":")

## Seperate the two alleles 

data<-separate(data=data, col=SNPAlleles, into = c("SNP1", "SNP2"), sep = "/")

## Set the effect and other alleles by comparing them to AlleleAssessed

data$A1<-data$SNPEffectAllele
data$A2<-ifelse(data$A1==data$SNP1, data$SNP2, data$SNP1)

## Filter the data set to information needed for MR analysis

data<-data[c(6,7,8,5,25,26,19,20,16,15,17)]

## Rename the columns appropriately

colnames(data)[1]<-"CHR"
colnames(data)[2]<-"BP"
colnames(data)[3]<-"SNP"
colnames(data)[4]<-"GENE"
colnames(data)[5]<-"A1"
colnames(data)[6]<-"A2"
colnames(data)[7]<-"BETA"
colnames(data)[8]<-"SE"
colnames(data)[9]<-"P"
colnames(data)[10]<-"EAF"
colnames(data)[11]<-"N"

## Save the new dataframe

write.table(data,"~/Documents/University/PhD Population Health Sciences/QTL/Exposure/MetaBrain/Data/Basal Ganglia/basal_ganglia_raw.txt", quote=F, row.names = F)

## Apply a genome-wide association filter

data<-data[data$P<5E-08,]

## Save the new dataframe

write.table(data,"basal_ganglia_data.txt", quote=F, row.names = F)
