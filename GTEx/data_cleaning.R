
#load modules

library(ieugwasr)
library(TwoSampleMR)

#load exposure data

data<-read.delim("amygdala_raw.txt", sep="\t")

#replace chromosome with just chromosome number

data$chr<-gsub("chr","",data$chr)
data$chr<-as.numeric(data$chr)

#make variant position numeric

data$variant_pos<-as.numeric(data$variant_pos)

#Filter for data needed then re-order

data<-data[c(14,15,19,2,17,16,25,26,24,22)]

#Rename columns

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

#add the sample size.

data$N<-114

#save cleaned raw data

write.table(data,"amygdala_raw.txt", quote = F, row.names = F)

#apply a genome-wide association filter

data<-data[data$P<5E-08,]

#extract a list of SNPs for extraction from outcome data

SNP<-data$SNP

#save the new dataframe

write.table(data,"amygdala_data.txt", quote = F, row.names = F)

#save the list of SNPs for extraction from outcome data

write.table(SNP, "am_SNP.txt", quote = F, row.names = F, col.names = F)
