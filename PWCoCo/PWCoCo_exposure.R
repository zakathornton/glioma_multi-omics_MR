## Load the raw data for the exposure of interest

data<-read.delim("brainQTL_IDH1.txt", sep = " ")

## Add the sample size

data$N<-(1277-data$NMISS)

## Reorder the columns for PWCoCo

data<-data[c(2,4,11,12,7,10,9,16)]

## Change the name of column from MAF to FREQA1

colnames(data)[4]<-"FREQA1"

## Save the formatted data

write.table(data, "brainQTL_IDH1_SNP.txt", quote = F, row.names = F, sep = " ")

## Extract the SNPs from the data to a variable

SNP<-data$SNP

## Save the list of SNPs to be extracted

write.table(SNP, "brainQTL_IDH1_SNPs.txt", quote = F, row.names = F, sep = " ")
