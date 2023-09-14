
#load modules

library(tidyr)

#load data

data<-read.delim("am_RTEL1_SNP.txt", sep=",")

#add sample size

data$N<-175

#filter and order for pwcoco

data<-data[c(1,5,4,10,8,9,7,11)]

#rename columns

colnames(data)[1]<-"SNP"
colnames(data)[2]<-"A1"
colnames(data)[3]<-"A2"
colnames(data)[4]<-"EAF"
colnames(data)[5]<-"BETA"
colnames(data)[6]<-"SE"
colnames(data)[7]<-"P"

#save the new dataset

write.table(data, "am_RTEL1_SNP.txt", quote = F, row.names = F)
