
#load modules

library(tidyr)

#load outcome data

data<-read.delim("~/Documents/University/PhD Population Health Sciences/QTL/Outcome/GTEx/Amygdala/PWCoCo/RTEL1/am_RTEL1_SNP_outcome.txt",
                 sep = "")

#add sample sizes

data$N<-ifelse(data$DISEASE=="all",30686,ifelse(data$DISEASE=="gbm",24381,24495))

#add case sizes

data$case<-ifelse(data$DISEASE=="all",12496,ifelse(data$DISEASE=="gbm",6191,6305))

#remove the disease phenotype

data<-data[-10]

#filter and order for pwcoco

data<-data[c(1,4,5,6,7,8,9,10,11)]

#separate by sample size into new objects

data_all<-data[data$N==30686,]
data_gbm<-data[data$N==24381,]
data_nongbm<-data[data$N==24495,]

#save the dataframe (all)

write.table(data_all,"~/Documents/University/PhD Population Health Sciences/QTL/Outcome/GTEx/PWCoCo/RTEL1/am_RTEL1_SNP_outcome_all.txt",
            sep=" ", row.names = F, quote=F)

#save the dataframe (gbm)

write.table(data_gbm,"~/Documents/University/PhD Population Health Sciences/QTL/Outcome/GTEx/PWCoCo/RTEL1/am_RTEL1_SNP_outcome_gbm.txt",
            sep=" ", row.names = F, quote=F)

#save the dataframe (non-gbm)

write.table(data_nongbm,"~/Documents/University/PhD Population Health Sciences/QTL/Outcome/GTEx/PWCoCo/RTEL1/am_RTEL1_SNP_outcome_nongbm.txt",
            sep=" ", row.names = F, quote=F)
