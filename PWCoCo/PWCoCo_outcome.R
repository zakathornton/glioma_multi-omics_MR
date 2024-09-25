## Load modules

library(tidyr)

## Load outcome data

data<-read.delim("brainQTL_D2HGDH_SNP_outcome.txt", sep = "")

## Add sample sizes

data$N<-ifelse(data$DISEASE=="all",30686,ifelse(data$DISEASE=="gbm",24381,24495))

## Add case sizes

data$case<-ifelse(data$DISEASE=="all",12496,ifelse(data$DISEASE=="gbm",6191,6305))

## Remove the disease phenotype column

data<-data[-10]

## Filter and order for PWCoCo

data<-data[c(1,4,5,6,7,8,9,10,11)]

## Separate by sample size into new objects

data_all<-data[data$N==30686,]
data_gbm<-data[data$N==24381,]
data_nongbm<-data[data$N==24495,]

## Save the dataframe (all)

write.table(data_all,"brainQTL_D2HGDH_SNP_outcome_all.txt", sep=" ", row.names = F, quote=F)

## Save the dataframe (gbm)

write.table(data_all,"brainQTL_D2HGDH_SNP_outcome_gbm.txt", sep=" ", row.names = F, quote=F)

## Save the dataframe (non-gbm)

write.table(data_all,"brainQTL_D2HGDH_SNP_outcome_nongbm.txt", sep=" ", row.names = F, quote=F)
