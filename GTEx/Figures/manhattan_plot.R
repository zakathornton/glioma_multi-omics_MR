
#load modules

library(biomaRt)
library(ggrepel)
library(dplyr)
library(qqman)

#load compiled MR results

data<-read.delim("comp_mr_res.txt")

#make p-value numeric

data$pval<-as.numeric(data$pval)

#create variable for the results which meet bonferroni correction

signif<-data[data$pval<1.11E-05,]

#create a vector of all the unique genes in the dataset

hgnc<-unique(data$exposure)

#load a biomart for finding missing data

grch37=useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               host = "grch37.ensembl.org",
               path = "/biomart/martservice",
               dataset = "hsapiens_gene_ensembl")

#retrieve the chromosome number, and start and end position of each gene

table<-getBM(attributes = c("hgnc_symbol","start_position", "end_position","chromosome_name"),
             filters = "hgnc_symbol",
             values = hgnc,
             mart = grch37)

#change the name of of the variable to exposure (to merge with the main dataset)

colnames(table)[1]<-"exposure"

#merge the new data with the MR results

data<-merge(data,table,by="exposure")

#make chromosome and start position numeric class, and remove any SNPs with missing info

data$start_position<-as.numeric(data$start_position)
data$chromosome_name<-as.numeric(data$chromosome_name)
data<-data[!is.na(data$chromosome_name),]

#remove the major histocompatability complex in chromosome 6

data<-data[!(data$chromosome_name==6 & data$start_position > 29000000 & data$end_position < 33000000),]

#create Manhattan plot and save

png("~/Documents/University/PhD Population Health Sciences/QTL/Exposure/GTEx/Figures/GTEx_manhattan_plot.png",
    width = 2000, height = 1000, pointsize = 18)
manhattan(data, chr = "chromosome_name", bp = "start_position", snp = "exposure", p = "pval", logp = T, suggestiveline = -log10(1.11E-05),
    genomewideline = F, annotateTop = TRUE, main="Manhattan Plot of Mendelian Randomisation Results in 13 Brain Tissues (sQTL)",
    cex=0.6, cex.axis=0.9, ylim=c(0,50), col = c("blue","orange"))
dev.off()
