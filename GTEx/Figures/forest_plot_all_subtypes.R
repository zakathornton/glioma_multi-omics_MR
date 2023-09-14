
#load modules

library(tidyverse)
library(readxl)

#read data

data<-read_excel("~/Documents/University/PhD Population Health Sciences/QTL/Documents/Tables/GTEx Table.xlsx", sheet = 3)

#create an object for each of the tissues in the results

data1<-data[data$Tissue=="Anterior Cingulate Cortex (BA24)",]
data2<-data[data$Tissue=="Cerebellar Hemisphere",]
data3<-data[data$Tissue=="Cerebellum",]
data4<-data[data$Tissue=="Frontal Cortex (BA9)",]
data5<-data[data$Tissue=="Nucleus Accumbens",]
data6<-data[data$Tissue=="Putamen",]

#subset the results to the gene(s) of interest in each tissue

data1<-data1[data1$Gene%in%c("IFT46"),]
data2<-data2[data2$Gene%in%c("IFT46"),]
data3<-data3[data3$Gene%in%c("RTEL1"),]
data4<-data4[data4$Gene%in%c("IFT46"),]
data5<-data5[data5$Gene%in%c("IFT46"),]
data6<-data6[data6$Gene%in%c("IFT46"),]

#combine the results into one dataframe

data<-rbind(data1,data2,data3,data4,data5,data6)

#create a label position for the exposure (gene) and tissue

data$labelpositionExposure<-0.2
data$labelpositionTissue<-0.02

#change the phenotype to a factor and order

data$Phenotype<-factor(data$Phenotype, ordered = is.ordered(data$Phenotype))

#order the data by gene and tissue (alphabetically)

data<-data[order(data$Gene, data$Tissue, data$Phenotype),]

#create table position for each entry

data$TablePosition<-seq(nrow(data),1,by=-1)

#alter gene labels to include the splice junctions

data$Gene<-c(rep("IFT46\nSplice junction 13", 15), rep("RTEL\nSplice junction 32",3))

#create forest plot

plot<-ggplot(dat=data, aes(x=TablePosition,y=`Odds Ratio`,ymin=`Lower 95% CI`,ymax=`Upper 95% CI`, colour=`Phenotype`)) +

      geom_pointrange(aes(colour=factor(`Phenotype`),group=`Gene`),shape=16,size=0.5) +
      geom_errorbar(width=0.5) +
  
      scale_y_log10(breaks=c(1,2,4,8,16,32,64), 
      labels = c("1","2","4","8","16","32","64"),
      expand = c(0,0)) +
  
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.4, ymax=0.005), fill="white", colour="black") +
      
      annotate("text", x=data$TablePosition, y=data$labelpositionExposure, label=data$Gene, size=2.5, fontface="italic") +
      annotate("text", x=data$TablePosition, y=data$labelpositionTissue, label=data$Tissue, size=2.5, fontface="italic") +
        
      annotate("text", x=nrow(data)+1, y=data$labelpositionExposure, label="Exposure", size=4, fontface="bold") +
      annotate("text", x=nrow(data)+1, y=data$labelpositionTissue, label="Tissue", size=4, fontface="bold") +
      
      geom_segment(aes(x=-Inf, y=64, xend=-Inf, yend=4)) +
      geom_segment(aes(x=Inf, y=64,  xend=Inf, yend=4)) +
  
      geom_hline(yintercept=Inf, linetype="solid") +
      geom_hline(yintercept=1, lty=2) +
      geom_vline(xintercept=Inf) +
      geom_vline(xintercept=-Inf) +
  
      coord_flip() +
        
      labs(y="Odds ratio (95% CI) for phenotype risk per standard deviation change in gene splicing",
           x="") +
  
      ggtitle("Mendelian randomisation results by phenotype for the two genes of interest \nin six brain tissues (sQTLs)") +

      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_text(size=13, face="bold"),
            axis.title.x = element_text(size=11),
            panel.background = element_rect(fill="grey95"),
            panel.grid.minor.y = element_line(linewidth = 1))
plot

#save the plot

plot
ggsave("GTEx_FP_supp.png", units = "in", width=12, height=8, dpi=1200)
