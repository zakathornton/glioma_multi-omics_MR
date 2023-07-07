
#load modules

library(tidyverse)
library(readxl)

#read data

data<-read_excel("~/QTL/Documents/Tables/GTEx Table.xlsx", sheet = 4)

#filter to only include data which passes colocalisation and steiger filtering

data<-data[data$`H4 (%)`>80,]
data<-data[data$`Steiger Flag`!=FALSE,]

#create a label position for the exposure (gene) and tissue

data$labelpositionExposure=rep(0.015,times=nrow(data))
data$labelpositionTissue=rep(0.002,times=nrow(data))

#order the data by gene and tissue (alphabetically)

data<-data[order(data$`Odds Ratio`, decreasing = TRUE),]

#create table position for each entry

data$TablePosition<-seq(nrow(data),1,by=-1)

#create forest plot

plot<-ggplot(dat=data, aes(x=TablePosition,y=`Odds Ratio`,ymin=`Lower 95% CI`,ymax=`Upper 95% CI`, colour=`Phenotype`)) +

      geom_pointrange(aes(colour=factor(`Phenotype`),group=`Gene`),shape=16,size=0.5) +
  
      scale_y_log10(breaks=c(0.0625,0.125,0.25,0.5,1,2), 
      labels = c("0.0625","0.125","0.25","0.50","1","2"),
      expand = c(0,0)) +
  
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.00048, ymax=0.03), fill="white", colour="black") +
      
      annotate("text", x=data$TablePosition, y=data$labelpositionExposure, label=data$Gene, size=3, fontface="italic") +
       annotate("text", x=data$TablePosition, y=data$labelpositionTissue, label=data$Tissue, size=3, fontface="italic") +
        
      annotate("text", x=nrow(data)+1, y=data$labelpositionExposure, label="Exposure", size=4, fontface="bold") +
      annotate("text", x=nrow(data)+1, y=data$labelpositionTissue, label="Tissue", size=4, fontface="bold") +
      
      geom_segment(aes(x=-Inf, xend=-Inf, yend=4)) +
      geom_segment(aes(x=Inf, xend=Inf, yend=4)) +
  
      geom_hline(yintercept=Inf, linetype="solid") +
      geom_hline(yintercept=1, lty=2) +
      geom_vline(xintercept=Inf) +
      geom_vline(xintercept=-Inf) +
  
      scale_color_manual(values=c("#27A829","#E03132")) + 
  
      coord_flip() +
        
      labs(y="Odds ratio (95% CI) for phenotype risk per standard deviation change in gene splicing",
           x="") +
  
      ggtitle("Mendelian randomisation results by phenotype for the two genes of interest in six brain tissues (sQTLs)") +

      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_text(size=10, face="bold"),
            axis.title.x = element_text(size=10),
            panel.background = element_rect(fill="grey95"),
            panel.grid.minor.y = element_line(size = 1))
plot

#save the plot

plot
ggsave("~/QTL/Exposure/GTEx/Figures/GTEx_forest_plot.png", units = "in", width=12, height=8, dpi=1200)
