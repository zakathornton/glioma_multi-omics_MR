
#load modules

library(tidyverse)
library(readxl)

#read data

data<-read_excel("GTEx Table.xlsx", sheet = 6)

#filter to only include data which passes colocalisation and steiger filtering

data<-data[data$`H4 (%)`=>80,]
data<-data[data$`Steiger Flag`!=FALSE,]

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

data$Gene<-c(rep("IFT46\nSplice junction 13", 5), "RTEL\nSplice junction 32")

#create forest plot

plot<-ggplot(dat=data, aes(x=TablePosition,y=`Odds Ratio`,ymin=`Lower 95% CI`,ymax=`Upper 95% CI`, colour=`Phenotype`)) +

      geom_pointrange(aes(colour=factor(`Phenotype`),group=`Gene`),shape=16,size=0.5) +
      geom_errorbar(width=0.1) +
  
      scale_y_log10(breaks=c(1,2,4,8,16,32), 
      labels = c("1","2","4","8","16","32"),
      expand = c(0,0)) +
  
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.5, ymax=0.005), fill="white", colour="black") +
      
      annotate("text", x=data$TablePosition, y=data$labelpositionExposure, label=data$Gene, size=3, fontface="italic") +
      annotate("text", x=data$TablePosition, y=data$labelpositionTissue, label=data$Tissue, size=3, fontface="italic") +
        
      annotate("text", x=nrow(data)+1, y=data$labelpositionExposure, label="Exposure", size=4, fontface="bold") +
      annotate("text", x=nrow(data)+1, y=data$labelpositionTissue, label="Tissue", size=4, fontface="bold") +
      
      geom_segment(aes(x=-Inf, y=32, xend=-Inf, yend=4)) +
      geom_segment(aes(x=Inf, y=32, xend=Inf, yend=4)) +
  
      geom_hline(yintercept=Inf, linetype="solid") +
      geom_hline(yintercept=1, lty=2) +
      geom_vline(xintercept=Inf) +
      geom_vline(xintercept=-Inf) +
  
      scale_color_manual(values=c("#F8766D","#619CFF")) + 
  
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
ggsave("GTEx_forest_plot.png", units = "in", width=12, height=8, dpi=1200)
