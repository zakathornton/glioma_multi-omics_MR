
## Load modules

library(tidyverse)
library(readxl)

## Load data

data<-read_excel("Summary Table.xlsx")

## Filter to only include data which passes colocalisation and steiger filtering

data<-data[data$`Pass All Analyses`=="PASS",]

## Create a label position for the exposure (gene) and tissue

data$labelpositionExposure<-0.085
data$labelpositionTissue<-0.035

## Make subtype a factor

data$Subtype<-factor(data$Subtype, ordered = is.ordered(data$Subtype))

## Order the data by gene, then tissue, then subtype

data<-data[order(data$Gene, data$Tissue, data$Subtype),]

## Create a table position variable for each row

data$TablePosition<-seq(nrow(data),1,by=-1)

## Create the plot

plot<-ggplot(dat=data, 
             aes(x=TablePosition,y=`Odds Ratio`,ymin=`Lower 95% CI`,ymax=`Upper 95% CI`, colour=`Subtype`)) +

      geom_pointrange(aes(colour=factor(`Subtype`),group=`Gene`),shape=16,size=0.5) +
      geom_errorbar(width=0.3) +
  
      scale_y_log10(breaks=c(0.25,0.5,1,2,4,8), 
      labels = c("0.25","0.50","1","2","4","8"),
      expand = c(0,0)) +
  
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.02, ymax=0.125), fill="white", colour="black") +
      
      annotate("text", x=data$TablePosition, y=data$labelpositionExposure, label=data$Gene, size=3, fontface="italic") +
      annotate("text", x=data$TablePosition, y=data$labelpositionTissue, label=data$Tissue, size=3) +
        
      annotate("text", x=nrow(data)+1, y=data$labelpositionExposure, label="Exposure", size=4, fontface="bold") +
      annotate("text", x=nrow(data)+1, y=data$labelpositionTissue, label="Tissue", size=4, fontface="bold") +
      
      geom_segment(aes(x=-Inf, xend=-Inf, yend=10)) +
      geom_segment(aes(x=Inf, xend=Inf, yend=10)) +
  
      geom_hline(yintercept=Inf, linetype="solid") +
      geom_hline(yintercept=1, lty=2) +
      geom_vline(xintercept=Inf) +
      geom_vline(xintercept=-Inf) +
  
      coord_flip() +
        
      labs(y="Odds ratio (95% CI) for subtype risk per standard deviation change in gene expression",
           x="") +
  
      ggtitle("Mendelian randomisation results by subtype for the ten genes of interest in five brain tissues (eQTLs)") +

      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_text(size=13, face="bold"),
            axis.title.x = element_text(size=11),
            panel.background = element_rect(fill="grey95"),
            panel.grid.minor.y = element_line(linewidth = 1))
plot

## Save the plot

plot
ggsave("MB_forest_plot.png", units = "in", width=12, height=8, dpi=1200)
