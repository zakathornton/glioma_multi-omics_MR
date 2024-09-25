## Load modules

library(dplyr)
library(readr)
library(readxl)

## Create the function for extracting important data

get_important_hpa_columns <- function(hpa_tb){
  hpa_tb <- hpa_tb |> 
    dplyr::select(Gene,
                  `Gene synonym`,
                  Ensembl,
                  `Gene description`,
                  Uniprot, 
                  `Protein class`,
                  `Molecular function`,
                  `RNA tissue specificity`, 
                  `RNA tissue specific nTPM`,
                  `RNA brain regional specificity`, 
                  `RNA brain regional specific nTPM`,
                  `RNA single cell type specificity`,
                  `RNA single cell type specific nTPM`) 
  return(hpa_tb)
}

## Read HPA table of gene annotations

hpa_data <- read_tsv("proteinatlas.tsv") |> 
  get_important_hpa_columns()

## Keep only genes that are expressed in the brain 

flt_data <- hpa_data |> 
  filter(`RNA brain regional specificity` != "Not detected")

## Remove broadly expressed genes

flt_data <- flt_data |> 
  filter(`RNA tissue specificity` != "Low tissue specificity")

## Add brain specificity columns

flt_data <- flt_data |> 
  mutate(`Brain enriched/enhanced` = grepl("brain", `RNA tissue specific nTPM`))

## Pick out genes of interest

dat<-flt_data[flt_data$Gene%in%c("CDKN2B","EGFR","CEP192","FAIM","D2HGDH","IDH1","HBEGF","HEATR3","GALNT6","RTEL1","MDM4","JAK1","PICK1","RAVER2","TERT","SLC8A1","PHLDB1"),]

## Save the table

write.table(dat, "~/OneDrive - University of Bristol/QTL/Enrichment Analysis/HPA Lookup.csv", quote = F, row.names = F)
