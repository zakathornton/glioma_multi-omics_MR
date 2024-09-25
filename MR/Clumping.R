## Load modules

library(ieugwasr)
library(TwoSampleMR)
library(foreach)
library(dplyr)

## Load the oxposure data, and format for clumping/MR analysis

exp_data<-read.delim("brainQTL_data.txt", sep=" ")
exp_data<-format_data(exp_data,type="exposure", header=T, chr_col = "CHR", pos_col = "POS", snp_col="SNP", phenotype_col="GENE", effect_allele_col="A1", other_allele_col="A2",
                        beta_col="BETA", se_col="SE", pval_col="P", eaf_col="FREQA1", samplesize_col="N")
exp_data$rsid<-exp_data$SNP
exp_data$pval<-exp_data$pval.exposure
exp_data$id<-exp_data$id.exposure

## Clump the exposures one-by-one, and save to a file

idx <- unique(exp_data$exposure)
exp_data <- parallel::mclapply(1:length(idx), function(i)
{
        if(file.exists(paste0("brainQTL_clumping_0.001_", idx[i],".txt")))
        return(NA)

        attempt <- try(dat <- ieugwasr::ld_clump(exp_data[exp_data$exposure == idx[i], ],
                                             clump_r2 = 0.001,
                                             plink_bin = "/mnt/storage/software/apps/Plink-1.90/plink",
                                             bfile = "GRCh37_EUR_dupvar"),
                                             silent = T)

         if (class(attempt) == "try-error" || nrow(dat) < 1) {
        warning("Error with clumping, or no SNPs retained thereafter for exposure: ", idx[i])
        return(NULL)
    }
     	write.table(dat, file = paste0("brainQTL_clumping_0.001_", idx[i],".txt"), quote = F, row.names = F)
        return(dat)
}, mc.cores = 2)
