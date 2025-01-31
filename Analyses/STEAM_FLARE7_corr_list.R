args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])

require(gdsfmt)
require(data.table)
require(SNPRelate)
require(STEAM)
require(SeqArray)
require(tidyverse)
source("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/code/STEAM_customized_corr_FLARE7.R")

afr.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE7/FLARE7_afr_b2_chr', index, '.gds')
eur.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE7/FLARE7_eur_b2_chr', index, '.gds')
nam.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE7/FLARE7_amer_b2_chr', index, '.gds')

snps.dt <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE7/snps.dt_FLARE7_chr_', index, '.rds')

# run get_corr_chr
snps.dt_corr <- get_corr_chr_customized(chrom = index, 
                                        pop1.gds = afr.gds, pop2.gds = eur.gds, 
                                        pop3.gds = nam.gds, snps.dt = snps.dt)

save(snps.dt_corr, 
     file = paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/GDS_corr_list/FLARE7_corr/FLARE7_b2_chr_", index, ".Rdata"))

