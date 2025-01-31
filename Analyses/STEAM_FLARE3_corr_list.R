args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])

require(gdsfmt)
require(data.table)
require(SNPRelate)
require(STEAM)
require(SeqArray)
require(tidyverse)
source("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/code/20240802_STEAM_customized_corrfunc.R")

# FLARE3_map <- readRDS("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/genetic_dist_FLARE3.rds")

afr.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/FLARE3_afr_b2_chr', index, '.gds')
eur.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/FLARE3_eur_b2_chr', index, '.gds')
nam.gds <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/FLARE3_amer_b2_chr', index, '.gds')

snps.dt <- paste0('/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/snps.dt_FLARE3_chr_', index, '.rds')

# run get_corr_chr
# FLARE3_map_chr <-  FLARE3_map %>% filter(chr == index)
snps.dt_corr <- get_corr_chr_customized(chrom = index, 
                                          pop1.gds = afr.gds, pop2.gds = eur.gds, 
                                          pop3.gds = nam.gds, snps.dt = snps.dt)

save(snps.dt_corr, 
     file = paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/GDS_corr_list/FLARE3_corr/FLARE3_b2_chr_", index, ".Rdata"))