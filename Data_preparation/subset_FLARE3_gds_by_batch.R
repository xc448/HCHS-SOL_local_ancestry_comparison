library(tidyverse)
library(GENESIS)
library(gdsfmt)
library(GWASTools)

# The SoL IDs for metabolic batch 1 (discovery batch)
load("./flare_batch1_SoLids.Rdata")

# Read the FLARE 3 GDS file by chromosome 
# use gdsSubset function to extract individuals in batch 1 for each ancestry 
for(i in 1:22){
  print(paste0("working on chromosome "), i)
  flare3 <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\FLARE3_ancestry_counts_chr",i,".gds")
  gdsflare <- openfn.gds(flare3, readonly = FALSE)
  put.attr.gdsn(index.gdsn(gdsflare, "afr_counts"), "sample.order", val=1:nrow(flare_b1))
  put.attr.gdsn(index.gdsn(gdsflare, "amer_counts"), "sample.order", val=1:nrow(flare_b1))
  put.attr.gdsn(index.gdsn(gdsflare, "eur_counts"), "sample.order", val=1:nrow(flare_b1))
  put.attr.gdsn(flare_b1$flare_b1, "sample.id")
  closefn.gds(gdsflare)
  dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap\\flare3_b1_chr", i ,".gds")
  gdsSubset(flare3, sample.include=flare_b1$flare_b1, dir,
            compress = "LZMA_RA:1M")
}

#### Repeat the same procedure  for batch 2 ####
# The SoL IDs for metabolic batch 2 (replication batch)
load("./flare_batch2_SoLids.Rdata")

# Read the FLARE 3 GDS file -  by chromosome
for(i in 1:22){
  print(paste0("working on chromosome ", i))
  flare3 <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\FLARE3_ancestry_counts_chr",i,".gds")
  gdsflare <- openfn.gds(flare3, readonly = FALSE)
  put.attr.gdsn(index.gdsn(gdsflare, "afr_counts"), "sample.order", val=1:nrow(flare_b2))
  put.attr.gdsn(index.gdsn(gdsflare, "amer_counts"), "sample.order", val=1:nrow(flare_b2))
  put.attr.gdsn(index.gdsn(gdsflare, "eur_counts"), "sample.order", val=1:nrow(flare_b2))
  closefn.gds(gdsflare)
  dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap\\flare3_b2_chr", i ,".gds")
  gdsSubset(flare3, sample.include=flare_b2$flare_b2, dir,
            compress = "LZMA_RA:1M")
}