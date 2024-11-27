# Create local ancestry GDS file for FLARE3 dataset from the original .txt file
# Load packages
library(tidyverse)
library(gdsfmt)
library(GWASTools)

setwd("R:\\Sofer Lab\\HCHS_SOL\\2024_FLARE_3pop")

counts <- function(ancestry_number, anc1, anc2) {
  # Output: a vector containing calculated ancestry counts for each 
  # individual (columns) for each SNP given the ancestry of interest.
  
  # Input: ancestry_number: number representing different ancestries (0-2)
  # anc1: copy #1 showing the ancestry information for each SNP (rows) on 
  # each individual(columns)
  # anc2:  copy #2 showing the ancestry information for each SNP (rows) on 
  # each individual(columns)
  tmp <- as.matrix((anc1 == ancestry_number) + (anc2 == ancestry_number)) # convert it into a matrix 
  # create a df with ancestry counts 
  return(tmp) 
}

# loop through the .txt local ancestry file for each chromosome
for(i in seq(22, 1)){
  chr_snp <- read.table(paste0("R:\\Sofer Lab\\HCHS_SOL\\2024_FLARE_3pop\\SOL_",
                               i, ".short.snp Iris Broce.txt"))
  id <- read.table("R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop\\SOL.id.txt",header = FALSE)
  snp <- read.table(header = FALSE, "R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop\\SOL.snp.txt") # snp info
  snpid <- match(chr_snp$V3,  snp$V3) 
  # snps are different from flare7 and flare3 bc the use of diff reference panels
  chrnum <- chr_snp$V1
  snppos <- chr_snp$V2
  sample_nums_id <- gsub( "SoL", x=id$V1, "")
  a1 <- chr_snp$V4 
  a2 <- chr_snp$V5
  
  # create a new gds file
  gfile <- createfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\FLARE3_ancestry_counts_chr", i ,".gds"))
  
  # adding snp information and alleles information to the GDS file
  add.gdsn(gfile, "snp.id", snpid, compress = "LZMA_RA:1M")
  add.gdsn(gfile, "snp.chromosome", chrnum, compress = "LZMA_RA:1M" )
  add.gdsn(gfile, "snp.position", snppos, compress = "LZMA_RA:1M" )
  add.gdsn(gfile, "sample.id", sample_nums_id, compress = "LZMA_RA:1M")
  add.gdsn(gfile, "a1", a1, compress = "LZMA_RA:1M" )
  add.gdsn(gfile, "a2", a2, compress = "LZMA_RA:1M" )
  
  # using LZMA_RA:1M as the compressing method
  anc1 <- read.table(paste0("SOL_", i, ".short.anc1 Iris Broce.txt"), header = FALSE) 
  anc2 <- read.table(paste0("SOL_", i, ".short.anc2 Iris Broce.txt"), header = FALSE)
  print("reading done")
  afr_counts <- t(counts(0, anc1, anc2))
  eur_counts <- t(counts(1, anc1, anc2))
  amer_counts <- t(counts(2, anc1, anc2))
  node_afr <- add.gdsn(gfile, "afr_counts", afr_counts, compress = "LZMA_RA:1M")
  node_amer <- add.gdsn(gfile, "amer_counts", amer_counts, compress = "LZMA_RA:1M")
  node_eur <- add.gdsn(gfile, "eur_counts", eur_counts, compress = "LZMA_RA:1M")
  
  afr_sum <- rowSums(afr_counts)
  amer_sum <- rowSums(amer_counts)
  eur_sum <- rowSums(eur_counts)
  
  flare3_global_proporiton <- as.data.frame(cbind(chr = rep(unique(chrnum), 
                                                            nrow(id)), 
                                                  id = id[,1], 
                                                  afr_sum = afr_sum, amer_sum
                                                  = amer_sum, eur_sum = eur_sum)) 
  
  closefn.gds(gfile)
  print(paste0("finishing on chromosome", i))
  save(flare3_global_proporiton, file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\global_sum_chr", 
                                               i, ".Rdata"))
}

# Calculating total number of SNPs in the FLARE3 file
nsnps  <- 0
for(i in 1:22){

  path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\FLARE3_ancestry_counts_chr", i, ".gds")
  gds <- openfn.gds(path)
  genotype_node <- index.gdsn(gds, "afr_counts")
  node_description <- objdesp.gdsn(genotype_node)
  
  # Adding up the number of SNPs
  nsnps <- nsnps + node_description$dim[2]
  closefn.gds(gds)
}
