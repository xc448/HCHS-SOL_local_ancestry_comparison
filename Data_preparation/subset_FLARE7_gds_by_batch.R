library(tidyverse)
library(GENESIS)
library(gdsfmt)
library(GWASTools)

# Read the FLARE 7 GDS file 
flare <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\ancestry_counts_FLARE.gds"
gdsflare <- openfn.gds(flare)

# The SoL IDs for metabolic batch 1 (discovery batch)
load("./flare_batch1_SoLids.Rdata")

# Subsetting the GDS file by extracting only the individuals from batch 1
# 200,000 SNPs / subfile, total # of SNPs = 5,105,005, 26 subfiles in total 

# The first 25 subfiles using gdsSubset
for(i in 1:25){
  dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b1_", i ,".gds")
  gdsSubset(flare, sample.include=flare_b1$flare_b1, dir,
            compress = "LZMA_RA:1M",snp.include = seq(200*(i-1)+1, 200*i))
}

# The 26th subfile 
dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b1_", 26 ,".gds")
gdsSubset(flare, sample.include=flare_b1$flare_b1, dir,
          compress = "LZMA_RA:1M",snp.include = seq(5000001, 5105005))


# Merging all the subfiles from batch 1 in a compressed GDS file 
gfile <- createfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b1_all.gds")

for( i in 1:26){
  tmp <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b1_", i, ".gds"))
  print(paste0("Starting with: ", i))
  if(i == 1){
    node_afr <- add.gdsn(gfile, "afr_counts", read.gdsn(index.gdsn(tmp,  "afr_counts")), compress = "LZMA_RA:1M")
    node_easia <- add.gdsn(gfile, "easia_counts", read.gdsn(index.gdsn(tmp, "easia_counts")), compress = "LZMA_RA:1M")
    node_eur <- add.gdsn(gfile, "eur_counts", read.gdsn(index.gdsn(tmp, "eur_counts")), compress = "LZMA_RA:1M")
    node_csasia <- add.gdsn(gfile, "csasia_counts", read.gdsn(index.gdsn(tmp, "csasia_counts" )), compress = "LZMA_RA:1M")
    node_amer <- add.gdsn(gfile, "amer_counts", read.gdsn(index.gdsn(tmp, "amer_counts")), compress = "LZMA_RA:1M")
    node_ocea <- add.gdsn(gfile, "ocea_counts", read.gdsn(index.gdsn(tmp, "ocea_counts")), compress = "LZMA_RA:1M")
    node_me <- add.gdsn(gfile, "me_counts", read.gdsn(index.gdsn(tmp, "me_counts")), compress = "LZMA_RA:1M")
  }
  else{
    append.gdsn(node_afr, read.gdsn(index.gdsn(tmp,  "afr_counts")))
    append.gdsn(node_easia, read.gdsn(index.gdsn(tmp, "easia_counts")))
    append.gdsn(node_eur, read.gdsn(index.gdsn(tmp, "eur_counts" )))
    append.gdsn(node_csasia, read.gdsn(index.gdsn(tmp, "csasia_counts" )))
    append.gdsn(node_amer, read.gdsn(index.gdsn(tmp, "amer_counts" )))
    append.gdsn(node_ocea, read.gdsn(index.gdsn(tmp, "ocea_counts" )))
    append.gdsn(node_me, read.gdsn(index.gdsn(tmp, "me_counts" )))
  }
  print(paste0("finishing up with: ", i))
}

############ Do the same thing for batch 2 ###########

# The SoL IDs for metabolic batch 2 (replication batch)
load("./flare_batch2_SoLids.Rdata")

# The first 25 subfiles using gdsSubset
for(i in 1:25){
  dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b2_", i ,".gds")
  gdsSubset(flare, sample.include=flare_b2$flare_b2, dir,
            compress = "LZMA_RA:1M",snp.include = seq(200000*(i-1)+1, 200000*i))
}

# The 26th subfile 
dir <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b2_", 26 ,".gds")
gdsSubset(flare, sample.include=flare_b2$flare_b2, dir,
          compress = "LZMA_RA:1M",snp.include = seq(5000001, 5105005))


# Merging all the subfiles from batch 1 in a compressed GDS file 
gfile <- createfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b2_all.gds")

for( i in 1:26){
  tmp <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\flare_b2_", i, ".gds"))
  print(paste0("Starting with: ", i))
  if(i == 1){
    node_afr <- add.gdsn(gfile, "afr_counts", read.gdsn(index.gdsn(tmp,  "afr_counts")), compress = "LZMA_RA:1M")
    node_easia <- add.gdsn(gfile, "easia_counts", read.gdsn(index.gdsn(tmp, "easia_counts")), compress = "LZMA_RA:1M")
    node_eur <- add.gdsn(gfile, "eur_counts", read.gdsn(index.gdsn(tmp, "eur_counts")), compress = "LZMA_RA:1M")
    node_csasia <- add.gdsn(gfile, "csasia_counts", read.gdsn(index.gdsn(tmp, "csasia_counts" )), compress = "LZMA_RA:1M")
    node_amer <- add.gdsn(gfile, "amer_counts", read.gdsn(index.gdsn(tmp, "amer_counts")), compress = "LZMA_RA:1M")
    node_ocea <- add.gdsn(gfile, "ocea_counts", read.gdsn(index.gdsn(tmp, "ocea_counts")), compress = "LZMA_RA:1M")
    node_me <- add.gdsn(gfile, "me_counts", read.gdsn(index.gdsn(tmp, "me_counts")), compress = "LZMA_RA:1M")
  }
  else{
    append.gdsn(node_afr, read.gdsn(index.gdsn(tmp,  "afr_counts")))
    append.gdsn(node_easia, read.gdsn(index.gdsn(tmp, "easia_counts")))
    append.gdsn(node_eur, read.gdsn(index.gdsn(tmp, "eur_counts" )))
    append.gdsn(node_csasia, read.gdsn(index.gdsn(tmp, "csasia_counts" )))
    append.gdsn(node_amer, read.gdsn(index.gdsn(tmp, "amer_counts" )))
    append.gdsn(node_ocea, read.gdsn(index.gdsn(tmp, "ocea_counts" )))
    append.gdsn(node_me, read.gdsn(index.gdsn(tmp, "me_counts" )))
  }
  print(paste0("finishing up with: ", i))
}
