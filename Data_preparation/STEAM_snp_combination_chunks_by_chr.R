require(gdsfmt)
require(data.table)
require(SNPRelate)
require(STEAM)
require(SeqArray)
require(tidyverse)
library(RcppAlgos)
source("./20240802_STEAM_customized_corrfunc.R")


sample_by_chunks <- function(chr){
  map.df <-  readRDS(paste0("./AWS/FLARE3_corr_list_STEAM/FLARE3_map_chr",  chr, ".rds"))
  
  # list all possible pairs of loci
  map.df$snp.id <- rownames(map.df)
  
  chunk_size <- 10000
  snp_chunks <- split(map.df$snp.id, ceiling(seq_along(map.df$snp.id) / chunk_size))
  
  snps.pairs.list <- list()
  
  # Initialize a list to store data.tables for each bin
  bin_tables <- list()
  
  for (i in seq_along(snp_chunks)) {
    for (j in i:length(snp_chunks)) {
      gc()
      print(i)
      if (i == j) {
        print(paste0("j is ", j))
        snps.pairs <- comboGeneral(snp_chunks[[i]], 2)
        gc()
      } else {
        print(paste0("j is ", j))
        snps.pairs <- comboGrid(snp_chunks[[i]], snp_chunks[[j]])
        gc()
      }
      snps.dt.chunk <- data.table(snps.pairs)
      rm(snps.pairs)
      gc()
      # Convert the pairs to a data.table
      names(snps.dt.chunk) <- c('snp1', 'snp2')
      
      # Calculate genetic distances
      snps.dt.chunk[, cM := get_dist(snp1, snp2, map.df)]
      gc()
      
      # Convert to recombination fractions
      snps.dt.chunk[, theta := L_to_theta(cM)]
      gc()
      
      # Determine the distance bin
      snps.dt.chunk[, bin := round(cM / binsize, 0) * binsize]
      gc()
      
      # Split the chunk by bin and store in the list
      split_chunks <- split(snps.dt.chunk, by = "bin")
      
      # Append the split chunks to their respective bin tables
      for (bin in names(split_chunks)) {
        print(bin)
        if (!is.null(bin_tables[[bin]])) {
          bin_tables[[bin]] <- 1
          table <- readRDS(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/chr_", chr, "_", bin, "snps.rds"))
          #print("appending")
          t <- rbind(table, split_chunks[[bin]], fill = TRUE)
          rm(table)
          gc()
          saveRDS(t, file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/chr_", chr,"_",  bin, "snps.rds"))
          rm(t)
          gc()
        } else {
          #print("creating")
          bin_tables[[bin]] <- 1
          saveRDS(split_chunks[[bin]], file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/chr_", chr, "_", bin, "snps.rds"))
        }
        gc()
      }
      
      
      # Remove the chunk to free memory
      rm(snps.dt.chunk, split_chunks)
      gc()
    }
  }
  
  snps.dt <- data.table(snp1 = character(), snp2 = character())
  for(names in names(bin_tables)){
    snps_dt_bin <- readRDS(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/chr_", chr, "_", names, "snps.rds"))
    print(names)
    # Initialize an empty data.table for the final sampled results
    if (nrow(snps_dt_bin) > 0) {
      idx <- sample(nrow(snps_dt_bin), size = min(20, nrow(snps_dt_bin)), replace = FALSE)
      sampled_chunk <- snps_dt_bin[idx, .(snp1, snp2)]
      snps.dt <- rbind(snps.dt, sampled_chunk, fill = TRUE)
    }
    rm(snps_dt_bin)
    gc()
  }
  
  saveRDS(snps.dt, 
          file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/snps.dt_FLARE3_chr_", chr, ".rds" ))                                      
  gc()
  cat('done with chromosome', chr, '\n')
}


for(chromosome in 22:1){
  sample_by_chunks(chromosome)
}

