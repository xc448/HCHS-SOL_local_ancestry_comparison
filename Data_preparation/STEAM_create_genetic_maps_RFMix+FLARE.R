library(gdsfmt)
library(tidyverse)

# load the map by chromosome 
# matching the genetic maps with the snp positions 
# FLARE3
gen_dist_markers_FLARE3 <- c()
for(i in 1:22){
  print(paste0("working on chromosome ", i))
  path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/plink.GRCh38.map/plink.chr",
  i, ".GRCh38.map") 
  working_chr <- read.table(path)
  colnames(working_chr) <- c("chr", "snpid", "cM", "bppos")
  gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr", i, ".gds"))
  pos <- as.data.frame(read.gdsn(index.gdsn(gds, "snp.position")))
  snpid <- as.data.frame(read.gdsn(index.gdsn(gds,  "snp.id")))
  chr <- read.gdsn(index.gdsn(gds,  "snp.chromosome"))
  closefn.gds(gds)
  colnames(snpid) <- "snpid"
  colnames(pos) <- "snppos"
  
  mapped_snp_pos <- c()
  for(i in 1:nrow(pos)){
    ind <- which.min(abs(working_chr$bppos - pos[i,])) # for each marker(snp), find the entry with the closest physical distance 
    # on the genetic map and assign the cM value to the specific marker
    recom_dist <- working_chr$cM[ind]
    mapped_snp_pos <- rbind(mapped_snp_pos, cbind("cM" = recom_dist,  "bppos" = pos$snppos[i]))
  }
  rownames(mapped_snp_pos) <-  snpid$snpid
  mapped_snp_pos <- as.data.frame(mapped_snp_pos) |> mutate(chr = chr) |>
    select(chr,  bppos, cM)
  gen_dist_markers_FLARE3 <-  rbind(gen_dist_markers_FLARE3,  mapped_snp_pos)
}

saveRDS(gen_dist_markers_FLARE3, 
        file = "R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/genetic_dist_FLARE3.rds")

# FLARE7

gen_dist_markers_FLARE7 <- c()
for(i in 2:1){
  print(paste0("working on chromosome ", i))
  path <- paste0("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/plink.GRCh38.map/plink.chr",
                 i, ".GRCh38.map") 
  working_chr <- read.table(path)
  colnames(working_chr) <- c("chr", "snpid", "cM", "bppos")
  gds <- openfn.gds(paste0("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_SNPs_filtered_chr", i, ".gds"))
  pos <- as.data.frame(read.gdsn(index.gdsn(gds, "snp.position")))
  snpid <- as.data.frame(read.gdsn(index.gdsn(gds,  "snp.id")))
  chr <- read.gdsn(index.gdsn(gds,  "snp.chromosome"))
  closefn.gds(gds)
  startpos <- min(which(chr ==  i))
  endpos <- max(which(chr  == i))
  colnames(snpid) <- "snpid"
  colnames(pos) <- "snppos"
  
  mapped_snp_pos <- c()
  for(i in 1:nrow(pos)){
    ind <- which.min(abs(working_chr$bppos - pos[i,])) # for each marker(snp), find the entry with the closest physical distance 
    # on the genetic map and assign the cM value to the specific marker
    recom_dist <- working_chr$cM[ind]
    mapped_snp_pos <- rbind(mapped_snp_pos, cbind("cM" = recom_dist,  "bppos" = pos$snppos[i]))
  }
  rownames(mapped_snp_pos) <-  snpid$snpid
  mapped_snp_pos <- as.data.frame(mapped_snp_pos) |> mutate(chr = chr) |>
    select(chr,  bppos, cM)
  gen_dist_markers_FLARE7 <-  rbind(mapped_snp_pos, gen_dist_markers_FLARE7)
}

saveRDS(gen_dist_markers_FLARE7, "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE7_gds_steam/genetic_dist_FLARE7.rds")

# RFMix
# load the ancestry interval files after liftover  
intervals_liftover_nodup <- read.table("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\intervals_37-38.bed")
colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 
intervals_liftover_nodup <- intervals_liftover_nodup |> filter(chr != "chrX") 
# 14753 mapped blocks

gds <- openfn.gds("./uwgds_b2_subset_37")
snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))
snp.chromosome <- read.gdsn(index.gdsn(gds,"snp.chromosome"))
snp.position <- read.gdsn(index.gdsn(gds,"snp.position"))

intervals_liftover_nodup$midpoint <- intervals_liftover_nodup$pos_start + (intervals_liftover_nodup$pos_end-
  intervals_liftover_nodup$pos_start+1)/2

gen_dist_markers_rfmix <- c()
for(i in 1:22){
  print(paste0("working on chromosome ", i))
  interval <-  intervals_liftover_nodup[intervals_liftover_nodup$chr == paste0("chr",i),]
    path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/plink.GRCh38.map/plink.chr", i, ".GRCh38.map") 
  working_chr <- read.table(path)
  colnames(working_chr) <- c("chr", "snpid", "cM", "bppos")
  
  mapped_snp_pos <- c()
  for(j in 1:nrow(interval)){
    ind <- which.min(abs(working_chr$bppos - interval$midpoint[j])) 
    # for each marker(snp), find the entry with the closest physical distance 
    # on the genetic map and assign the cM value to the specific marker
    recom_dist <- working_chr$cM[ind]
    mapped_snp_pos <- rbind(mapped_snp_pos, cbind("cM" = recom_dist,  
                                                  "bppos" = interval$midpoint[j]))
  }
  rownames(mapped_snp_pos) <-  interval$snpid
  mapped_snp_pos <- as.data.frame(mapped_snp_pos) |> 
    mutate(chr = rep(i,  nrow(mapped_snp_pos))) |>
    select(chr,  bppos, cM)
  gen_dist_markers_rfmix <-  rbind(mapped_snp_pos, gen_dist_markers_rfmix)
}

saveRDS(gen_dist_markers_rfmix, 
        file = "R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/RFMix_gds_by_chr/genetic_dist_RFMix.rds" )


