args <- commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])

require(gdsfmt)
require(tidyverse)

# open the GDS file created from the FLARE3 ancestry mapping for the specified chromosome 
g <- openfn.gds(paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/FLARE3/FLARE3_ancestry_counts_chr",
index, ".gds"))

sampleid_FALRE <- as.data.frame(read.gdsn(index.gdsn(g, "sample.id")))
colnames(sampleid_FALRE) <- "FLARE_ID" 
sampleid_FALRE$FLARE_ID <- paste0("SoL", sampleid_FALRE$FLARE_ID)
rownames(sampleid_FALRE) <- sampleid_FALRE[,1]


# open the GDS file from local old ancestry inference 
gdsold <- openfn.gds("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/UW_GAC/lai_HGDP_1000G_comb_unique.gds") 

# Load the liftOver BED file for the old local ancestry interval annotations
intervals_liftover_nodup <- read.table("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/UW_GAC/intervals_37-38.bed")
colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 

print(paste0("working on ", index))

# Subsetting the intervals by chromosome
interval_old_lift_subset <- intervals_liftover_nodup %>% 
  dplyr::filter(chr == paste0("chr", index)) 

rownames(interval_old_lift_subset) <- interval_old_lift_subset$snpid

# reading SNP position (exact) file from the GDS node
snpnode <- index.gdsn(g, "snp.position")
snps_flare <- as.data.frame(read.gdsn(snpnode))
snps_flare$snpid <- read.gdsn(index.gdsn(g, "snp.id"))
colnames(snps_flare) <- c("snppos", "flare_snpid")

snp_interval_mapped <- c()
for(j in 1:nrow(interval_old_lift_subset)){
  indices <- which(snps_flare$snppos >= interval_old_lift_subset$pos_start[j] 
                   & snps_flare$snppos < interval_old_lift_subset$pos_end[j])
  tmp <- cbind(FLARE_snpid =  snps_flare$flare_snpid[indices], 
               FLARE_snppos = snps_flare$snppos[indices],
               old_snpid = rep(interval_old_lift_subset$snpid[j], 
                               length(indices)), 
               old_startpos = rep(interval_old_lift_subset$pos_start[j], 
                                  length(indices)), 
               old_endpos = rep(interval_old_lift_subset$pos_end[j], 
                                length(indices)))
  snp_interval_mapped <- rbind(snp_interval_mapped, tmp)
}
snp_interval_mapped <- as.data.frame(snp_interval_mapped)
          


# Compute the proportion of SNPs mapped onto the old intervals
print(nrow(snp_interval_mapped) / nrow(snps_flare)) 

readoldcounts <- function(ancestry){
  #Extract ancestry counts from the Old gds file (all counts from 15500 intervals)
  anc <- paste0("dosage_", ancestry)
  old_node <-  index.gdsn(gdsold, anc)
  old_counts <- as.data.frame(read.gdsn(old_node))
  colnames(old_counts) <- seq(1:ncol(old_counts)) #verified
  return(old_counts)
}

# take ancestry counts from the FLARE GDS file

readflarecounts <- function(ancestry){
  anc <- paste0(ancestry, "_counts")
  node <- index.gdsn(g, anc)
  counts <- as.data.frame(read.gdsn(node))
  colnames(counts) <- snps_flare$flare_snpid
  counts <- counts[ ,as.character(snp_interval_mapped$FLARE_snpid)]
  return(counts)
}      

#  250 intervals mapped for chr22

afr_old_counts_all <- readoldcounts("afr")
afr_old_counts_subset <- afr_old_counts_all[, unique(snp_interval_mapped$old_snpid)]

amer_old_counts_all <- readoldcounts("amer")
amer_old_counts_subset <- amer_old_counts_all[, unique(snp_interval_mapped$old_snpid)]

eur_old_counts_all <- readoldcounts("eur")
eur_old_counts_subset <- eur_old_counts_all[, unique(snp_interval_mapped$old_snpid)]

# indexing the SNPs with mapped intervals to the old inference 

afr_flare_counts_subset <- readflarecounts("afr")
amer_flare_counts_subset <- readflarecounts("amer")
eur_flare_counts_subset <- readflarecounts("eur")

closefn.gds(g)
closefn.gds(gdsold)

# loading the overlapped ID for FLARE
overlap_id <- read.csv("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/UW_GAC/old_new_overlap_ind.csv")[,2]
# loading the overlapped ID for old inference
merged_sample_old_clean <- read.csv("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/UW_GAC/merged_sample_old.csv")
load("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/UW_GAC/corr_matching_old_indices.Rdata")

# The corr_result was run separately by ancestries on AWS to save memory
# For different ancestry just comment out the code in other ancestries and change the file output names

mapped_afr_counts_old <- afr_old_counts_subset[old_ind, ]
rownames(mapped_afr_counts_old) <- merged_sample_old_clean$SUBJECT_ID
mapped_afr_counts_old <- mapped_afr_counts_old[overlap_id, ]

mapped_amer_counts_old <- amer_old_counts_subset[old_ind, ]
rownames(mapped_amer_counts_old) <- merged_sample_old_clean$SUBJECT_ID
mapped_amer_counts_old <- mapped_amer_counts_old[overlap_id, ]

mapped_eur_counts_old <- eur_old_counts_subset[old_ind, ]
rownames(mapped_eur_counts_old) <- merged_sample_old_clean$SUBJECT_ID
mapped_eur_counts_old <- mapped_eur_counts_old[overlap_id, ]

# FLARE counts
rownames(afr_flare_counts_subset) <- sampleid_FALRE$FLARE_ID
mapped_afr_counts_FLARE <- afr_flare_counts_subset[overlap_id,]

rownames(amer_flare_counts_subset) <- sampleid_FALRE$FLARE_ID
mapped_amer_counts_FLARE <- amer_flare_counts_subset[overlap_id,]

rownames(eur_flare_counts_subset) <- sampleid_FALRE$FLARE_ID
mapped_eur_counts_FLARE <- eur_flare_counts_subset[overlap_id,]
# sum(rownames(mapped_eur_counts_FLARE) != overlap_id) # should be 0, which means mapped correctly 

# Compute correlation 
corr_result <- c()

prop <- function(x){
  mean(x)/2
}

newvar <- function(x){
  var(x/2)
}

print("working on correlations")

for(i in colnames(mapped_afr_counts_old)){
  rowi <- c()
  counts_old_afr <- mapped_afr_counts_old[,i]
  counts_old_amer <- mapped_amer_counts_old[,i]
  counts_old_eur <- mapped_eur_counts_old[,i]
  
  ind <- which(as.character(snp_interval_mapped$old_snpid) == i)
  print(ind)
  counts_flare_afr <- mapped_afr_counts_FLARE[,as.character(snp_interval_mapped$FLARE_snpid[ind])]
  counts_flare_amer <- mapped_amer_counts_FLARE[,as.character(snp_interval_mapped$FLARE_snpid[ind])]
  counts_flare_eur <- mapped_eur_counts_FLARE[,as.character(snp_interval_mapped$FLARE_snpid[ind])]
  
  corr_afr_p <- as.data.frame(cor(counts_flare_afr, counts_old_afr, method = "pearson"))[, 1]
  corr_afr_s <- as.data.frame(cor(counts_flare_afr, counts_old_afr, method = "spearman"))[, 1]
  corr_amer_p <- as.data.frame(cor(counts_flare_amer, counts_old_amer, method = "pearson"))[, 1]
  corr_amer_s <- as.data.frame(cor(counts_flare_amer, counts_old_amer, method = "spearman"))[, 1]
  corr_eur_p <- as.data.frame(cor(counts_flare_eur, counts_old_eur, method = "pearson"))[, 1]
  corr_eur_s <- as.data.frame(cor(counts_flare_eur, counts_old_eur, method = "spearman"))[, 1]
  
  rowi <- cbind(Chromosome = rep(index, length(ind)), 
                Position = snp_interval_mapped$FLARE_snppos[ind], 
                Corr_afr_p = corr_afr_p,
                Corr_afr_s = corr_afr_s,
                Corr_amer_p = corr_amer_p,
                Corr_amer_s = corr_amer_s,
                Corr_eur_p = corr_eur_p,
                Corr_eur_s = corr_eur_s,
                Prop_afr_old = mean(counts_old_afr)/2,
                Prop_afr_new = apply(as.data.frame(counts_flare_afr), 2, prop),
                Prop_amer_old = mean(counts_old_amer)/2,
                Prop_amer_new = apply(as.data.frame(counts_flare_amer), 2, prop),
                Prop_eur_old = mean(counts_old_eur)/2,
                Prop_eur_new = apply(as.data.frame(counts_flare_eur), 2, prop),
                Var_afr_old = var(counts_old_afr/2),
                Var_afr_new = apply(as.data.frame(counts_flare_afr), 2, newvar),
                Var_amer_old = var(counts_old_amer/2),
                Var_amer_new = apply(as.data.frame(counts_flare_amer), 2, newvar),
                Var_eur_old = var(counts_old_eur/2),
                Var_eur_new = apply(as.data.frame(counts_flare_eur), 2, newvar)
  )
  corr_result <- rbind(corr_result, rowi)
  
}

corr_result <- as.data.frame(corr_result)

# customize file name to generate separate output files by ancestries
save(corr_result, file = paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/Correlation/corr_old_flare3_afr_chr", index, ".Rdata"))
#save(corr_result, file = paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/Correlation/corr_old_flare3_amer_chr", index, ".Rdata"))
#save(corr_result, file = paste0("/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/Correlation/corr_old_flare3_eur_chr", index, ".Rdata"))

