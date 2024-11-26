library(rtracklayer)
# Aligned regions of ENCODE blacklist or annotations of gene clusters 
# from UCSC genome browser to the high correlation (>0.97 in all ancestries)
# SNP-LAI pair between FLARE3 and RFMix

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

load("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_RFMix_corr_by_chr/correlation_old_flare3_all_SNPs.Rdata")

# load the combined files for SNPs post-QC (R2 > 0.8, MAF >= 0.005)
# filtered_snps <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_combined.rds")
filtered_snps <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_combined.rds")

# filter the SNPs for the correlation file
dat_tot <-  dat_tot |> inner_join(filtered_snps, by =  c("Chromosome" = "chr", 
                                                         "Position" = "snppos"))

# Load ENCODE hg38 black list and gene clusters to R 
blacklist <- read.table("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/hg38-blacklist.bed" , header = FALSE, sep = "\t")
colnames(blacklist) <- c("Chromosome", "startpos", "endpos", "note")
blacklist$Chromosome <- gsub( "chr", "", blacklist$Chromosome)
#blacklist$description <- rep("ENCODE Blacklist", nrow(blacklist))
ucsc <- read.table("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/UCSC_problematic.bed" , header = FALSE, sep = "\t")
colnames(ucsc) <- c("Chromosome", "startpos", "endpos", "note", "score", "strand", 
                    "thickStart", "thickEnd","reserved", "description", "otherLoc")
ucsc <- ucsc |> filter(!grepl("_", Chromosome)) |>
  mutate(Chromosome = gsub("chr", "", Chromosome)) |> 
  filter(Chromosome %in% 1:22) |>
  select(Chromosome, startpos, endpos, note)# only include gene cluster info, 40 rows in total

problematic_regions <- rbind(blacklist, ucsc)
problematic_regions # 659 rows


# random sampling equal number of FLARE3-RFMix correlation pairs from high correlation regions (> 0.97)
# filter first 
tot_high_region <- c()
for(chr in 1:22){
  dat_chr <- dat_tot |> 
    filter(Chromosome == chr) |>
    filter(Corr_afr_p > 0.97 & Corr_amer_p > 0.97 & Corr_eur_p > 0.97)
  tot_high_region <- rbind(tot_high_region, dat_chr)
} # 43657 pairs 

# sample 
set.seed(42) # for replication
sampled_pairs_ind <- sample(nrow(tot_high_region), 31148)
sampled_pairs <- tot_high_region[sampled_pairs_ind, ]


# Helper functions for filtering intervals 
filter_intervals <- function(intervals, snp_positions, distance = 500000) {
  # Merge intervals and SNP positions based on chromosome (chr)
  merged <- merge(intervals, snp_positions, by = "Chromosome")
  # Calculate if SNP is within the predefined distance of start or end of each interval
  within_distance <- with(merged, 
                          (Position >= startpos - distance & Position <= endpos + distance))
  
  # Filter intervals that satisfy the condition
  filtered_intervals <- merged[within_distance, c("Chromosome", "startpos", "endpos", "note", "Position", "rsID")]
  
  return(unique(filtered_intervals))
}

# extract the FLARE regions (SNPs) that has a low correlation with RFMix above a certain correlation threshold (0.97)
# for a certain chromosome and plot both the SNPs and the mapped intervals from the ENCODE blacklist. 
plot_mapped_high_correlation <- function(dat, chrnum, distance = 500000, cor_thresh = 0.97){
  dat_chr <- dat |> 
    filter(Chromosome == chrnum) |>
    filter(Corr_afr_p > cor_thresh & Corr_amer_p > cor_thresh & Corr_eur_p > cor_thresh)
  
  problematic_chr <- problematic_regions |> filter(Chromosome == chrnum) 
  if(nrow(dat_chr) > 0) {
    filtered_intervals <- filter_intervals(problematic_chr, dat_chr)
    g <- dat_chr |> # this generate plots to see where the problematic regions are with respect to the high correlation regions
      ggplot(aes(x = Position)) +
      geom_point(aes(y = Corr_afr_p, color = "African"), size = 0.3, alpha = 0.4) +
      geom_point(aes(y = Corr_amer_p, color = "Amerindian"), size = 0.3, alpha = 0.4) +
      geom_point(aes(y = Corr_eur_p, color = "European"), size = 0.3, alpha = 0.4) +
      
      # Add vertical lines for blacklisted regions
      geom_vline(xintercept = c(filtered_intervals$startpos, filtered_intervals$endpos), 
                 color = "red", linetype = "dashed", alpha = 0.5) +
      
      theme_minimal() +
      labs(y = "Correlation", 
           title = paste0("High local ancstry correlations mapped to blacklisted regions, chr", chrnum)) +
      scale_color_manual(values = c("Amerindian" = "#D95F02", 
                                    "African" = "#1B9E77",
                                    "European" = "#7570B3"),
                         name = "Ancestry")
    return(list(g = g, filtered_intervals = filtered_intervals))
  }
}


problematic_intervals_highcor <- c()
plots <- list()
for(chr in 1:22){
  res <-  plot_mapped_high_correlation(dat_tot, chr)
  plots[[chr]] <- res$g 
  problematic_intervals_highcor <- rbind(problematic_intervals_highcor, res$filtered_intervals)
}

length(unique(in_blacklist$rsID))/length(unique(tot_high_region$rsID))


# how many high correlation pairs in ENCODE blacklist? --- 100% 
in_blacklist <- problematic_intervals_highcor|>
  filter(!note %in% c("IGH", "MHC", "IGL")) 

length(unique(in_blacklist$rsID))/length(unique(problematic_intervals_highcor$rsID))

# how many high correlation pairs in gene clusters? -- None
in_clusters <- problematic_intervals_highcor|>
  filter(note %in% c("IGH", "MHC", "IGL")) 

length(unique(in_clusters$rsID))/length(unique(problematic_intervals_highcor$rsID))
