library(gdsfmt)
library(tidyverse)

# load FLARE7 SNP information 
flare7 <- read.table("R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop/SOL.snp.txt")
colnames(flare7) <- c("chr", "snppos", "rsID", "a1", "a2")
flare7 <- flare7 |> dplyr::select(chr, snppos, rsID)

# record the proportion of SNPs left  
flare3_precentage <- c()
flare7_precentage <- c()

# iterate through all chromosomes to filter out SNPs with MAF < 0.005 and keep
# the SNPs with high imputation quality R2 >= 0.95
for(i in 22:1){
  print(paste0("working on chromosome ", i))
  imputation_info <- read.table(paste0("R:\\Sofer Lab\\HCHS_SOL\\2024_FLARE_3pop/MEGA_HCHS_SOL.chr", 
                                       i, ".qc Iris Broce.info")) |>
    select(V1, V2, V7)
  colnames(imputation_info) <- c("chr", "snppos", "quality")
  
  imputation_splitted <- separate(imputation_info, sep = ";", 
                                  col = quality, into  = c("passaf", 
                                                           "maf", "R2", 
                                                           "imputed", "R2_hat", 
                                                           "rsID")) 
  imputation_splitted <- imputation_splitted |> 
    mutate(R2 = as.numeric(gsub(".*=", "", R2))) |>
    filter(R2 >= 0.95) |>
    mutate(maf = as.numeric(gsub(".*=", "", maf)))|>
    filter(maf >= 0.005)
  
  
  #  match with FLARE3 and FLARE7 SNP info
  flare3 <- read.table(paste0("R:\\Sofer Lab\\HCHS_SOL\\2024_FLARE_3pop/SOL_",
                              i, ".short.snp Iris Broce.txt")) |>
    select(V1, V2, V3)
  colnames(flare3) <- c("chr",  "snppos",  "rsID")
  matched_snps_flare3  <- match(imputation_splitted$snppos,  flare3$snppos)
  matched_snps_flare3 <- matched_snps_flare3[!is.na(matched_snps_flare3)]
  flare3_filtered <- flare3[matched_snps_flare3,]
  flare3_precentage  <- c(flare3_precentage,
                          (nrow(flare3_filtered)/nrow(flare3)))
  saveRDS(flare3_filtered,
          file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_r.95_chr",
                        i, ".rds"))
  
  
  flare7_bychr <-  flare7 |> filter(chr ==  i)
  matched_snps_flare7 <- match(imputation_splitted$snppos,  flare7_bychr$snppos)
  matched_snps_flare7 <- matched_snps_flare7[!is.na(matched_snps_flare7)]
  flare7_filtered <- flare7_bychr[matched_snps_flare7,]
  flare7_precentage  <- c(flare7_precentage, nrow(flare7_filtered)/nrow(flare7_bychr))
  
  saveRDS(flare7_filtered, 
          file = paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_SNPs_filtered_r.95_chr", 
                        i, ".rds"))
  
}

# Create a combined file with all SNPs as well for FLARE3 and FLARE7
# flare_version takes only 3 or 7 as integers
combine_snpinfo <-  function(flare_version){
  combined_file  <- c()
  for(i in 1:22){
    snp_file <- readRDS(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE",
                               flare_version, "/FLARE", flare_version, "_SNPs_filtered_r.95_chr", 
                               i, ".rds"))
    
    combined_file  <- rbind(combined_file, snp_file)
    
  }
  return(as.data.frame(combined_file))
}

combined_flare3 <-  combine_snpinfo(3) 
combined_flare7 <-  combine_snpinfo(7) 

# Save the combined files
saveRDS(combined_flare3, file  =  "R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_r.95_combined.rds")
saveRDS(combined_flare7, file  =  "R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_SNPs_filtered_r.95_combined.rds")


# Compute the average and sd after filtering the SNPs based on MAF > 0.005 and R2 >= 0.95.
mean(flare3_precentage) 
sd(flare3_precentage) 
mean(flare7_precentage) 
sd(flare7_precentage) 
