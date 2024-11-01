# The script creates covariates and metabolite data set for looking at the correlations 
# Extracts propyl-4-hydroxybenzoate sulfate from both batch1 and batch2 data 
# for admixture mapping analysis and the final result in the supplement 

# Load neccessary packages
library(haven)
library(readxl)
library(tidyverse)

# create two datasets with variables needed for association analysis: 
# ID, 5 PCs, age, sex, eGFR, recruitment center, genetic analysis group
# and the 4 metabolites (need to choose them) based on each metabolomics batch. 
# (we need to remove individuals from batch 2 who have metabolomics in batch 1). 

pc_file <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/subject_annotation_2017-09-05.csv"
pheno_file <- "/Volumes/Sofer Lab/HCHS_SOL/Datasets/ms968_covariates_20200220.csv"
hh_file <- "/Volumes/Sofer Lab/HCHS_SOL/Datasets/household.sas7bdat"

## metabolomics data file:
## batch1 
metab_file_b1 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/2batch_combined_data_V1_only.xlsx"
# file in raw peaks
metab_file_b1_raw <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_Batch1/SOL_metabolomics_std_10202017.csv"
## batch2
metab_file_b2 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/batch2_data.xlsx"

dat1 <- read.csv(pc_file)
dat1$ID <- dat1$HCHS_ID
dat1 <- dat1[,c("GAC_scanID", "SUBJECT_ID", "ID", paste0("EV", 1:5), "gengrp6")]

dat2 <- read.csv(pheno_file)
dat2 <- dat2[,c("ID", "PSU_ID", "AGE", "GENDER", "GFRSCYS", "CENTER")]
head(dat2)

dat <- merge(dat1, dat2, by = "ID") # find overlapping individuals from both the pc file and the
# pheno file with other covariates needed

hh <- read_sas(hh_file) # household information 
hh$ID <- as.numeric(hh$ID)

dat <- merge(dat, hh, by = "ID") # merge the covariates and household info into one dataset 

## Prepare and add metabolomics data:

batch1_metab_sheet_name <- "data"
batch2_metab_sheet_name <- "batch2_batchnormalized_v1"
batch1_sample_info_sheet_name <- "sample.info"
batch2_sample_info_sheet_name <- "sample.info_batch2_v1"
metabolites_info_sheet_name <- "metabolites.info"

metab_vals_batch1 <- read_excel(metab_file_b1, sheet = batch1_metab_sheet_name)
metab_vals_batch1_raw <- read_csv(metab_file_b1_raw)
metab_vals_batch2 <- read_excel(metab_file_b2, sheet = batch2_metab_sheet_name)

sample_info_batch1 <- read_excel(metab_file_b1, sheet = batch1_sample_info_sheet_name)
sample_info_batch2 <- read_excel(metab_file_b2, sheet = batch2_sample_info_sheet_name)

# add sol Id and lab id for batch 1 
metab_vals_batch1 <- merge(metab_vals_batch1, 
                           sample_info_batch1[,c("PARENT_SAMPLE_NAME", "SOL_ID", "LABID", "VISIT",
                                                  "BATCH")], 
                           by = "PARENT_SAMPLE_NAME")

# add the metabolites info of propyl4hydroxybenzoatesulfate_std from the raw 
# csv file to metab_vals_batch1
metab_vals_batch1 <- merge(metab_vals_batch1, 
                           metab_vals_batch1_raw[,c("propyl4hydroxybenzoatesulfate_std", 
                                                    "LAB_ID")], 
                           by.x = "LABID", by.y = "LAB_ID")

metab_vals_batch1 <- metab_vals_batch1 |> filter(BATCH == "B01")

# rank-normalize the peaks and center it around the 1.0 to avoid non-positive values
metab_vals_batch1 <- metab_vals_batch1 |>
  mutate(propyl4hydroxybenzoatesulfate_std = if_else(is.na(propyl4hydroxybenzoatesulfate_std), 
                                                     min(propyl4hydroxybenzoatesulfate_std, na.rm= TRUE),
                                                     propyl4hydroxybenzoatesulfate_std))|>
  mutate(ranks = rank(propyl4hydroxybenzoatesulfate_std, ties.method = "average"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         `100006264` = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE)) + 0.5) |>
  select(-c(ranks, int_norm))
 
hist(metab_vals_batch1$`100006264`) # make sure the metabolite propyl4hydroxybenzoatesulfate 
# is normalized and cetered around 1. 
 
# add sol Id and lab id for batch 2
metab_vals_batch2 <- merge(metab_vals_batch2, 
                           sample_info_batch2[,c("PARENT_SAMPLE_NAME", "SOL_ID")], 
                           by = "PARENT_SAMPLE_NAME")

set.seed(1997)
# choose at random samples that are from the same individual (to have one sample per individual) 
random_index_b1 <- data.frame(index = 1:nrow(metab_vals_batch1 ), SOL_ID = metab_vals_batch1$SOL_ID)
for (id in unique(random_index_b1$SOL_ID)){
  row_inds <- which(random_index_b1$SOL_ID == id)
  if (length(row_inds) == 1) next
  selected_ind <- sample(row_inds, 1)
  random_index_b1 <- random_index_b1[-setdiff(row_inds, selected_ind),]
}
metab_vals_batch1 <- metab_vals_batch1[random_index_b1$index,]

# choose at random samples that are from the same individual (to have one sample per individual)
random_index_b2 <- data.frame(index = 1:nrow(metab_vals_batch2 ), SOL_ID = metab_vals_batch2$SOL_ID)
for (id in unique(random_index_b2$SOL_ID)){
  row_inds <- which(random_index_b2$SOL_ID == id)
  if (length(row_inds) == 1) next
  selected_ind <- sample(row_inds, 1)
  random_index_b2 <- random_index_b2[-setdiff(row_inds, selected_ind),]
}
metab_vals_batch2 <- metab_vals_batch2[random_index_b2$index,]

# remove overlapped individuals in both batches 
metab_vals_batch2 <- metab_vals_batch2[-which(is.element(metab_vals_batch2$SOL_ID, metab_vals_batch1$SOL_ID)),]

length(intersect(metab_vals_batch1$SOL_ID, metab_vals_batch2$SOL_ID)) # no overlaps, should be 0

# now we have the same set of metabolites in both batches, and different set of people

# now we want to choose metabolites with admixture mapping associations with Reynolds paper
# For chromosome 16 we have reported significant association for propyl 4-hydroxybenzoate sulfate 
# ChemID = 100006264 

metab_selected_batch1 <- metab_vals_batch1[,colnames(metab_vals_batch1) %in% c("SOL_ID", 100006264)] |>
  rename(ID = SOL_ID)

metab_selected_batch2 <- metab_vals_batch2[,colnames(metab_vals_batch2) %in% c("SOL_ID", 100006264)] |>
  rename(ID = SOL_ID)

# Combine all the information we need into two dataframes for each batch separately
# load overlapped individuals in RFMix and FLARE for both batches 
load("./admix_map_all/flare_batch1_SoLids.Rdata")
load("./admix_map_all/flare_batch2_SoLids.Rdata")

# join all covariate information that will be needed for admixture mapping 
# batch1 
colnames(flare_b1) <- c("SUBJECT_ID")
flare_b1$SUBJECT_ID <- paste0("SoL", flare_b1$SUBJECT_ID)
allinfo_batch1 <- merge(metab_selected_batch1, dat, by = "ID") # add metaboliate information 
allinfo_model_b1 <- merge(allinfo_batch1, flare_b1, by = "SUBJECT_ID") # only include overlapped individual in FLARE and RFMix

# batch2 
colnames(flare_b2) <- c("SUBJECT_ID")
flare_b2$SUBJECT_ID <- paste0("SoL", flare_b2$SUBJECT_ID)
allinfo_batch2 <- merge(metab_selected_batch2, dat, by = "ID") # add metaboliate information 
allinfo_model_b2 <- merge(allinfo_batch2, flare_b2, by = "SUBJECT_ID") # only include overlapped individual in FLARE and RFMix
                        
# Save the results into two separate RDS file
saveRDS(allinfo_model_b1, file = "./admix_map_all/low_corr_covariates_b1.RDS")
saveRDS(allinfo_model_b2, file = "./admix_map_all/low_corr_covariates_b2.RDS")
