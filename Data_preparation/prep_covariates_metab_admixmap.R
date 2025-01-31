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
metab_file_b1 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/2batch_combined_data_V1_only.xlsx"
metab_file_b2 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/batch2_data.xlsx"

dat1 <- read.csv(pc_file)
dat1$ID <- dat1$HCHS_ID
dat1 <- dat1[,c("GAC_scanID", "SUBJECT_ID", "ID", paste0("EV", 1:5), "gengrp6")]

dat2 <- read.csv(pheno_file)
dat2 <- dat2[,c("ID", "PSU_ID", "AGE", "GENDER", "GFRSCYS", "CENTER")]
head(dat2)

dat <- merge(dat1, dat2, by = "ID")

hh <- read_sas(hh_file)
hh$ID <- as.numeric(hh$ID)

dat <- merge(dat, hh, by = "ID")


## Prepare and add metabolomics data:

batch1_metab_sheet_name <- "data"
batch2_metab_sheet_name <- "batch2_batchnormalized_v1"
batch1_sample_info_sheet_name <- "sample.info"
batch2_sample_info_sheet_name <- "sample.info_batch2_v1"
metabolites_info_sheet_name <- "metabolites.info"


metab_vals_batch1 <- read_excel(metab_file_b1, sheet = batch1_metab_sheet_name)
metab_vals_batch2 <- read_excel(metab_file_b2, sheet = batch2_metab_sheet_name)

# use metabolites that were present in batch 1
metab_vals_batch2 <- metab_vals_batch2[, colnames(metab_vals_batch1)]

sample_info_batch1 <- read_excel(metab_file_b1, sheet = batch1_sample_info_sheet_name)
sample_info_batch2 <- read_excel(metab_file_b2, sheet = batch2_sample_info_sheet_name)


metab_vals_batch1 <- merge(metab_vals_batch1, 
                           sample_info_batch1[,c("PARENT_SAMPLE_NAME", "SOL_ID")], 
                           by = "PARENT_SAMPLE_NAME")

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

metab_vals_batch2 <- metab_vals_batch2[-which(is.element(metab_vals_batch2$SOL_ID, metab_vals_batch1$SOL_ID)),]

length(intersect(metab_vals_batch1$SOL_ID, metab_vals_batch2$SOL_ID)) #0
# now we have the same set of metabolites in both batches, and different set of people

# now we want to choose 4 metabolites with admixture mapping associations with Reynolds paper
# The four metabolites: PC 16:0/20:4 (chem ID 100008914), PE 16:0/20:4 (chem ID 100008990), 
# N-acetylarginine (chem ID 100001266), 3-aminoisobutyrate (chem ID 1114)

metab_selected_batch1 <- metab_vals_batch1[,colnames(metab_vals_batch1) %in% c("SOL_ID", 1114, 100001266, 100008914, 100008990)] |>
  rename(ID = SOL_ID)

metab_selected_batch2 <- metab_vals_batch2[,colnames(metab_vals_batch2) %in% c("SOL_ID", 1114, 100001266, 100008914, 100008990)] |>
  rename(ID = SOL_ID)

# Combine all the information we need into two dataframes for each batch separately
allinfo_batch1 <- merge(metab_selected_batch1, dat, by = "ID") 

allinfo_batch2 <- merge(metab_selected_batch2, dat, by = "ID")


# Finally create HH matrix using HH ID 
# batch1 
hh.matrix_b1 <- outer(allinfo_batch1$HH_ID, allinfo_batch1$HH_ID, FUN = "==") * 1
colnames(hh.matrix_b1) <- allinfo_batch1$GAC_scanID
rownames(hh.matrix_b1) <- allinfo_batch1$GAC_scanID

sum(rowSums(hh.matrix_b1) > 1) # 1044 individuals have the shared household with others in batch1 

# batch2
hh.matrix_b2 <- outer(allinfo_batch2$HH_ID, allinfo_batch2$HH_ID, FUN = "==") * 1

sum(rowSums(hh.matrix_b2) > 1) 
# 366 individuals have the shared household with others in batch2
colnames(hh.matrix_b2) <- allinfo_batch2$GAC_scanID
rownames(hh.matrix_b2) <- allinfo_batch2$GAC_scanID

# output matrices
save(hh.matrix_b1, file = "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/hh.matrix_b1.Rdata")
save(hh.matrix_b2, file = "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/hh.matrix_b2.Rdata")

# Create block matrix matrix
# batch1 
block.matrix_b1 <- outer(allinfo_batch1$PSU_ID, allinfo_batch1$PSU_ID, FUN = "==") * 1
colnames(block.matrix_b1) <- allinfo_batch1$GAC_scanID
rownames(block.matrix_b1) <- allinfo_batch1$GAC_scanID

# batch2
block.matrix_b2 <- outer(allinfo_batch2$PSU_ID, allinfo_batch2$PSU_ID, FUN = "==") * 1 
colnames(block.matrix_b2) <- allinfo_batch2$GAC_scanID
rownames(block.matrix_b2) <- allinfo_batch2$GAC_scanID

# output matrices
save(block.matrix_b1, file = "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/block.matrix_b1.Rdata")
save(block.matrix_b2, file = "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/block.matrix_b2.Rdata")

# Save prepared data
write.csv(allinfo_batch1, "../Data/metab_admixture_association_batch1.csv" )
write.csv(allinfo_batch2, "../Data/metab_admixture_association_batch2.csv" )

