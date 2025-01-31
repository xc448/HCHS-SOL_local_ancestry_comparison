# load neccessary packages 
library(tidyverse) 
library(GENESIS)
library(gdsfmt)
library(GWASTools)
library(SNPRelate) 

# load the kinship matrix -- already generated previously 
#kinpath <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files"
kinpath <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files"
kin.mat <- getobj(file.path(kinpath, "pcreap_phi_matrix.RData"))

#######################BATCH 1#############################

# load the household matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch1
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\hh.matrix_b1.Rdata")

# load block matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch1
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\block.matrix_b1.Rdata")

# Load SNP Annotations - RFMix
snpannot <-  "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries\\lai_HGDP_1000G_comb_snpAnnot_unique.RData"
#snpannot <-  "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries/lai_HGDP_1000G_comb_snpAnnot_unique.RData"
load(snpannot)

# Load covariates for b1

load("./admix_map_all/covariates_b1_model.Rdata")
original_scanid_b1 <- covariate_b1$GAC_scanID
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=original_scanid_b1, stringsAsFactors=FALSE))

covariate_b1 <- readRDS("./admix_map_all/low_corr_covariates_b1.RDS")
colnames(covariate_b1)[3:6] <- paste0("X", colnames(covariate_b1)[3:6])

# match individuals that are in the filtered batch 1, with correct order based on the GDS file
# the order of the samples is ascending by SoLID
tmp <- covariate_b1$GAC_scanID

# subset and order kinship matrix by both rows and column names
kin.mat_b1 <-kin.mat[match((as.character(tmp)),rownames(kin.mat)), ] 
kin.mat_b1 <-kin.mat_b1[, match(as.character(tmp), colnames(kin.mat_b1))]

# subset and order household matrix by both rows and column names
hh.matrix_b1 <- hh.matrix_b1[match(as.character(tmp),rownames(hh.matrix_b1)),]
hh.matrix_b1 <-hh.matrix_b1[, match(as.character(tmp), colnames(hh.matrix_b1))]

# subset and order block matrix by both rows and column names
block.matrix_b1 <- block.matrix_b1[match(as.character(tmp),rownames(block.matrix_b1)),]
block.matrix_b1 <-block.matrix_b1[, match(as.character(tmp), colnames(block.matrix_b1))]

# used the null models to run the joint test for associations between ancestries 
# at each locus and metabolites using a Wald test (default).
rownames(covariate_b1) <- as.character(covariate_b1$GAC_scanID)
x_b1 <- covariate_b1[as.character(tmp),c(seq(8,13), seq(15, 18))] # make sure the 

# define a function for fitting the null model 
nullmod <- function(outcome, covdat, kinmatrix, hhmatrix, blockmatrix){
  covMatList <- list(HH = hhmatrix, kinship = kinmatrix, block = blockmatrix)
  x <- cbind(covariate_b1[as.character(tmp), outcome], covdat)
  colnames(x)[1] <- outcome
  mod <- fitNullModel(x = x, outcome = outcome, 
                      covars=c("AGE","GENDER","EV1","EV2","EV3","EV4","EV5","gengrp6", "GFRSCYS", "CENTER"),
                      cov.mat = covMatList)
  return(mod)
}
# Metabolite propyl 4-hydroxybenzoate sulfate - for chromosome 16, 100006264 
nullmod_b1_x6264 <- nullmod("X100006264", x_b1, kin.mat_b1, hh.matrix_b1, block.matrix_b1)


# Save the null model
save(nullmod_b1_x6264, file = "./admix_map_all/admixmap_low_corr/nullmod_b1_x6264.Rdata")

# open the gds file for RFMix 
gds_old_b1 <- openfn.gds("./uwgds_b1_subset_37", readonly = FALSE) # This would also have 3861 individuals

# runassoc function for RFMix gds file with local ancestry counts
runassoc <- function(ancestry, gdsfile, nullmod, inference = "old"){
  gds.reader <- GdsGenotypeReader(gdsfile, 
                                  genotypeVar=paste0("dosage_", ancestry))
  genodata <- GenotypeData(gds.reader, scanAnnot=scanAnnot)
  iterator <- GenotypeBlockIterator(genodata, snpBlock = 100)
  geno <- getGenotypeSelection(iterator)
  colnames(geno) <- original_scanid_b1
  geno <- geno[, colnames(geno) %in% tmp]
  testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
  # this one uses Score test 
  res <- testres
  
  # loop through all the snp blocks
  while(iterateFilter(iterator)){
    geno <- getGenotypeSelection(iterator)
    colnames(geno) <- original_scanid_b1
    geno <- geno[, colnames(geno) %in% tmp]
    # table(geno) only gives 0, 1, 2
    testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
    res <- rbind(testres, res)
  }
  return(res)
}

# FLARE scan annotation 
# Sort by subject ID now
gdsflare7_b1 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b1_all.gds")
flareid <- read.gdsn(index.gdsn(gdsflare7_b1, "sample.id"))
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=tmp, stringsAsFactors=FALSE))

scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
  scanID=flareid, stringsAsFactors=FALSE)) 
closefn.gds(gdsflare7_b1)

# flare3 admixture mapping function 
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare3 <- function(ancestry, nullmod, metab, from, to){
  for(i in seq(from=from, to = to)){
    print(paste0("working on chromosome", i))
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    gdsflare_b1 <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap/flare3_b1_chr", 
                                     i, ".gds" ))
    gds <- GdsGenotypeReader(gdsflare_b1,  genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
    chrom <- getChromosome(gds)[1:length(filtered_id)]
    snppos <- as.data.frame(getPosition(gds))
    rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_b1, "snp.id"))
    snppos <- snppos[rownames(snppos) %in% filtered_id,]
    geno <- getGenotype(genodata)
    rownames(geno) <- read.gdsn(index.gdsn(gdsflare_b1, "snp.id"))
    geno <- geno[rownames(geno) %in% filtered_id,]
    colnames(geno) <- paste0("SoL", flareid)
    geno <- geno[, colnames(geno) %in% covariate_b1$SUBJECT_ID]
    
    closefn.gds(gdsflare_b1)
    closefn.gds(filtered_gds)
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_low_corr/", metab, "FLARE3_", 
                            ancestry ,"_b1_fixed_chr", i, ".Rdata"))
    gc()
  }
}


# flare7 admixture mapping function 
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare7 <- function(ancestry, from, to, nullmod, metab){
  for(i in seq(from=from, to = to) ){
    print(paste0("working on chromosome", i))
    gdsflare_b1 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b1_all.gds")
    gds <- GdsGenotypeReader(gdsflare_b1,  genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
    chrom <- getChromosome(gds)
    chromid <- which(chrom == i)
    snppos <- as.data.frame(getPosition(gds)[chromid])
    
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    
    chrom <- getChromosome(gds)[1:length(filtered_id)]
    rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_b1, "snp.id"))[chromid]
    snppos <- snppos[rownames(snppos) %in% filtered_id,]
    
    geno <- getGenotype(genodata, snp = c(min(chromid), 
                                          max(chromid) - min(chromid) +1))
    
    rownames(geno) <- read.gdsn(index.gdsn(gdsflare_b1, "snp.id"))[chromid]
    geno <- geno[rownames(geno) %in% filtered_id,]
    colnames(geno) <- paste0("SoL", flareid)
    geno <- geno[, colnames(geno) %in% covariate_b1$SUBJECT_ID]
    closefn.gds(gdsflare_b1)
    closefn.gds(filtered_gds)
    
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_low_corr/", metab, "FLARE7_", 
                            ancestry ,"_b1_fixed_chr", i, ".Rdata"))
    gc()
  }
}

# Run admixture mapping
# Metabolite propyl 4-hydroxybenzoate sulfate - for chromosome 16, 100006264
# Using RFMix local ancestry inference 
res_afr_x6264 <-runassoc("afr", gds_old_b1, nullmod_b1_x6264)
res_amer_x6264 <- runassoc("amer", gds_old_b1, nullmod_b1_x6264)

# Using FLARE3 local ancestry inference
runassoc_flare3("afr", from = 16, to = 16, nullmod = nullmod_b1_x6264,
                metab = "x6264")
runassoc_flare3("amer", from = 16, to = 16, nullmod = nullmod_b1_x6264,
                metab = "x6264")

# Using FLARE7 local ancestry inference 
runassoc_flare7("afr", from = 16, to = 16, nullmod = nullmod_b1_x6264,
                metab = "x6264")
runassoc_flare7("amer", from = 16, to = 16, nullmod = nullmod_b1_x6264,
                metab = "x6264")


# Save RFMix admixture mapping results from batch 1:
saveRDS(res_afr_x6264, file = "./admix_map_all/admixmap_low_corr/RFMix_res_afr_x6264_b1_fixed.rds")
saveRDS(res_amer_x6264, file = "./admix_map_all/admixmap_low_corr/RFMix_res_amer_x6264_b1_fixed.rds")


#######################BATCH 2#############################
# load the household matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch2
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\hh.matrix_b2.Rdata")

# load block matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch2
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\block.matrix_b2.Rdata")

# Load covariates for b2
covariate_b2 <- readRDS("./admix_map_all/low_corr_covariates_b2.RDS")
colnames(covariate_b2)[3:6] <- paste0("X", colnames(covariate_b2)[3:6])
original_scanid <- covariate_b2$GAC_scanID
covariate_b2 <- covariate_b2[!is.na(covariate_b2$X100001208),] # remove NA from one individual, 1643 complete cases

# match individuals that are in the filtered batch 2, with correct order based on the GDS file
# the order of the samples is ascending by SoLID
tmp2 <- covariate_b2$GAC_scanID # tmp2 is the GAC_scanID indices

# subset and order kinship matrix by both rows and column names
kin.mat_b2 <-kin.mat[match((as.character(tmp2)),rownames(kin.mat)), ] 
kin.mat_b2 <-kin.mat_b2[, match(as.character(tmp2), colnames(kin.mat_b2))]

# subset and order household matrix by both rows and column names
hh.matrix_b2 <- hh.matrix_b2[match(as.character(tmp2),rownames(hh.matrix_b2)),]
hh.matrix_b2 <-hh.matrix_b2[, match(as.character(tmp2), colnames(hh.matrix_b2))]

# subset and order block matrix by both rows and column names
block.matrix_b2 <- block.matrix_b2[match(as.character(tmp2),rownames(block.matrix_b2)),]
block.matrix_b2 <-block.matrix_b2[, match(as.character(tmp2), colnames(block.matrix_b2))]

# used the null models to run the joint test for associations between ancestries 
# at each locus and metabolites using a Wald test (default).
rownames(covariate_b2) <- as.character(covariate_b2$GAC_scanID)
x_b2 <- covariate_b2[as.character(tmp2),c(seq(8,13), seq(15, 18))] 

# define a function for fitting the null model 
nullmod <- function(outcome, covdat, kinmatrix, hhmatrix, blockmatrix){
  covMatList <- list(HH = hhmatrix, kinship = kinmatrix, block = blockmatrix)
  x <- cbind(covariate_b2[as.character(tmp2), outcome], covdat)
  colnames(x)[1] <- outcome
  mod <- fitNullModel(x = x, outcome = outcome,
                      covars=c("AGE","GENDER","EV1","EV2","EV3","EV4","EV5","gengrp6", "GFRSCYS", "CENTER"),
                      cov.mat = covMatList)
  return(mod)
}


# Metabolite propyl 4-hydroxybenzoate sulfate - for chromosome 16, 100006264 
nullmod_b2_x6264 <- nullmod("X100006264", x_b2, kin.mat_b2, hh.matrix_b2, block.matrix_b2)


# Save the null model
save(nullmod_b2_x6264, file = "./admix_map_all/admixmap_low_corr/nullmod_b2_x6264.Rdata")

# open the gds file for RFMix 
gds_old_b2 <- openfn.gds("./uwgds_b2_subset_37", readonly = FALSE) # This would also have 3861 individuals

# FLARE scan annotation 
gdsflare7_b2 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b2_all.gds")

original_flareid <- read.gdsn(index.gdsn(gdsflare7_b2, "sample.id"))
flareid2 <- as.data.frame(original_flareid)
colnames(flareid2) <- "SUBJECT_ID"
covariate_b2$SUBJECT_ID <- gsub("SoL", "", covariate_b2$SUBJECT_ID)
flareid2 <- flareid2 |> inner_join(covariate_b2,
                                   by = c("SUBJECT_ID")) |>
  select(SUBJECT_ID)
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=original_scanid, stringsAsFactors=FALSE))

scanAnnot_flare2 <- ScanAnnotationDataFrame(data.frame(
  scanID=read.gdsn(index.gdsn(gdsflare7_b2, "sample.id")), stringsAsFactors=FALSE)) 

# RFMix admixture mapping function
runassoc <- function(ancestry, gdsfile, nullmod, inference = "old"){
  gds.reader <- GdsGenotypeReader(gdsfile, 
                                  genotypeVar=paste0("dosage_", ancestry))
  
  genodata <- GenotypeData(gds.reader, scanAnnot=scanAnnot)
  iterator <- GenotypeBlockIterator(genodata, snpBlock = 100)
  geno <- getGenotypeSelection(iterator)
  colnames(geno) <- original_scanid
  geno <- geno[, colnames(geno) %in% tmp2]
  testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
  # this one uses Score test 
  res <- testres
  
  # loop through all the snp blocks
  while(iterateFilter(iterator)){
    geno <- getGenotypeSelection(iterator)
    colnames(geno) <- original_scanid
    geno <- geno[, colnames(geno) %in% tmp2]
    # table(geno) only gives 0, 1, 2
    testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
    res <- rbind(testres, res)
  }
  return(res)
}

# flare3 admixture mapping function 
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare3 <- function(ancestry, nullmod, metab, from, to){
  for(i in seq(from=from, to = to)){
    print(paste0("working on chromosome", i))
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    gdsflare_b2 <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap/flare3_b2_chr", 
                                     i, ".gds" ))
    gds <- GdsGenotypeReader(gdsflare_b2,  genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare2)
    chrom <- getChromosome(gds)[1:length(filtered_id)]
    snppos <- as.data.frame(getPosition(gds))
    rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))
    snppos <- snppos[rownames(snppos) %in% filtered_id,]
    geno <- getGenotype(genodata)
    rownames(geno) <- read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))
    colnames(geno) <- original_flareid
    geno <- geno[rownames(geno) %in% filtered_id, colnames(geno) %in% flareid2$SUBJECT_ID]
    print(dim(geno))
    closefn.gds(gdsflare_b2)
    closefn.gds(filtered_gds)
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_low_corr/", metab, "FLARE3_", 
                            ancestry ,"_b2_fixed_chr", i, ".Rdata"))
    gc()
  }
}


# flare7 admixture mapping function 
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare7 <- function(ancestry, from, to, nullmod, metab){
  for(i in seq(from=from, to = to) ){
    print(paste0("working on chromosome", i))
    gdsflare_b2 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b2_all.gds")
    gds <- GdsGenotypeReader(gdsflare_b2,  genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare2)
    chrom <- getChromosome(gds)
    chromid <- which(chrom == i)
    snppos <- as.data.frame(getPosition(gds)[chromid])
    
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    
    chrom <- getChromosome(gds)[1:length(filtered_id)]
    rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))[chromid]
    snppos <- snppos[rownames(snppos) %in% filtered_id,]
    
    geno <- getGenotype(genodata, snp = c(min(chromid), 
                                          max(chromid) - min(chromid) +1))
    rownames(geno) <- read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))[chromid]
    colnames(geno) <- original_flareid
    geno <- geno[rownames(geno) %in% filtered_id, colnames(geno) %in% flareid2$SUBJECT_ID]
    closefn.gds(gdsflare_b2)
    closefn.gds(filtered_gds)
    
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_low_corr/", metab, "FLARE7_", 
                            ancestry ,"_b2_fixed_chr", i, ".Rdata"))
    gc()
  }
}

# Run admixture mapping
# Metabolite 4-hydroxybenzoate sulfate - from chromosome 16, 100006264
res_afr_x6264 <-  runassoc("afr", gds_old_b2, nullmod_b2_x6264)
res_amer_x6264 <-  runassoc("amer", gds_old_b2, nullmod_b2_x6264)

# Using FLARE3 local ancestry inference
runassoc_flare3("afr", from = 16, to = 16, nullmod = nullmod_b2_x6264,
                metab = "x6264")
runassoc_flare3("amer", from = 16, to = 16, nullmod = nullmod_b2_x6264,
                metab = "x6264")

# Using FLARE7 local ancestry inference 
runassoc_flare7("afr", from = 16, to = 16, nullmod = nullmod_b2_x6264,
                metab = "x6264")
runassoc_flare7("amer", from = 16, to = 16, nullmod = nullmod_b2_x6264,
                metab = "x6264")

# Save RFMix admixture mapping results from batch2:
saveRDS(res_afr_x6264, file = "./admix_map_all/admixmap_low_corr/RFMix_res_afr_x6264_b2_fixed.rds")
saveRDS(res_amer_x6264, file = "./admix_map_all/admixmap_low_corr/RFMix_res_amer_x6264_b2_fixed.rds")
