# load neccessary packages 
library(tidyverse)
library(GENESIS)
library(gdsfmt)
library(GWASTools)
library(SNPRelate)

# Load the old-RFMix GDS file
uwgds <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries\\lai_HGDP_1000G_comb_unique.gds"
#uwgds <-"/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries/lai_HGDP_1000G_comb_unique.gds"
gdsold <- openfn.gds(uwgds)

# Load the FLARE7 GDS file
flare <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7\\ancestry_counts_FLARE.gds"
gdsflare <- openfn.gds(flare)

covariate_b2_path <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\metab_admixture_association_batch2.csv"
#covariate_b2_path <- "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/metab_admixture_association_batch2.csv"

covariate_b2 <- read.csv(covariate_b2_path)

# load the kinship matrix -- already generated previously 
#kinpath <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files"
kinpath <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files"
kin.mat <- getobj(file.path(kinpath, "pcreap_phi_matrix.RData"))

# load the household matrices -- generated from 20240418prep_data_metab_admix_map.R
# Household matrix for batch2  
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\hh.matrix_b2.Rdata")

# load block matrices -- generated from 20240418prep_data_metab_admix_map.R
# Block matrices for batch2 
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\block.matrix_b2.Rdata")

# Load SNP Annotations 
snpannot <-  "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries\\lai_HGDP_1000G_comb_snpAnnot_unique.RData"
# snpannot <-  "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries/lai_HGDP_1000G_comb_snpAnnot_unique.RData"
load(snpannot)

# Load ID information for FLARE batch 2
load("./admix_map_all/flare_batch2_SoLids.Rdata")

################## BATCH 2 -  Replication batch #######################

# subset gds file for batch 2
sum(is.na(covariate_b2))
#sum(is.na(covariate_b2$X1114))
covariate_b2 <- covariate_b2[complete.cases(covariate_b2),]
# 1836 complete cases for batch 2

# subset gds file for batch 2
gdsSubset(uwgds , "./uwgds_b2",
          sample.include=covariate_b2$GAC_scanID,
          sub.storage=NULL)
gdsSubsetCheck(uwgds, "./uwgds_b2", sample.include=covariate_b2$GAC_scanID)
gds_old_b2 <- openfn.gds("./uwgds_b2")

tmp2 <- read.gdsn(index.gdsn(gds_old_b2, "sample.id")) 
# the subject order of the subset gds file for batch2

# get rid of the 3s in the local ancestry counts 
# first identify which individuals have ancestry counts = 3
# Check 3s in the local ancestry counts -- afr
# afr <- getGenotype(genoDataList[["afr"]]) #2470 NAs, which are 3s
afr <- t(as.data.frame(read.gdsn(index.gdsn(gds_old_b2, "dosage_afr")))) 
# sum(is.na(afr)) == 0
colnames(afr) <- tmp2
columns_to_remove <- apply(afr, 2, function(col) any(col == 3))
# Remove these columns
afr <- afr[, !columns_to_remove] # removed 2 individuals
# table(afr) check using this

# resubsetting the gds file, matricies, and covariates 
covariate_b2 <- covariate_b2 |> filter(GAC_scanID %in% colnames(afr))
gds_old_b2 <- openfn.gds("./uwgds_b2_subset_37") #load subsetted RFMix GDS for batch 2

# First get the sample ids from the FLARE file to find overlappings
flare_ids <- as.data.frame(read.gdsn(index.gdsn(gdsflare, "sample.id"))) 
colnames(flare_ids) <- "sample.id"
covariate_b2$SUBJECT_ID <- gsub(".*SoL", "",  covariate_b2$SUBJECT_ID )
flare_b2 <- flare_ids$sample.id[match(covariate_b2$SUBJECT_ID, flare_ids$sample.id)]
flare_b2 <- flare_b2[!is.na(flare_b2)] # 1644 individuals found, removed 192 NAs
covariate_b2 <- covariate_b2 |> filter(SUBJECT_ID %in% flare_b2) # ONLY take the overlapping individuals

flare_b2 <- as.data.frame(flare_b2)
save(flare_b2, file = "flare_batch2_SoLids.Rdata") # saving the SoL IDs for batch 2

# Load ID information for FLARE batch 2
load("./admix_map_all/flare_batch2_SoLids.Rdata")

# match individuals that are in batch 2, with correct order based on the GDS file
tmp2 <- read.gdsn(index.gdsn(gds_old_b2, "sample.id")) # the order of the subset gds file for batch2

# Sort by subject ID now
gdsflare_b2 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b2_all.gds")
flareid <- read.gdsn(index.gdsn(gdsflare_b2, "sample.id"))
covariate_b2 <- covariate_b2[match(as.character(flareid), covariate_b2$SUBJECT_ID),]

# match individuals that are in the filtered batch 2, with correct order based on the GDS file
# the order of the samples is ascending by SoLID
tmp2 <- covariate_b2$GAC_scanID
write.csv(tmp2, file = "./GAC_scanID_order_b2_admixmap.csv")

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
x_b2 <- covariate_b2[as.character(tmp2),c(seq(9,14), seq(16, 19))] # make sure the 
# order of the individuals are the same as in the GDS file

# reindexing the old GDS file
# Load the liftOver BED file for the old local ancestry interval annotations
# Overlifted from hg37 to hg38
intervals_liftover_nodup <- read.table("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\intervals_37-38.bed")

colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 
intervals_liftover_nodup <- intervals_liftover_nodup |> filter(chr != "chrX") # 14753 mapped blocks

# oldind is the  original sample order from both the gdsold files
# in ascending order 
oldind <- as.data.frame(read.gdsn(index.gdsn(gdsold, "sample.id")))
oldind <- oldind[oldind[,1] %in% tmp2,] # keep the ones in the tmp2, still ascending
# 
# gdsSubset(uwgds, "./uwgds_b2_subset_37",
#           sample.include=covariate_b2$GAC_scanID,
#           sub.storage=NULL)
# 
# # Add neccessary information the subsetted gdsfile specifically for batch2 RFMix admixture mapping 
# gds_old_b2 <- openfn.gds("./uwgds_b2_subset_37", readonly = FALSE) # This would also have 1644 individuals
# # sum(read.gdsn(index.gdsn(gds_old_b2, "sample.id")) == oldind) = 1644
# add.gdsn(gds_old_b2, name = "sample.id", val = tmp2, replace = TRUE, 
#          compress = "LZMA_ra")
# add.gdsn(gds_old_b2, name = "snp.id", val = intervals_liftover_nodup$snpid, 
#          replace = TRUE, compress = "LZMA_ra")
# val = read.gdsn(index.gdsn(gds_old_b2, "snp.position"))[intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b2, name = "snp.position", val = val, replace = TRUE, 
#          compress = "LZMA_ra")
# 
# val = read.gdsn(index.gdsn(gds_old_b2, "snp.chromosome"))[intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b2, name = "snp.chromosome", val = val, replace = TRUE, 
#          compress = "LZMA_ra")
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b2, "dosage_afr")))
# rownames(val) <- oldind
# val <- val[as.character(tmp2), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b2, name = "dosage_afr", val = as.matrix(val), replace = TRUE, 
#          compress = "LZMA_ra")
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b2, "dosage_amer")))
# rownames(val) <- oldind
# val <- val[as.character(tmp2), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b2, name = "dosage_amer", val = as.matrix(val), replace = TRUE, 
#          compress = "LZMA_ra")
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b2, "dosage_eur")))
# rownames(val) <- oldind
# val <- val[as.character(tmp2), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b2, name = "dosage_eur",val = as.matrix(val), replace = TRUE, 
#          compress = "LZMA_ra")

########## Run admixture mapping, focusing on African and Amerindian ancestry ########

# first fit the models under the null hypothesis of no genetic ancestry effect 
# while including multiple random (pairwise kinship coefficients, household, and census block group) 
# and fixed (age, sex, eGFR, recruitment center, genetic analysis group, and the first five PCs) effects.

# define the nullmod function
nullmod <- function(outcome, covdat, kinmatrix, hhmatrix, blockmatrix){
  covMatList <- list(HH = hhmatrix, kinship = kinmatrix, block = blockmatrix)
  x <- cbind(covariate_b1[as.character(tmp), outcome], covdat)
  colnames(x)[1] <- outcome
  mod <- fitNullModel(x = x, outcome = outcome, 
                      covars=c("AGE","GENDER","EV1","EV2","EV3","EV4","EV5","gengrp6", "GFRSCYS", "CENTER"),
                      cov.mat = covMatList, 
                      verbose = TRUE)
  return(mod)
}

# Metabolite 3-aminoisobutyrate
nullmod_b2_x1114 <- nullmod("X1114", x_b2, kin.mat_b2, hh.matrix_b2, 
                            block.matrix_b2)

# Metabolite N-acetylarginine
nullmod_b2_x1266 <- nullmod("X100001266", x_b2, kin.mat_b2, hh.matrix_b2, 
                            block.matrix_b2)

# Metabolite PE 16:0/20:4
nullmod_b2_x8990 <- nullmod("X100008990", x_b2, kin.mat_b2, hh.matrix_b2, 
                            block.matrix_b2)

# Metabolite PC 16:0/20:4
nullmod_b2_x8914 <- nullmod("X100008914", x_b2, kin.mat_b2, hh.matrix_b2, 
                            block.matrix_b2)

# Saving all the null models
save(nullmod_b2_x1114, file = "nullmod_b2_x1114_fixed.Rdata")
save(nullmod_b2_x1266, file = "nullmod_b2_x1266_fixed.Rdata")
save(nullmod_b2_x18990, file = "nullmod_b2_x18990_fixed.Rdata")
save(nullmod_b2_x18914, file = "nullmod_b2_x18914_fixed.Rdata")

runassoc <- function(ancestry, gdsfile, nullmod, inference = "old"){
  gds.reader <- GdsGenotypeReader(gdsfile, 
                                    genotypeVar=paste0("dosage_", ancestry))
  
  genodata <- GenotypeData(gds.reader, scanAnnot=scanAnnot)
  
  iterator <- GenotypeBlockIterator(genodata, snpBlock = 100)
  geno <- getGenotypeSelection(iterator)
  testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) # this one uses Score test 
  res <- testres
  
  # loop through all the snp blocks
  while(iterateFilter(iterator)){
    geno <- getGenotypeSelection(iterator)
    # table(geno) only gives 0, 1, 2
    testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
    res <- rbind(testres, res)
  }
  return(res)
}


# ancestry takes abbreviations such as "afr", "amer", or "eur"
# Returns admixture mapping results for each chromosome and each ancestry in a 
# .Rdata format stored in the local directory
# Parameter ancestry: the ancestry to be examined at & extracting local ancestry counts from
# from: starting chromosome to do admixture mapping in a for-loop
# to: ending chromosome to do admixture mapping in a for-loop
# nullmod: the nulllmod to be used in admixture mapping computed above 
# metab: the metabolites (outcome) to be tested for
runassoc_flare <- function(ancestry, from = 22, to = 1, nullmod, metab){
  for(i in seq(from=from, to = to) ){
    print(paste0("working on chromosome", i))
    gdsflare_b2 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b2_all.gds")
    gds <- GdsGenotypeReader(gdsflare_b2,  genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
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
    geno <- geno[rownames(geno) %in% filtered_id,]
    closefn.gds(gdsflare_b2)
    closefn.gds(filtered_gds)
    
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_", metab, 
                            "_all_FLARE7/admix", 
                            ancestry ,"_b2_filtered_chr", i, ".Rdata"))
    gc()
  }
}

# Loading/creating Scan annotation objects for both the old and new inference
# RFMix inference scan annotation
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=tmp2, stringsAsFactors=FALSE))

# FLARE7 scan annotation 
scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
  scanID=flareid, stringsAsFactors=FALSE)) # flareid from gdsflare_b2, smaple.id

# Save the scanAnnotation object
save(scanAnnot_flare, file = "scanAnnot_b2_admixmap.Rdata")

################# Run admixture mapping on 3-aminoisobutyrate ##################

res_afr_x1114_b2 <-  runassoc("afr", gds_old_b2, 
                              nullmod_b2_x1114) #RFMix, afr ancestry
res_amer_x1114_b2 <- runassoc("amer", gds_old_b2, 
                              nullmod_b2_x1114) #RFMix, amer ancestry

runassoc_flare("afr", nullmod = nullmod_b2_x1114, 
               metab ="x1114" ) #FLARE7, afr ancestry
runassoc_flare("amer", nullmod = nullmod_b2_x1114,
               metab = "x1114") #FLARE7, amer ancestry

################# Run admixture mapping on metabolite N-acetylarginine #########

res_afr_x1266_b2 <-  runassoc("afr", gds_old_b2, 
                              nullmod_b2_x1266) #RFMix, afr ancestry
res_amer_x1266_b2 <- runassoc("amer", gds_old_b2, 
                              nullmod_b2_x1266) #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b2_x1266, 
               metab = "x1266") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b2_x1266,
               metab = "x1266") #FLARE7, amer ancestry

################# Run admixture mapping on metabolite PE 16:0/20:4 #############

res_afr_x8990_b2 <-  runassoc("afr", gds_old_b2, 
                              nullmod_b2_x8990)  #RFMix, afr ancestry
res_amer_x8990_b2 <- runassoc("amer", gds_old_b2, 
                              nullmod_b2_x8990) #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b2_x8990, 
               metab = "x8990") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b2_x8990,
               metab = "x8990") #FLARE7, amer ancestry

################# Run admixture mapping on metabolite PC 16:0/20:4 #############

res_afr_x8914_b2 <-  runassoc("afr", gds_old_b2, 
                              nullmod_b2_x8914) #RFMix, afr ancestry
res_amer_x8914_b2 <-  runassoc("amer", gds_old_b2, 
                               nullmod_b2_x8914) #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b2_x8914, 
               metab = "x8914") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b2_x8914, 
               metab = "x8914") #FLARE7, amer ancestry

######### Saving the RFMix Admixture mapping results ###########################
saveRDS(res_afr_x1114_b2, file = "./admix_map_all/RFMix_res_afr_x1114_b2.rds")
saveRDS(res_amer_x1114_b2, file = "./admix_map_all/RFMix_res_amer_x1114_b2.rds")
saveRDS(res_afr_x1266_b2, file = "./admix_map_all/RFMix_res_afr_x1266_b2.rds")
saveRDS(res_amer_x1266_b2, file = "./admix_map_all/RFMix_res_amer_x1266_b2.rds")
saveRDS(res_afr_x8990_b2, file = "./admix_map_all/RFMix_res_afr_x8990_b2.rds")
saveRDS(res_amer_x8990_b2, file = "./admix_map_all/RFMix_res_amer_x8990_b2.rds")
saveRDS(res_afr_x8990_b2, file = "./admix_map_all/RFMix_res_afr_x8990_b2.rds")
saveRDS(res_amer_x8990_b2, file = "./admix_map_all/RFMix_res_amer_x8990_b2.rds")
