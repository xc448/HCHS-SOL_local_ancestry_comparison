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

# Load the covariate dataset -- output from 20240418prep_data_metab_admix_map.R
covariate_b1_path <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\metab_admixture_association_batch1.csv"
#covariate_b1_path <- "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/metab_admixture_association_batch1.csv"

covariate_b1 <- read.csv(covariate_b1_path)

# load the kinship matrix -- already generated previously 
#kinpath <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files"
kinpath <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files"

kin.mat <- getobj(file.path(kinpath, "pcreap_phi_matrix.RData"))

# load the household matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch1
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\hh.matrix_b1.Rdata")

# load block matrices -- generated from 20240418prep_data_metab_admix_map.R
# batch1
load("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\block.matrix_b1.Rdata")

# Load SNP Annotations 
snpannot <-  "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries\\lai_HGDP_1000G_comb_snpAnnot_unique.RData"
#snpannot <-  "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries/lai_HGDP_1000G_comb_snpAnnot_unique.RData"
load(snpannot)

########### Dataset cleaning/subsetting for batch  1 ##################

# remove NAs in the covariate datasets -- only keep the complete cases
sum(is.na(covariate_b1)) 
covariate_b1 <- covariate_b1[complete.cases(covariate_b1),]
# 97 individuals are dropped

# subset gds file for batch 1
closefn.gds(gdsold)
gdsSubset(uwgds , "./uwgds_b1",
          sample.include=covariate_b1$GAC_scanID,
          sub.storage=NULL)
gdsSubsetCheck(uwgds, "./uwgds_b1", sample.include=covariate_b1$GAC_scanID)

gds_old_b1 <- openfn.gds("./uwgds_b1") 


# Get the sample ID order of the subsetted gds file for batch 1
tmp <- read.gdsn(index.gdsn(gds_old_b1, "sample.id")) 

# get rid of the 3s in the  old local ancestry counts -- Missing values 
# first identify which individuals have ancestry counts == 3
# Check 3s in the local ancestry counts -- afr
#afr <- getGenotype(genoDataList[["afr"]]) #2470 NAs, which are 3s
afr <- t(as.data.frame(read.gdsn(index.gdsn(gds_old_b1, "dosage_afr")))) # sum(is.na(afr)) == 0
colnames(afr) <- tmp
columns_to_remove <- apply(afr, 2, function(col) any(col == 3))

# Remove these columns
afr <- afr[, !columns_to_remove] # removed 4 individual
table(afr) # verify that this only gives 0, 1, and 2

# only subsetting from afr ancestry is sufficient for NA removal
# amer <- t(as.data.frame(read.gdsn(index.gdsn(gds_old_b1, "dosage_amer"))))
# colnames(amer) <- tmp
# columns_to_remove_amer <- apply(amer, 2, function(col) any(col == 3))
# sum(columns_to_remove == columns_to_remove_amer )

## resubsetting the gds files, covariate datasets, and matricies for batch 1
# subsetting the covariate dataset for batch 1
covariate_b1 <- covariate_b1 |> filter(GAC_scanID %in% colnames(afr))

# match individuals in the batch 1 from both old and FLARE gds file

# First get the sample ids from the FLARE file to find overlappings
flare_ids <- as.data.frame(read.gdsn(index.gdsn(gdsflare, "sample.id"))) 
colnames(flare_ids) <- "sample.id"
covariate_b1$SUBJECT_ID <- gsub(".*SoL", "",  covariate_b1$SUBJECT_ID )
flare_b1 <- flare_ids$sample.id[match(covariate_b1$SUBJECT_ID, (flare_ids$sample.id))]
flare_b1 <- flare_b1[!is.na(flare_b1)] # 3861 individuals found, removed 37 NAs
covariate_b1 <- covariate_b1 |> filter(SUBJECT_ID %in% flare_b1) # ONLY take the overlapping individuals

flare_b1 <- as.data.frame(flare_b1)
save(flare_b1, file = "./flare_batch1_SoLids.Rdata") # saving the SoL IDs for batch 1

# loading the subetted flare gds file 
gdsflare_b1 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b1_all.gds")

# Sort by subject ID now
flareid <- read.gdsn(index.gdsn(gdsflare_b1, "sample.id"))
covariate_b1 <- covariate_b1[match(as.character(flareid), covariate_b1$SUBJECT_ID),]

# match individuals that are in the filtered batch 1, with correct order based on the GDS file
# the order of the samples is ascending by SoLID
tmp <- covariate_b1$GAC_scanID
write.csv(tmp, file = "./GAC_scanID_order_b1_admixmap.csv")

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
x_b1 <- covariate_b1[as.character(tmp),c(seq(9,14), seq(16, 19))] # make sure the 
# order of the individuals are the same as in the GDS file

# reindexing the old GDS file
# Load the liftOver BED file for the old local ancestry interval annotations
# Overlifted from hg37 to hg38, nodup meand no duplications when  mapping intervals
intervals_liftover_nodup <- read.table("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\intervals_37-38.bed")

colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 
intervals_liftover_nodup <- intervals_liftover_nodup |> filter(chr != "chrX") # 14506 mapped blocks

gds_old_b1 <- openfn.gds("./uwgds_b1_subset_37", readonly = FALSE) # This would also have 3861 individuals

# # Add neccessary information the subsetted gdsfile specifically for batch1 RFMix admixture mapping 
# add.gdsn(gds_old_b1, name = "sample.id", val = tmp, replace = TRUE, compress = "LZMA_ra")
# add.gdsn(gds_old_b1, name = "snp.id", val = intervals_liftover_nodup$snpid, replace = TRUE, compress = "LZMA_ra")
# val = read.gdsn(index.gdsn(gds_old_b1, "snp.position"))[intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b1, name = "snp.position", val = val, replace = TRUE, compress = "LZMA_ra")
# 
# val = read.gdsn(index.gdsn(gds_old_b1, "snp.chromosome"))[intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b1, name = "snp.chromosome", val = val, replace = TRUE, compress = "LZMA_ra")
# 
# oldind <- as.data.frame(read.gdsn(index.gdsn(gdsold, "sample.id")))
# oldind <-oldind[oldind[,1] %in% tmp,]
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b1, "dosage_afr")))
# rownames(val) <- oldind
# val <- val[as.character(tmp), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b1, name = "dosage_afr", val = as.matrix(val), replace = TRUE, compress = "LZMA_ra")
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b1, "dosage_amer")))
# rownames(val) <- oldind
# val <- val[as.character(tmp), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b1, name = "dosage_amer", val = as.matrix(val), replace = TRUE, compress = "LZMA_ra")
# 
# val = as.data.frame(read.gdsn(index.gdsn(gds_old_b1, "dosage_eur")))
# rownames(val) <- oldind
# val <- val[as.character(tmp), intervals_liftover_nodup$snpid]
# add.gdsn(gds_old_b1, name = "dosage_eur",val = as.matrix(val), replace = TRUE, compress = "LZMA_ra")


########## Run admixture mapping, focusing on African and Amerindian ancestry ########

# first fit the models under the null hypothesis of no genetic ancestry effect 
# while including multiple random (pairwise kinship coefficients, household, and census block group) 
# and fixed (age, sex, eGFR, recruitment center, genetic analysis group, and the first five PCs) effects.
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
nullmod_b1_x1114 <- nullmod("X1114", x_b1, kin.mat_b1, hh.matrix_b1, block.matrix_b1)

# Metabolite N-acetylarginine
nullmod_b1_x1266 <- nullmod("X100001266", x_b1, kin.mat_b1, hh.matrix_b1, block.matrix_b1)

# Metabolite PE 16:0/20:4
nullmod_b1_x8990 <- nullmod("X100008990", x_b1, kin.mat_b1, hh.matrix_b1, block.matrix_b1)

# Metabolite PC 16:0/20:4
nullmod_b1_x8914 <- nullmod("X100008914", x_b1, kin.mat_b1, hh.matrix_b1, block.matrix_b1)

# Saving all the null models
save(nullmod_b1_x1114, file = "nullmod_b1_x1114_fixed.Rdata")
save(nullmod_b1_x1266, file = "nullmod_b1_x1266_fixed.Rdata")
save(nullmod_b1_x8990, file = "nullmod_b1_x8990_fixed.Rdata")
save(nullmod_b1_x8914, file = "nullmod_b1_x8914_fixed.Rdata")

# runassoc function for RFMix gds file
runassoc <- function(ancestry, gdsfile, nullmod, inference = "old"){
  gds.reader <- GdsGenotypeReader(gdsfile, 
                                  genotypeVar=paste0("dosage_", ancestry))

  genodata <- GenotypeData(gds.reader, scanAnnot=scanAnnot)
  
  iterator <- GenotypeBlockIterator(genodata, snpBlock = 100)
  geno <- getGenotypeSelection(iterator)
  testres <- GENESIS:::testGenoSingleVar(nullmod, t(geno)) 
  # this one uses Score test 
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

# Loading/creating Scan annotation objects for both the old and new inference
# RFMix inference scan annotation
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=read.gdsn(index.gdsn(gds_old_b1, "sample.id")), stringsAsFactors=FALSE))
save(scanAnnot, file = "scanAnnot_b1")

# FLARE7 scan annotation 
scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
  scanID=flareid, stringsAsFactors=FALSE))
# Save the scanAnnotation object
save(scanAnnot_flare, file = "scanAnnot_b1_admixmap.Rdata")
# load("scanAnnot_b1_admixmap.Rdata") # the flare annotation data

# ancestry takes abbreviations such as "afr", "amer", or "eur"
# Returns admixture mapping results for each chromosome and each ancestry in a 
# .Rdata format stored in the local directory
# Parameter ancestry: the ancestry to be examined at & extracting local ancestry counts from
# from: starting chromosome to do admixture mapping in a for-loop
# to: ending chromosome to do admixture mapping in a for-loop
# nullmod: the nulllmod to be used in admixture mapping computed above 
# metab: the metabolites (outcome) to be tested for
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare <- function(ancestry, from, to, nullmod, metab){
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
    closefn.gds(gdsflare_b1)
    closefn.gds(filtered_gds)
    
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_", metab, "_all_FLARE7/admix", 
                            ancestry ,"_b1_filtered_fixed_chr", i, ".Rdata"))
    gc()
  }
}

################# Run admixture mapping on 3-aminoisobutyrate ##################

res_afr_x1114 <-  runassoc("afr", gds_old_b1, 
                           nullmod_b1_x1114) #RFMix, afr ancestry
res_amer_x1114 <- runassoc("amer", gds_old_b1, 
                           nullmod_b1_x1114)  #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b1_x1114,
               metab = "x1114") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b1_x1114,
               metab = "x1114") #FLARE7, amer ancestry 


################# Run admixture mapping on metabolite N-acetylarginine #########
res_afr_x1266 <-  runassoc("afr", gds_old_b1, 
                           nullmod_b1_x1266) #RFMix, afr ancestry
res_amer_x1266 <- runassoc("amer", gds_old_b1, 
                           nullmod_b1_x1266) #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b1_x1266, 
               metab = "x1266") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b1_x1266,
               metab = "x1266") #FLARE7, amer ancestry


################# Run admixture mapping on metabolite PE 16:0/20:4 #############

res_afr_x8990 <-  runassoc("afr", gds_old_b1, 
                           nullmod_b1_x8990) #RFMix, afr ancestry
res_amer_x8990 <- runassoc("amer", gds_old_b1, 
                           nullmod_b1_x8990) #RFMix, amer ancestry

runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b1_x8990,
               metab = "x8990") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b1_x8990,
               metab = "x8990") #FLARE7, amer ancestry


################# Run admixture mapping on metabolite PC 16:0/20:4 #############

res_afr_x8914 <-  runassoc("afr", gds_old_b1, 
                           nullmod_b1_x8914) #RFMix, afr ancestry
res_amer_x8914 <- runassoc("amer", gds_old_b1, 
                           nullmod_b1_x8914) #RFMix, amer ancestry


runassoc_flare("afr", from = 22, to = 1, nullmod = nullmod_b1_x8914,
               metab = "x8914") #FLARE7, afr ancestry
runassoc_flare("amer", from = 22, to = 1, nullmod = nullmod_b1_x8914,
               metab = "x8914") #FLARE7, amer ancestry

###########Save the old fitted results from the old (RFMix) inference  #########

save(res_afr_x1114, file  =  "RFMix_res_afr_x1114_b1_fixed.Rdata")
save(res_amer_x1114, file  =  "RFMix_res_amer_x1114_b1_fixed.Rdata")
save(res_afr_x1266, file  =  "RFMix_res_afr_x1266_b1_fixed.Rdata")
save(res_amer_x1266, file  =  "RFMix_res_amer_x1266_b1_fixed.Rdata")
save(res_afr_x8914, file  =  "RFMix_res_afr_x8914_b1_fixed.Rdata")
save(res_amer_x8914, file  =  "RFMix_res_amer_x8914_b1_fixed.Rdata")
save(res_afr_x8990, file  =  "RFMix_res_afr_x8990_b1_fixed.Rdata")
save(res_amer_x8990, file  =  "RFMix_res_amer_x8990_b1_fixed.Rdata")