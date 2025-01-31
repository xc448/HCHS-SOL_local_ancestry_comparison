# load neccessary packages 
library(tidyverse)
library(GENESIS)
library(gdsfmt)
library(GWASTools)
library(SNPRelate)

# Sort by subject ID to get ID info
# gdsflare_b2_flare7 <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_b2_all.gds")
# flareid <- read.gdsn(index.gdsn(gdsflare_b2_flare7, "sample.id"))

# Load ID information for FLARE batch 2
load("./admix_map_all/flare_batch2_SoLids.Rdata")

################ Admixture mapping analysis starts  here ####################
#  Load null models 
load("nullmod_b2_x1114_fixed.Rdata")
load("nullmod_b2_x1266_fixed.Rdata")
load("nullmod_b2_x8990_fixed.Rdata")
load("nullmod_b2_x8914_fixed.Rdata")

# FLARE scan annotation 
scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
  scanID=flareid, stringsAsFactors=FALSE)) 

# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare <- function(ancestry, nullmod, metab, from = 22, to = 1){
  for(i in from:to){
    print(paste0("working on chromosome", i))
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    gdsflare_b2 <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap/flare3_b2_chr", 
                                     i, ".gds" ))
    gds <- GdsGenotypeReader(gdsflare_b2,  
                             genotypeVar=paste0(ancestry, "_counts"))
    genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
    chrom <- getChromosome(gds)[1:length(filtered_id)]
    snppos <- as.data.frame(getPosition(gds))
    rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))
    snppos <- snppos[rownames(snppos) %in% filtered_id,]
    geno <- getGenotype(genodata)
    rownames(geno) <- read.gdsn(index.gdsn(gdsflare_b2, "snp.id"))
    geno <- geno[rownames(geno) %in% filtered_id,]
    closefn.gds(gdsflare_b2)
    closefn.gds(filtered_gds)
    res <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
    res <- cbind(snppos, res)
    save(res, file = paste0("./admix_map_all/admixmap_", metab, "_all_FLARE3/admix", 
                            ancestry ,"_b2_filtered_fixed_chr", i, ".Rdata"))
    gc()
  }
}

#  x1114 -- 3-aminoisobutyrate
runassoc_flare("afr", nullmod = nullmod_b2_x1114, metab = "x1114")
runassoc_flare("amer", nullmod = nullmod_b2_x1114, metab = "x1114")

#  x1266 -- N-acetylarginine
runassoc_flare("afr", nullmod = nullmod_b2_x1266, metab = "x1266")
runassoc_flare("amer", nullmod = nullmod_b2_x1266, metab = "x1266")

#  x8990 -- PE 16:0/20:4
runassoc_flare("afr", nullmod = nullmod_b2_x8990, metab = "x8990")
runassoc_flare("amer", nullmod = nullmod_b2_x8990, metab = "x8990")

#  x8914 -- PC 16:0/20:4
runassoc_flare("afr", nullmod = nullmod_b2_x8914, metab = "x8914")
runassoc_flare("amer", nullmod = nullmod_b2_x8914, metab = "x8914")
