library(gdsfmt)
library(tidyverse)

### Create global ancestry data ###

## FLARE7 ##
afr_tot <- c()
amer_tot <- c()
eur_tot <- c()
easia_tot <- c()
csasia_tot <- c()
ocea_tot <- c()
me_tot <- c()

for(i in 22:1){
  print(paste0("working on chromosome", i))
  subsetted_path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr", i, ".gds")
  subset_gds <- openfn.gds(subsetted_path,  readonly  = FALSE)
  afr <- rowSums(as.data.frame(read.gdsn(index.gdsn(subset_gds, "afr_counts"))))
  amer <- rowSums(as.data.frame(read.gdsn(index.gdsn(subset_gds, "amer_counts"))))
  eur <- rowSums(as.data.frame(read.gdsn(index.gdsn(subset_gds, "eur_counts"))))
  easia <- rowSums(as.data.frame(read.gdsn(index.gdsn(subset_gds, "easia_counts"))))
  csasia <- rowSums(read.gdsn(index.gdsn(subset_gds, "csasia_counts")))
  ocea <- rowSums(read.gdsn(index.gdsn(subset_gds, "ocea_counts")))
  me <- rowSums(read.gdsn(index.gdsn(subset_gds, "me_counts")))
  gc()
  if(i == 22){
    afr_tot <- afr
    amer_tot <- amer
    eur_tot  <- eur
    easia_tot <- easia
    csasia_tot <- csasia
    ocea_tot <- ocea
    me_tot <- me
  }
  else{
    afr_tot <- afr_tot + as.data.frame(afr)
    amer_tot <- amer_tot + as.data.frame(amer)
    eur_tot <- eur_tot + as.data.frame(eur)
    easia_tot <- easia_tot + as.data.frame(easia)
    csasia_tot <- csasia_tot + as.data.frame(csasia)
    ocea_tot <- ocea_tot  + as.data.frame(ocea)
    me_tot <- me_tot + as.data.frame(me)
  }
  # gc()
  # print("do gds nodes")
  # add.gdsn(subset_gds, name = "easia_counts", 
  #          val = counts_easia,
  #          compress = "LZMA_RA", replace = TRUE)
  # add.gdsn(subset_gds, name = "csasia_counts", 
  #          val = counts_csasia,
  #          compress = "LZMA_RA", replace = TRUE)
  # add.gdsn(subset_gds, name = "ocea_counts", 
  #          val = counts_ocea,
  #          compress = "LZMA_RA", replace = TRUE)
  # add.gdsn(subset_gds, name = "me_counts", 
  #          val = counts_me,
  #          compress = "LZMA_RA", replace = TRUE)
  # 
  # rm(counts_easia, counts_csasia, counts_ocea, counts_me)
  closefn.gds(subset_gds)
  gc()
  
}

global_tot <- as.data.frame(sapply(c(afr_tot, amer_tot, eur_tot,
                                     easia_tot, csasia_tot, 
                                     ocea_tot, me_tot),  function(x){
                                       x/5091756/2
                                     }))

saveRDS(global_tot, 
        file = "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_global_prop_QC_tot.rds")
