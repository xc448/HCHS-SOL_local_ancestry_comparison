library(gdsfmt)
library(tidyverse)

### Create global ancestry data ###

## FLARE3 ##
# load Gds file
afr_tot <- c()
amer_tot <- c()
eur_tot <- c()

for(i in 22:1){
  print(paste0("working on chr ", i))
  gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr", i, ".gds"))
  afr <- rowSums(as.data.frame(read.gdsn(index.gdsn(gds, "afr_counts"))))
  amer <- rowSums(as.data.frame(read.gdsn(index.gdsn(gds, "amer_counts"))))
  eur <- rowSums(as.data.frame(read.gdsn(index.gdsn(gds, "eur_counts"))))
  closefn.gds(gds)
  if(i == 22){
    afr_tot <- afr
    amer_tot <- amer
    eur_tot  <- eur
  }
  else{
    afr_tot <- afr_tot + as.data.frame(afr)
    amer_tot <- amer_tot + as.data.frame(amer)
    eur_tot <- eur_tot + as.data.frame(eur)
  }
}

global_tot <- as.data.frame(sapply(c(afr_tot, amer_tot, eur_tot),  function(x){
  x/3610937/2
}))

# put the rownames as the subject SOL ID
gds <-  openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE7/")
solid <- read.gdsn(index.gdsn(gds, "sample.id"))
rownames(global_tot) <-  solid
closefn.gds(gds)

# save the global proportion file as RDS
saveRDS(global_tot, file  =  "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_global_prop_QC_tot.rds")

# global_tot <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_global_prop_QC_tot.rds")
