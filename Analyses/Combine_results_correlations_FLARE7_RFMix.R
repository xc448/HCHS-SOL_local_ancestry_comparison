# Combining Rdata files downloaded from the cluster
library(tidyverse)

################ Combine Rdata files from AWS ###########################

# combine data from the first 12 chromosomes by ancestry -- The first 12 chromosomes were run separately by ancestry as in FLARE3
# Afr ancestry 
dat_afr <- c()
for(i in seq(1, 12)){
  load(paste0("./AWS/FLARE7_RFMix_correlation_data/corr_old_flare7_afr_chr",i, ".Rdata"))
  dat_afr <- rbind(dat_afr ,corr_result)
}

# Amer ancestry 
dat_amer <- c()
for(i in seq(1, 12)){
  load(paste0("./AWS/FLARE7_RFMix_correlation_data/corr_old_flare7_amer_chr",i, ".Rdata"))
  dat_amer <- rbind(dat_amer ,corr_result)
}

# Eur ancestry
dat_eur <- c()
for(i in seq(1, 12)){
  load(paste0("./AWS/FLARE7_RFMix_correlation_data/corr_old_flare7_eur_chr",i, ".Rdata"))
  dat_eur <- rbind(dat_eur ,corr_result)
}


# merge all ancestries to together. 
dat_tot <-  dat_afr |> inner_join(dat_amer, by = (c("Chromosome" = "Chromosome", "Position" = "Position")))
dat_tot <-  dat_tot |> inner_join(dat_eur, by = (c("Chromosome" = "Chromosome", "Position" = "Position")))

# Reordering the colulumns
dat_tot <-  dat_tot |> dplyr::select(Chromosome, Position, Corr_afr_p,
                                     Corr_afr_s, Corr_amer_p, Corr_amer_s,
                                     Corr_eur_p, Corr_eur_s, Prop_afr_old, 
                                     Prop_afr_new, Prop_amer_old,
                                     Prop_amer_new, Prop_eur_old, Prop_eur_new,
                                     Var_afr_old, Var_afr_new, Var_amer_old,
                                     Var_amer_new, Var_eur_old, Var_eur_new)


# Chromosome 13-22 - ALL ancestries -- this was done with three ancestries all together.
dat_all <- c()
for(i in seq(13, 22)){
  load(paste0("./AWS/FLARE7_RFMix_correlation_data/corr_old_flare7_chr",i, ".Rdata"))
  dat_all <- rbind(dat_all, corr_result)
}

# Combine data from all chromosomes 1-22
dat_tot <- rbind(dat_tot, dat_all)
dim(dat_tot) # 4704075 out of 5105005 -- 92.15% SNPs mapped to the old intervals

# save dataframe
save(dat_tot, file = "correlation_old_flare7_all_SNPs.Rdata")

  