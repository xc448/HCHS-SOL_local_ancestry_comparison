# Combining Rdata files run in the Correlation_ancestry_mean_var_proportion_FLARE3_RFMix file, which was done
# by ancestry by chromosome. 
library(tidyverse)

################ Combine Rdata files from AWS ###########################

# combine data from the 22 chromosomes by ancestry 

# Afr ancestry 
dat_afr <- c()
for(i in seq(1, 22)){
  load(paste0("./AWS/FLARE3_RFMix_correlation_data/corr_old_flare3_afr_chr",i, ".Rdata"))
  dat_afr <- rbind(dat_afr ,corr_result)
}

# Amer ancestry 
dat_amer <- c()
for(i in seq(1, 22)){
  load(paste0("./AWS/FLARE3_RFMix_correlation_data/corr_old_flare3_amer_chr",i, ".Rdata"))
  dat_amer <- rbind(dat_amer ,corr_result)
}

# Eur ancestry
dat_eur <- c()
for(i in seq(1, 22)){
  load(paste0("./AWS/FLARE3_RFMix_correlation_data/corr_old_flare3_eur_chr",i, ".Rdata"))
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

# save dataframe
save(dat_tot, file = "correlation_old_flare3_all_SNPs.Rdata")
