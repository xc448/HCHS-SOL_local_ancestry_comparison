library(gdsfmt)
library(STEAM)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(GWASTools)
library(data.table)

############### RFMix ##################
# Load genetic map of RFMix on hg38 build
RFMix_map <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/RFMix_gds_by_chr/genetic_dist_RFMix.rds")

# Load global ancestry information from RFMix(old)
global_rfmix <- read.csv("R:\\Sofer Lab\\HCHS_SOL/Projects\\2023_local_ancestry_comparison_sol/Data/global_proportion_old.csv")

# Load individual SoLIDs from FLARE dataset
load("./admix_map_all/flare_batch1_SoLids.Rdata")
flare_b1$flare_b1 <- paste0("SoL", flare_b1$flare_b1)
global_prop_b1  <- global_rfmix[global_rfmix$SUBJECT_ID %in% flare_b1$flare_b1,]
# flare_b2$flare_b2 <- paste0("SoL", flare_b2$flare_b2)
# global_prop_b2  <- global_rfmix[global_rfmix$SUBJECT_ID %in% flare_b2$flare_b2,]

# Initialize, only compute threshold on b1
corr.list_rfmix_b1 <- list()

# binsize = 0.5 as default 
binsize <- 0.5

# load subsampled local ancestry correlation list for each chromosome and append into the list
for(i in 1:22){
  load(paste0("../Code/AWS/RFMix_corr_list_STEAM/RFMix_b1_chr_", i, ".Rdata"))
  corr.list_rfmix_b1[[i]] <- corr
}

# combine snp correlations into one data
corr_K3_rfmix <- combine_corr_chr(corr.list_rfmix_b1)
saveRDS(corr.list_rfmix_b1, file = "corr.list_b1_RFMix.rds")
saveRDS(corr_K3_rfmix, file = "corr_K3_b1_RFMix.rds")
get_g(corr_K3_rfmix) # 10.26769 -- number of generation estiamted since admixture

colnames(global_prop_b1)[2:4] <- c("afr",  "eur", "nam")
global_prop_b1 <- global_prop_b1 |>
  dplyr::select(afr, eur, nam)

RFMix_map  <- RFMix_map|>
  select(cM, chr) |>
  arrange(chr)

# Computing threshold 
dlt <- c(0,STEAM:::get_deltas(RFMix_map))

beta <- 0.01*get_g(corr_K3_rfmix)
a <- exp(-beta*dlt)
b <- sqrt(1-exp(-2*beta*dlt)) 
b <- if_else(is.na(b),  0 , b) # handling na 
ab <- list(a=a,b=b)

avg_props <- apply(global_prop_b1[1:3], 2, mean, na.rm=T)
L <- STEAM:::get_L(avg_props)

set.seed(42) # set seed for replication 
# replicate 10000 times, same as FLARE
max_stats <- replicate(10000, STEAMcpp:::simstatSingle(m = nrow(RFMix_map), K = 3, 
                                                      as = ab$a, bs = ab$b, L = L))
# get upper alpha quantile
zstar <- STEAM:::upper_alpha(max_stats, 0.05)
# bootstrapping, from source code, same as FLARE
z_ci <- quantile(bootstrap::bootstrap(max_stats, nboot = 5000, theta = STEAM:::upper_alpha, 0.05)$thetastar, c(0.025,0.975))


thresh <- 2 * pnorm(zstar, lower.tail = F) # 2.045687e-06
thresh_ci <- 2 * pnorm(z_ci, lower.tail = F) # 2.336054e-06 1.832638e-06    

############# FLARE3 ##################

# load genetic distance map from FLARE3 
FLARE3_map <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/genetic_dist_FLARE3.rds")
FLARE3_map <- FLARE3_map |> 
  select(cM, chr)

# set up list to store results
load("./admix_map_all/flare_batch1_SoLids.Rdata")
global_tot <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_global_prop_QC_tot.rds")
flare3_global_b1 <- global_tot[rownames(global_tot)  %in% flare_b1$flare_b1,]
#flare3_global_b2 <- global_tot[rownames(global_tot)  %in% flare_b2$flare_b2,]

corr.list_flare3 <- list()

# load subsampled local ancestry correlation list for each chromosome and append into the list
for(i in 22:1){
  # set up file names
  load(paste0("R:/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/FLARE3_b1_chr_", i, ".Rdata"))
  corr.list_flare3[[i]] <- snps.dt_corr
  # print update
  cat('done with chromosome', i, '\n')
  
}
saveRDS(corr.list_flare3,  
        file = "R:/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE3_gds_steam/corr.list_b1_FLARE3.rds" )

binsize <- 0.5
# computing variables needed for threshold estimation
for(i in 22:1){
  print(i)
  map.df <- readRDS(paste0("./AWS/FLARE3_corr_list_STEAM/FLARE3_map_chr",  i, ".rds"))
  map.df$snp.id <- rownames(map.df)
  corr.list_flare3[[i]]$cM = get_dist(corr.list_flare3[[i]]$snp1, corr.list_flare3[[i]]$snp2, map.df)
  corr.list_flare3[[i]]$theta = L_to_theta(corr.list_flare3[[i]]$cM)
  corr.list_flare3[[i]]$bin = round(corr.list_flare3[[i]]$cM / binsize, 0) * binsize
}

colnames(flare3_global_b1)[1:3] <- c("afr", "nam", "eur")
flare3_global_b1 <- flare3_global_b1 |> dplyr::select(afr, eur, nam)

corr_K3_FLARE3 <- combine_corr_chr(corr.list_flare3)
saveRDS(corr.list_flare3, file = "corr.list_b1_FLARE3.rds")
saveRDS(corr_K3_FLARE3, file = "corr_K3_b1_FLARE3.rds")
get_g(corr_K3_FLARE3) # 10.21922 -- number of generation estiamted since admixture


set.seed(1)
thresh_FLARE3_b1 <- get_thresh_simstat(g = get_g(corr_K3_FLARE3),
                                       map = FLARE3_map, 
                                       props = flare3_global_b1, 
                                       nreps = 10000)

############# FLARE7 ##################
# load FLARE7 global ancestry and genetic map
# set up list to store results
global_tot_flare7 <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_global_prop_QC_tot.rds")
flare7gds <- openfn.gds("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE7/ancestry_counts_FLARE7_QC.gds")
FLARE7_sample_id <- read.gdsn(index.gdsn(flare7gds, "sample.id"))
rownames(global_tot_flare7) = paste0("SoL", FLARE7_sample_id)
flare_b1$flare_b1 <- paste0("SoL", flare_b1$flare_b1)

# get global ancestry, match individuals by SoL ID
flare7_global_b1 <- global_tot_flare7[rownames(global_tot_flare7)  %in% flare_b1$flare_b1,]

# combine the other 4 ancestries: easia, csasia, me, and ocea into one
other4 <- flare7_global_b1 |> 
  select(easia, csasia, ocea, me) |>
  mutate(other4 = rowSums(across(everything()))) |>
  dplyr::select(other4) 

flare7_global_b1$other4 <- other4
flare7_global_b1 <- flare7_global_b1 |> select(afr, amer, eur, other4)

# load genetic distance map from FLARE7 
FLARE7_map <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects/2023_local_ancestry_comparison_sol/Data/FLARE7_gds_steam/genetic_dist_FLARE7.rds")
FLARE7_map <- FLARE7_map |> 
  select(cM, chr)

corr.list_b1_flare7 <- list()

# loop through chromosomes
for(i in 22:1){
  # set up file names
  load(paste0("R:/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE7_gds_steam/FLARE7_b1_chr_", i, ".Rdata"))
  corr.list_b1_flare7[[i]] <- snps.dt_corr
  # print update
  cat('done with chromosome', i, '\n')
  
}
saveRDS(corr.list_b1_flare7,  file = "corr.list_b1_FLARE7.rds" )

binsize <- 0.5
# computing variables needed for threshold estimation
for(i in 22:1){
  print(i)
  map.df <- FLARE7_map |> filter(chr == i)
  map.df$snp.id <- rownames(map.df)
  corr.list_b1_flare7[[i]]$cM = get_dist(corr.list_b1_flare7[[i]]$snp1, 
                                         corr.list_b1_flare7[[i]]$snp2, map.df)
  corr.list_b1_flare7[[i]]$theta = L_to_theta(corr.list_b1_flare7[[i]]$cM)
  corr.list_b1_flare7[[i]]$bin = round(corr.list_b1_flare7[[i]]$cM / binsize, 0) * binsize
}

corr_K4_FLARE7 <- combine_corr_chr(corr.list_b1_flare7)
saveRDS(corr_K4_FLARE7, file = "corr_K4_b1_FLARE7.rds")
get_g(corr_K4_FLARE7) #10.33934 -- number of generation estiamted since admixture

colnames(flare7_global_b1)[1:4] <- c("afr", "nam", "eur", "other4")
flare7_global_b1 <- flare7_global_b1 |> dplyr::select(afr, eur, nam, other4)
set.seed(1) # set seed for replication 
thresh_FLARE7_b1 <- get_thresh_simstat(g = get_g(corr_K4_FLARE7),
                                       map = FLARE7_map, 
                                       props = flare7_global_b1 , nreps = 10000)
