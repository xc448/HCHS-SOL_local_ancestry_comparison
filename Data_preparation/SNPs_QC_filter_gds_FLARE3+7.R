library(gdsfmt)
library(tidyverse)
library(GWASTools)

# load SNP information for FLARE7 dataset
flare7 <- read.table("R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop/SOL.snp.txt")
colnames(flare7) <- c("chr", "snppos", "rsID", "a1", "a2")

# load gds for FLARE7
flare7_gds <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/ancestry_counts_FLARE.gds"
snpids_flare7 <- readRDS("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_combined.rds")
filtered_id <- match(snpids_flare7$rsID,  flare7$rsID)

# Extracting snp position and chromosome information from FLARE7
all_snps  <- as.data.frame(read.gdsn(index.gdsn(flare7_gds,  "snp.position")))
all_chr <- as.data.frame(read.gdsn(index.gdsn(flare7_gds,  "snp.chromosome")))
all_snps <- cbind(all_chr, all_snps)
colnames(all_snps)  <- c("chr", "snppos")

# filter SNPs that pass the QC criteria (R2 > 0.8 if imputed, MAF >= 0.005 )
snp_qc <- all_snps |> mutate(original_id = rownames(all_snps)) |>
  inner_join(snpids_flare7, by = c("chr" = "chr", "snppos" = "snppos"))

subsetted_path <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/ancestry_counts_FLARE7_QC.gds"

# Subsetting the FLARE7 based on the QC file
path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7\\ancestry_counts_FLARE.gds")
gds <- openfn.gds(path)
chrinfo <- as.data.frame(read.gdsn(index.gdsn(gds,  "snp.chromosome")))
snppos <- read.gdsn(index.gdsn(gds,  "snp.position"))
snpid  <- as.data.frame(read.gdsn(index.gdsn(gds,  "snp.id")))
rownames(snpid) <- snpid[,1]
closefn.gds(gds)

for(i in 22:1){
  print(paste0("working on chromosome", i))
  startind <- min(which(chrinfo[,1] == i)) # starting to read ancestry counts from this index
  endind <- max(which(chrinfo[,1] == i))
  filtered_snps <- readRDS(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr",i, ".rds"))
  all_snps  <- as.data.frame(cbind(snpid[seq(startind, endind),],
                                   snppos[seq(startind, endind)]))
  colnames(all_snps) <- c("snp.id", "snppos")
  filtered_snps <- filtered_snps |> inner_join(all_snps, by = c("snppos"))
  subsetted_path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr", i, ".gds")
  gdsSubset(path, subsetted_path, snp.include=filtered_snps$snp.id, compress = "LZMA_RA")
  subset_gds <- openfn.gds(subsetted_path,  readonly  = FALSE)
  add.gdsn(subset_gds, name = "afr_counts", 
           val = read.gdsn(index.gdsn(subset_gds, "afr_counts")),
           compress = "LZMA_RA", replace = TRUE)
  add.gdsn(subset_gds, name = "amer_counts", 
           val = read.gdsn(index.gdsn(subset_gds, "amer_counts")),
           compress = "LZMA_RA", replace = TRUE)
  add.gdsn(subset_gds, name = "eur_counts", 
           val = read.gdsn(index.gdsn(subset_gds, "eur_counts")),
           compress = "LZMA_RA", replace = TRUE)
  closefn.gds(subset_gds)
  
}

# Subsetting the FLARE3 based on the QC file
for(i in 1:22){ # iterate through each chromosome
  print(paste0("working on chr ", i))
  path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3\\FLARE3_ancestry_counts_chr", i, ".gds")
  gds <- openfn.gds(path)
  filtered_snps <- readRDS(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr",i, ".rds"))
  all_snps  <- as.data.frame(cbind(read.gdsn(index.gdsn(gds, "snp.id")),
                                   read.gdsn(index.gdsn(gds, "snp.position"))))
  closefn.gds(gds)
  colnames(all_snps) <- c("snp.id", "snppos")
  filtered_snps <- filtered_snps |> inner_join(all_snps, by = c("snppos"))
  subsetted_path <- paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr", i, ".gds")
  gdsSubset(path, subsetted_path, snp.include=filtered_snps$snp.id, compress = "LZMA_RA")
  subset_gds <- openfn.gds(subsetted_path,  readonly  = FALSE)
  compression.gdsn(index.gdsn(subset_gds, "afr_counts"), compress = "LZMA_RA")
  compression.gdsn(index.gdsn(subset_gds, "amer_counts"), compress = "LZMA_RA")
  compression.gdsn(index.gdsn(subset_gds, "eur_counts"), compress = "LZMA_RA")
  closefn.gds(subset_gds)
}


