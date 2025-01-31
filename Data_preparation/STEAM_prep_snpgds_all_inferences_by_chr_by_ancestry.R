library(gdsfmt)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(GWASTools)


######## Create snpGDS for batch1 and batch 2 #########

# snpGDS --rfmix

for(j in c("b1","b2")){
  g <- openfn.gds(paste0("./uwgds_", j, "_subset_37"))
  chrinfo <- as.data.frame(read.gdsn(index.gdsn(g, "snp.chromosome")))
  colnames(chrinfo) <- c("chr")
  for(i in 22:1){
    print(paste0("working on chromosome", i))
    startind <- min(which(chrinfo$chr == i)) # starting to read ancestry counts from this index
    endind <- max(which(chrinfo$chr == i))
    for(ancestry in c("afr",  "amer", "eur")){
      snpgdsCreateGeno(paste0("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/RFMix_gds_by_chr/RFMix_", 
                              ancestry,"_",j, "_chr", i, ".gds"), 
                       genmat = read.gdsn(index.gdsn(g, paste0("dosage_", ancestry)), 
                                          start = c(1,startind), 
                                          count = c(-1, endind-startind+1)), 
                       sample.id = read.gdsn(index.gdsn(g, "sample.id")),
                       snp.id = read.gdsn(index.gdsn(g, "snp.id"))[startind:endind],
                       snp.chromosome=rep(i, endind-startind+1),
                       snp.position = read.gdsn(index.gdsn(g, "snp.position"))[startind:endind],
                       compress.annotation="LZMA_ra",
                       compress.geno =  "LZMA_ra",
                       snpfirstdim=FALSE)
    }
  }
  
}



# snpGDS --flare7
for(j in c("b1", "b2")){
  g <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_admixmap/flare_", j, "_all.gds"))
  chrinfo <- as.data.frame(read.gdsn(index.gdsn(g, "snp.chromosome")))
  colnames(chrinfo) <- c("chr")
  
  for(i in 22:1){
    print(paste0("working on chromosome", i))
    startind <- min(which(chrinfo$chr == i)) # starting to read ancestry counts from this index
    endind <- max(which(chrinfo$chr == i))
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7/FLARE7_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    chr <- as.data.frame(read.gdsn(index.gdsn(g, "snp.chromosome")))[filtered_id,]
    for(ancestry in c("afr",  "amer", "eur")){
      genmat <-  read.gdsn(index.gdsn(g, paste0(ancestry, "_counts")), 
                           start = c(1,startind), 
                           count = c(-1, endind-startind+1))
      snppos <- read.gdsn(index.gdsn(g, "snp.position"))[filtered_id]
      colnames(genmat) <- read.gdsn(index.gdsn(g, "snp.id"))[startind:endind]
      genmat <- genmat[, colnames(genmat) %in% filtered_id]
      snpgdsCreateGeno(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE7_gds_steam/FLARE7_", 
                              ancestry,"_",j, "_chr", i, ".gds"), 
                       genmat = genmat, 
                       sample.id = read.gdsn(index.gdsn(g, "sample.id")),
                       snp.id = colnames(genmat),
                       snp.chromosome=chr,
                       snp.position = snppos,
                       compress.annotation="LZMA_ra",
                       compress.geno =  "LZMA_ra",
                       snpfirstdim=FALSE)
    }
    closefn.gds(filtered_gds)
  }
  closefn.gds(g)
}


#snpgds flare3
for(j in c("b1", "b2")){
  for(i in c(22:1)){
    g <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_admixmap\\flare3_", j, "_chr", i, ".gds"))
    chrinfo <- as.data.frame(read.gdsn(index.gdsn(g, "snp.chromosome")))
    colnames(chrinfo) <- c("chr")
    print(paste0("working on chromosome", i))
    filtered_gds <- openfn.gds(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3/FLARE3_SNPs_filtered_chr", 
                                      i, ".gds" ))
    filtered_id  <-  read.gdsn(index.gdsn(filtered_gds, "snp.id"))
    for(ancestry in c("afr",  "amer", "eur")){
      genmat <-  read.gdsn(index.gdsn(g, paste0(ancestry, "_counts")))
      snppos <- as.data.frame(read.gdsn(index.gdsn(g, "snp.position")))
      rownames(snppos) <-  read.gdsn(index.gdsn(g, "snp.id"))
      snppos <- snppos[rownames(snppos)  %in% filtered_id,]
      colnames(genmat) <- read.gdsn(index.gdsn(g, "snp.id"))
      genmat <- genmat[, colnames(genmat) %in% filtered_id]
      snpgdsCreateGeno(paste0("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\FLARE3_gds_steam\\FLARE3_", 
                              ancestry,"_",j, "_chr", i, ".gds"), 
                       genmat = genmat, 
                       sample.id = read.gdsn(index.gdsn(g, "sample.id")),
                       snp.id = colnames(genmat),
                       snp.chromosome=rep(i, length(snppos)),
                       snp.position = snppos,
                       compress.annotation="LZMA_ra",
                       compress.geno =  "LZMA_ra",
                       snpfirstdim=FALSE)
    }
    closefn.gds(filtered_gds)
    closefn.gds(g)
  }
}
