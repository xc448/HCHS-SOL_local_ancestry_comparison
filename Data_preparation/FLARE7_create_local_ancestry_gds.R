# Create local ancestry GDS file for FLARE7 dataset from the original .txt file

path <- "/Volumes/Sofer Lab/HCHS_SOL/2023_FLARE_7pop" # Mac
path <- "R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop" # Windows
setwd(path)
library(GWASTools)
library(gdsfmt)
library(data.table)

# reading files
id <- read.delim("SOL.id.txt",header = FALSE) # subject id
snp <- read.delim(header = FALSE, "SOL.snp.txt") # snp information, chr, position
chrnum <- snp$V1
snppos <- snp$V2
snpid <- c(1:5105005)
sample_nums_id <- c(1:11929)
a1 <- snp$V4 # allele 1
a2 <- snp$V5 # allele 2

# adding snp information and alleles information to the GDS file
add.gdsn(gfile, "snp.id", snpid, compress = "LZMA_RA:1M")
add.gdsn(gfile, "snp.chromosome", chrnum, compress = "LZMA_RA:1M" )
add.gdsn(gfile, "snp.position", snppos, compress = "LZMA_RA:1M" )
add.gdsn(gfile, "sample.id", sample_nums_id, compress = "LZMA_RA:1M")
add.gdsn(gfile, "a1", a1, compress = "LZMA_RA:1M" )
add.gdsn(gfile, "a2", a2, compress = "LZMA_RA:1M" )

counts <- function(ancestry_number, anc1, anc2) {
  # Output: a vector containing calculated ancestry counts for each 
  # individual (columns) for each SNP given the ancestry of interest.
  
  # Input: ancestry_number: number representing different ancestries (0-6)
  # anc1: copy #1 showing the ancestry information for each SNP (rows) on 
  # each individual(columns)
  # anc2:  copy #2 showing the ancestry information for each SNP (rows) on 
  # each individual(columns)
  tmp <- as.matrix((anc1 == ancestry_number) + (anc2 == ancestry_number)) # convert it into a matrix 
  # create a df with ancestry counts 
  return(tmp) 
}

# Use gdsfmt package:
# using LZMA_RA:1M as compressing method to calculate the ancestry counts for the first 
# 5000000 rows in the anc file using a for loop and append them to the gds file
gfile <- createfn.gds("global_ancestry_final.gds")
for(i in 1:10){
  print(i)
  anc1 <- read.delim("SOL.anc1.txt", header = FALSE, nrows = 500000, skip = 500000*(i-1)) 
  anc1 <- anc1[,-c(11929)]
  anc2 <- read.delim("SOL.anc2.txt", header = FALSE, nrows = 500000, skip =  500000*(i-1))
  anc2 <- anc2[,-c(11929)]
  geno <- read.table("R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop\\SOL.geno.txt",  header = FALSE, nrows = 500000, skip = 500000*(i-1))
  print("reading done")
  
  # compute local ancestry counts (0, 1, or 2 copies at each SNP position)
  afr_counts <- counts(0, anc1, anc2)
  easia_counts <- counts(1, anc1, anc2)
  eur_counts <- counts(2, anc1, anc2)
  csasia_counts <- counts(3, anc1, anc2)
  amer_counts <- counts(4, anc1, anc2)
  ocea_counts <- counts(5, anc1, anc2)
  me_counts <- counts(6, anc1, anc2)
  
  if(i == 1){
    # adding local ancestry count for each ancestry as nodes into the gds file 
    node_afr <- add.gdsn(gfile, "afr_counts", t(afr_counts), compress = "LZMA_RA:1M")
    node_easia <- add.gdsn(gfile, "easia_counts", t(easia_counts), compress = "LZMA_RA:1M")
    node_eur <- add.gdsn(gfile, "eur_counts", t(eur_counts), compress = "LZMA_RA:1M")
    node_csasia <- add.gdsn(gfile, "csasia_counts", t(csasia_counts), compress = "LZMA_RA:1M")
    node_amer <- add.gdsn(gfile, "amer_counts", t(amer_counts), compress = "LZMA_RA:1M")
    node_ocea <- add.gdsn(gfile, "ocea_counts", t(ocea_counts), compress = "LZMA_RA:1M")
    node_me <- add.gdsn(gfile, "me_counts", t(me_counts), compress = "LZMA_RA:1M")
    
    # genotype information 
    node_geno <- add.gdsn(Flare_gds, "genotype", t(geno), compress = "LZMA_RA:1M")
  }
  else{
    append.gdsn(node_afr, t(afr_counts))
    append.gdsn(node_easia, t(easia_counts))
    append.gdsn(node_eur, t(eur_counts))
    append.gdsn(node_csasia, t(csasia_counts))
    append.gdsn(node_amer, t(amer_counts))
    append.gdsn(node_ocea, t(ocea_counts))
    append.gdsn(node_me, t(me_counts))
    append.gdsn(node_geno, t(geno))
  }
}

# taking care of the last 105005 rows in the anc file
anc1 <- read.delim("SOL.anc1.txt", header = FALSE, nrows = 105005, skip = 5000000) 
anc2 <- read.delim("SOL.anc2.txt", header = FALSE, nrows = 105005, skip = 5000000)
afr_counts <- counts(0, anc1, anc2)
easia_counts <- counts(1, anc1, anc2)
eur_counts <- counts(2, anc1, anc2)
csasia_counts <- counts(3, anc1, anc2)
amer_counts <- counts(4, anc1, anc2)
ocea_counts <- counts(5, anc1, anc2)
me_counts <- counts(6, anc1, anc2)

# genotype information 
geno <-  read.table("R:\\Sofer Lab\\HCHS_SOL\\2023_FLARE_7pop\\SOL.geno.txt",  header = FALSE, nrows = 105005, skip = 5000000)
append.gdsn(node_geno, t(geno))

# accessing the stored file
node <- index.gdsn(gfile, node_afr)
readmode.gdsn(node) # have to switch to read mode first 
read.gdsn(node, start=c(1000, 1000), count=c(2, 3)) # read from specified starting position

# close the gds file
closefn.gds(gfile)
