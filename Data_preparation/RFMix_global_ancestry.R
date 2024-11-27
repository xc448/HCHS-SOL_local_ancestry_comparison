# path <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries"
path <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries"
setwd(path)
library(GWASTools)
library(gdsfmt)
library(GenomicRanges) 

# read annotation 
annot <- getobj("lai_HGDP_1000G_comb_snpAnnot_unique.RData")
# open the gds file
gds <- openfn.gds("lai_HGDP_1000G_comb_unique.gds")
# the dosage_afr, dosage_amer, and dosage_eur nodes store the local ancestry calls
# set up GdsGenotypeReader objects for the three different ancestries
# this call is slightly different than for normal genotype gds files,
# since the nodes have different names
genoDataList <- list()
for (ancestry in c("afr", "eur", "amer")){
  genotype_variable <- sprintf("dosage_%s", ancestry)
  genoDataList[[ancestry]] <- GdsGenotypeReader(gds, genotypeVar=genotype_variable)
}

# computing interval lengths, which is used in computing global proportions
intervals <- annot$pos.end - annot$pos.start
afr <- getGenotype(genoDataList[["afr"]])
eur <- getGenotype(genoDataList[["eur"]])
amer <- getGenotype(genoDataList[["amer"]])

global_proportion <- function(ancestry){
  tmp <- c()
  for(i in 1:ncol(afr)){
    tmp <- cbind(tmp, sum(ancestry[,i]*intervals, na.rm =  TRUE)/sum((ancestry[,i]*0+2) * intervals, na.rm =  TRUE))
  } # to handle NA values by dropping them: for example, scanID 84 has 685 NAs (4.419%), the formula above dropped all the intervals with NAs
  # in both numerator and denominator. 
  return(t(tmp))
}

# sum(is.na(afr[,84]))
global_afr <- global_proportion(afr)
global_eur <- global_proportion(eur)
global_amer <- global_proportion(amer)
ids <- as.data.frame(read.gdsn(index.gdsn(gds, "sample.id"))) 
colnames(ids) <- "ID"
sample_info <- read.csv("/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/subject_annotation_2017-09-05.csv")
global_proportion_old <- as.data.frame(cbind(sample_ID = ids, afr = global_afr, eur = global_eur, amer = global_amer))
global_proportion_solid <- right_join(sample_info[,1:2], global_proportion_old, by = c("GAC_scanID" = "ID"))|>
  group_by(SUBJECT_ID) |>
  summarize(mean_afr = mean(afr), mean_eur = mean(eur), mean_amer = mean(amer)) |>
  filter(!is.na(SUBJECT_ID))

write.csv(global_proportion_solid, "../RFMix_global_ancestry_from_LAI.csv", row.names = FALSE )
# total = 12689 individuals 

# making sure everything adds up equal to 1:
unique(global_afr + global_eur + global_amer)
