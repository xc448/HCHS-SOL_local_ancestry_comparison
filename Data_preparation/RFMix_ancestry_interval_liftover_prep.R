path <- "R:\\Sofer Lab\\HCHS_SOL\\Ancestry_files\\UW_GAC_DMDA_20180516_local_ancestries"
#path <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/UW_GAC_DMDA_20180516_local_ancestries"

setwd(path)

annot <- read.csv("lai_HGDP_1000G_comb_snpAnnot_unique.csv")
intervalid <-annot$snpID
pos <- as.data.frame(paste("chr", annot$chromosome, ":",  annot$pos.start+1,"-" , annot$pos.end, sep = ""))
bed <- as.data.frame(cbind(paste0("chr", annot$chromosome), annot$pos.start, annot$pos.end, intervalid))


write.table(bed, 
            quote = FALSE,
            file = "/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/liftover37.bed",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

# Notes from the annot file:
# column "position" is simply the midpoint in the interval
# column snpName is the ID of the SNP right at the beginning of the interval 