library(tidyverse)

# For the four main metabolites prepare the tab-limited files separately 

createLocusZoomdf <- function(inference, batch, metabolite, chr, ancestry){
  snpinfo <- readRDS(paste0("./admix_map_all/admixmap_x", metabolite, 
                            "_all_FLARE", inference, "/admix", ancestry,
                            "_filtered_", batch, "_combined.rds"))
  load(paste0("./admix_map_all/admixmap_x", metabolite, 
                            "_all_FLARE", inference, "/admix", ancestry,"_", batch, 
              "_filtered_fixed_chr", chr, ".Rdata"))
  df <- cbind(rep(chr, nrow(res)), res$snppos, -log(res$Score.pval, 10),
              res$Est)
  colnames(df) <- c("Chromosome", "Position", "-log10pval", "Effect_size")
  write.table(df, file = paste0("./admix_map_all/admixmap_x", metabolite, 
                                       "_all_FLARE", inference, "/", ancestry,"_", batch, 
                                       "chr", chr, "_locuszoom.txt"), sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# 3-Aminoisobutyrate at chromosome 5 
# Load file from FLARE3
# Amerindian ancestry (the driving ancestry)
createLocusZoomdf(3, "b1", 1114, 5, "amer")
# batch2
createLocusZoomdf(3, "b2", 1114, 5, "amer")

# Load file from FLARE7
# Amerindian ancestry (the driving ancestry)
createLocusZoomdf(7, "b1", 1114, 5, "amer")
# batch2
createLocusZoomdf(7, "b2", 1114, 5, "amer")

# N-acetylarginine at chromosome 2
# African ancestry (the driving ancestry)
# Load file from FLARE3
# batch1
createLocusZoomdf(3, "b1", 1266, 2, "afr")
# batch2
createLocusZoomdf(3, "b2", 1266, 2, "afr")

# Load file from FLARE7
# batch1
createLocusZoomdf(7, "b1", 1266, 2, "afr")
# batch2
createLocusZoomdf(7, "b2", 1266, 2, "afr")

# PE 16:0/20:4 at chromosome 15
# Amerindian ancestry (the driving ancestry)
# Load file from FLARE3
# batch1
createLocusZoomdf(3, "b1", 8990, 15, "amer")
# batch2
createLocusZoomdf(3, "b2", 8990, 15, "amer")

# Load file from FLARE7
# batch1
createLocusZoomdf(7, "b1", 8990, 15, "amer")
# batch2
createLocusZoomdf(7, "b2", 8990, 15, "amer")

# PC 16:0/20:4 at chromosome 11
# African ancestry (significant)
# batch1
createLocusZoomdf(3, "b1", 8914, 11, "afr")
# batch2
createLocusZoomdf(3, "b2", 8914, 11, "afr")

# Amerindian ancestry (the driving ancestry)
createLocusZoomdf(3, "b1", 8914, 11, "amer")
# batch2
createLocusZoomdf(3, "b2", 8914, 11, "amer")

# Load file from FLARE7
# African ancestry (significant)
# batch1
createLocusZoomdf(7, "b1", 8914, 11, "afr")
# batch2
createLocusZoomdf(7, "b2", 8914, 11, "afr")

# Amerindian ancestry (the driving ancestry)
createLocusZoomdf(7, "b1", 8914, 11, "amer")
# batch2
createLocusZoomdf(7, "b2", 8914, 11, "amer")