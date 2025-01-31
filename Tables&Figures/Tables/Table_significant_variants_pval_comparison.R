# Comparison of p-values from the most significant variants (p*) detected in the
# most significant genomic regions from RFMix admixture mapping results. 
# The p* values are compared across the three inferences (RFMix, FLARE3, and FLARE7) 
# as well as between the discovery batch 1 (b1) and replication batch 2 (b2) 

# load libraries
library(gdsfmt)
library(tidyverse)
library(kableExtra)
library(knitr)
library(htmltools)

#load the liftover file from hg37 to hg38 in RFMix intervals
#intervals_liftover_nodup <- read.table("R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\intervals_37-38.bed")
intervals_liftover_nodup <- read.table("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/intervals_37-38.bed")
colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 
intervals_liftover_nodup <- intervals_liftover_nodup |> filter(chr != "chrX") # 14753 mapped blocks

### find the most significant region in RFMix results ###
# arguments: metabolite, ancestry (African or Amerindian), and batch number to
# locate the RFMix admixture mapping results from the local directory 
sigRegion_RFMix <- function(metab, ancestry, batch){
  res_admixmap <- readRDS(paste0("../Code/admix_map_all/admixmap_", metab, "_RFMix/", "RFMix_res_", ancestry,
                                 "_", metab, "_b", batch, "_fixed.rds"))
  ind <- rownames(res_admixmap[which.min(res_admixmap$Score.pval),])
  sig_region <- cbind(intervals_liftover_nodup[intervals_liftover_nodup$snpid == ind,], res_admixmap[which.min(res_admixmap$Score.pval),])
  return(sig_region)
}

# flare7 <- "R:\\Sofer Lab\\HCHS_SOL\\Projects\\2023_local_ancestry_comparison_sol\\Data\\ancestry_counts_FLARE.gds"
flare7 <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE7/FLARE7_SNPs_filtered_combined.rds") 
chromid <- flare7$chr
rsID <- flare7$rsID
 
### find the most significant variant in FLARE7 results given a genomic region##
combine_find_variant_FLARE7 <- function(metab, ancestry, batch, rfmix_interval,
                                        chr_region, sigsnpb1 = "",
                                        from = 1, to = 22){
  res_flare7 <- c()
  for(i in from:to){
    load(paste0("./admix_map_all/admixmap_", metab, "_all_FLARE7/admix", 
                ancestry, "_b", batch, "_filtered_fixed_chr" , i , ".Rdata"))
    res_flare7 <-  rbind(res_flare7, as.data.frame(cbind(res$snppos, res$Score.pval)))
  }
  res_flare7 <- cbind(chromid, rsID, res_flare7)
  colnames(res_flare7) <- c("chr", "rsID", "snppos", "Score.pval")
  saveRDS(res_flare7, paste0("./admix_map_all/admixmap_", metab, 
                             "_all_FLARE7/admix", ancestry, "_filtered_b", batch, 
                             "_combined.rds"))
  
  chr <- rfmix_interval$chr
  if(batch ==  1){
    Var_filtered <-  res_flare7 |> filter(chr == chr) |>
      filter(snppos >= chr_region[1] & snppos < chr_region[2])
    return(Var_filtered[which.min(Var_filtered$Score.pval),])
  }
  else{
    return(res_flare7[res_flare7$rsID == sigsnpb1,])
  }
  
}

### find the most significant variant in FLARE3 results given a genomic region##
combine_find_variant_FLARE3 <- function(metab, ancestry, batch, rfmix_interval, 
                                        chr_region, sigsnpb1 = "",
                                        from = 1, to = 22){
  res_flare3 <- c()
  totrows <- 0
  for(i in from:to){
    load(paste0("./admix_map_all/admixmap_", metab, "_all_FLARE3/admix", ancestry, 
                "_b", batch, "_filtered_fixed_chr" , i , ".Rdata"))
    pvals <-  as.data.frame(cbind(res$snppos, res$Score.pval))
    
    snpinfo <- readRDS(paste0("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr",
                                i, ".rds"))
    chromid_flare3 <- snpinfo$chr
    rsid_flare3 <- snpinfo$rsID
    res_flare3 <- rbind(res_flare3, cbind(chromid_flare3, rsid_flare3, pvals))
    totrows <- totrows + nrow(snpinfo)
  }
  colnames(res_flare3) <- c("chr", "rsID", "snppos", "Score.pval")
  saveRDS(res_flare3, paste0("./admix_map_all/admixmap_", metab,
                             "_all_FLARE3/admix", ancestry, "_filtered_b", batch, "_combined.rds"))
  
  if(batch == 1){
    chr <- rfmix_interval$chr
    Var_filtered <-  res_flare3 |> filter(chr == chr) |>
      filter(snppos >= chr_region[1] & snppos < chr_region[2])
    return(Var_filtered[which.min(Var_filtered$Score.pval),])
  }else{
    return(res_flare3[res_flare3$rsID == sigsnpb1,])
  }
  }

# Helper function to create rows for the table, which includes the genomic 
# regions where the most significant ancestral difference was observed in RFMix,
# the most significant variant rsID and p-value in that genomic region for FLARE3 and FLARE7 results
createRows <- function(metab, ancestry, chr_region, from = 1, to = 22){
  row <- c()
  for(i in 1:2){
    rfmix_interval <-sigRegion_RFMix(metab, ancestry, i)
    if(i == 1){
      sigFLARE7 <- combine_find_variant_FLARE7(metab, ancestry, i, 
                                               rfmix_interval, chr_region,
                                               from = from, to = to)
      sigFLARE3 <- combine_find_variant_FLARE3(metab, ancestry, i, 
                                               rfmix_interval, chr_region,
                                               from = from, to = to)
    } else {
      sigFLARE3 <- combine_find_variant_FLARE3(metab, ancestry, i, rfmix_interval, 
                                               chr_region, row[1,2],
                                               from = from, to = to)
      sigFLARE7 <- combine_find_variant_FLARE7(metab, ancestry, i, rfmix_interval, 
                                               chr_region, row[1,4],
                                               from = from, to = to)
      
    }
    row <- cbind(row, cbind(rfmix_interval$Score.pval, sigFLARE3$rsID, 
                              sigFLARE3$Score.pval,
                              sigFLARE7$rsID, sigFLARE7$Score.pval))
    print(row) 
    
  }
  return(row)
}

# The genomic regions were obtained from the UCSC genome brower
x1114_amer <- createRows("x1114","amer", c(33800001, 38400000)) # chr5(p13.2)
x1266_afr <- createRows("x1266","afr", c(71300001, 73300000)) # chr2(p13.2)
x8990_amer <- createRows("x8990","amer", c(52600001, 58800000)) # chr15(q21.3)
x8914_afr <- createRows("x8914","afr", c(60100001, 61900000)) # chr11(q12.2)
x8914_amer <- createRows("x8914","amer", c(60100001, 61900000)) # chr11(q12.2)


## Make the final table by adjusting notations and variables to include ########
final_table <- c()
metab <- c("N-acetylarginine", "3-aminoisobutyrate", rep("PC 16:0/20:4",2 ), "PE 16:0/20:4")
ancestry <- c("African",  "Amerindian",  "African", "Amerindian", "Amerindian")
region <- c("chr2(p13.2)", "chr5(p13.2)", rep("chr11(q12.2)", 2), "chr15(q21.3)") 
final_table  <- as.data.frame(cbind(metab, ancestry, region, 
                                    rbind( x1266_afr, x1114_amer,
                                           x8914_afr, x8914_amer,
                                           x8990_amer)))


final_table[, c(4,6,8,9,11,13)] <- apply(final_table[, c(4,6,8,9,11,13)], 
                                         2, function(x){x = as.numeric(x)
                                         formatC(x, format = "e", digits = 3) })


colnames(final_table) <- c("Metabolite", "Ancestry",  "Region(RFMix)",
                           "p* RFMix b1", "rsID FLARE3 b1", "p* FLARE3 b1",
                           "rsID FLARE7 b1", "p* FLARE7 b1", "p* RFMix b2", 
                           "rsID FLARE3 b2", "p* FLARE3 b2",
                           "rsID FLARE7 b2", "p* FLARE7 b2" )

final_table <- final_table |> dplyr::select(Metabolite, Ancestry, 
                                             `Region(RFMix)`,  `p* RFMix b1`,
                                             `p* RFMix b2`, `rsID FLARE3 b1`,
                                             `p* FLARE3 b1`, `p* FLARE3 b2`,
                                             `rsID FLARE7 b1`, `p* FLARE7 b1`,
                                             `p* FLARE7 b2`)

colnames(final_table)[c(6, 9)] <- c("rsID FLARE3", "rsID FLARE7")
  
# Visualize the table and save it as a html file
kable(final_table) %>%
  kable_classic( full_width = F,
                html_font = "Arial") %>%
  row_spec(0, bold = TRUE, align = "center")  %>%
  column_spec(1,bold = T, width =  "4cm")

