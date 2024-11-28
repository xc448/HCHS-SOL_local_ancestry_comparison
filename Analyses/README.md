## Analyses conducted

**SNPs_QC_by_imputation_quality.R**:
Quality control (QC) by imputation quality and MAF for FLARE3 and FLARE7 gds files

**Correlation+ancestry_mean_var_proportion_FLARE3_RFMix.R**
**Combine_results_correlations_FLARE3_RFMix.R**
- Local ancestry correlation computation between FLARE3 and RFMix + computing mean & variance of the ancestry proportions
- Combine results from all chromosomes into one dataframe 

**Correlation+ancestry_mean_var_proportion_FLARE7_RFMix.R**
**Combine_results_correlations_FLARE7_RFMix.R**
- Local ancestry correlation computation between FLARE7 and RFMix + computing mean & variance of the ancestry proportions
- Combine results from all chromosomes into one dataframe

**Admixture_mapping_metabolites_RFMix+FLARE7_batch1.R** 
**Admixture_mapping_metabolites_RFMix+FLARE7_batch1.R** 
Admixture mapping for all four metabolites using either RFMix or FLARE7 local ancestry inferences
in batch 1 and batch 2: 3-aminoaceyltate, N-acetylarginine, PE 16:0/20:4, and PC 16:0/20:4, for African and Amerindian ancestries respectively

**Admixture_mapping_metabolites_FLARE3_batch1.R**
**Admixture_mapping_metabolites_FLARE3_batch2.R**
Admixture mapping for all four metabolites using the FLARE3 local ancestry inference 
in batch 1 and batch 2: 3-aminoaceyltate, N-acetylarginine, PE 16:0/20:4, and PC 16:0/20:4, for African and Amerindian ancestries respectively

**Admixture_mapping_low_corr_all_batches+inferencess.R**
Admixture mapping for metabolite propyl 4-hydroxybenzoate sulfate using all of the local ancestry inferences (RFMix, FLARE3, and FLARE7)
in both batch 1 and batch 2, for Africen ancestry at chromosome 16 only 

**High_corr_sampled_mapping_blacklist**
Randomly sampled high correlation SNP-LAI pairs between FLARE3 and RFMix, 
and mapped these pairs into either ENCODE blacklist or gene clusters annotated by UCSC genome brower.
Serves as a comparison to test the enrichment of the blacklisted regions in low-correlation pairs 


