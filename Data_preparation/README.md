## Data Prepraration Code ##

**subset_FLARE7_gds_by_batch.R**:
Subsetting the FLARE7 gds file by metabolomic batches (batch 1 and 2 respectively), used in admixture mapping.

**subset_FLARE3_gds_by_batch.R**:
Subsetting the FLARE3 gds file by metabolomic batches (batch 1 and 2 respectively), used in admixture mapping.

**STEAM_create_genetic_maps_RFMix+FLARE3+7.R**:
Creating genetic maps by matching the snp position information in RFMix (using midpoint of a local ancestry interval), FLARE3, and FLARE7 to the GRCh38 map by chromosome. 

**SNPs_QC_filtergds_FLARE3+7.R**:
Code for filtering the gds local ancestry files from FLARE3 and FLARE7 using the QC metrics
(R^2 >= 0.8) and MAF (>= 0.005).

**SNPs_QC_by_imputation_quality_FLARE3+7.R**:
Quality control (QC) by imputation quality (R^2 >= 0.8) and MAF (>= 0.005) for FLARE3 and FLARE7 gds files

**SNPs_QC_by_imputation_quality_r2_0.95_FLARE3+7.R**:
Harsh quality control (QC) by imputation quality (R^2 >= 0.95) and MAF (>= 0.005) for FLARE3 and FLARE7 gds files, which is used for plotting FLARE vs. RFMix correlations in supplementary figure.

**RFMix_global_ancestry.R**
Code for computing global ancestry using the gds file from RFMix inference

**RFMix_ancestry_interal_liftover_prep.R**
Code for formatting RFMix local ancestry intervals into inputs to be lifted over in UCSC genome browser, from hg19 to hg38. 

**FLARE7_global_ancestry.R**
Code for computing global ancestry using the gds file from the FLARE7 inference.

**FLARE7_create_local_ancestry_gds.R**
Code for creating the FLARE7 local ancestry gds file from raw txt file. 

**FLARE3_low_corr_region_admixmap_cov_prep.R**
Code for preparing covariate file for a genomic region at chr16 with low FLARE3-RFMix local ancestry correlations, which was used for admixture mapping with propyl 4-hydroxybenzoate sulfate. 

**FLARE3_global_ancestry.R**
Code for computing global ancestry using the gds file from the FLARE3 inference.

**FLARE3_create_local_ancestry_gds.R**
Code for creating the FLARE3 local ancestry gds file from raw txt file. 

**FLARE+RFMix_create_SNPgds.R**
Code for converting gds to snpgds files for RFMix, FLARE3, and FLARE7, which were used in estimating multiple testing burden using STEAM. 

