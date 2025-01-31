# extract object from .RData file
get_obj <- function(Rdata){
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

# convert genetic distance into recombination fraction
L_to_theta <- function(cM){
  L <- cM/100
  theta <- 0.5*(1-exp(-2*L))
  return(theta)
}

# function to get distance between pairs of SNPs
get_dist <- function(s1, s2, gen.map){
  pos1 <- gen.map$cM[match(s1, gen.map$snp.id)]
  pos2 <- gen.map$cM[match(s2, gen.map$snp.id)]
  return(abs(pos2-pos1))
}

#' Calculate Correlation in Local Ancestry (K = 3)
#'
#' Calculate the correlation of local ancestry vectors
#' for a single chromosome. Code is currently only applicable to
#' admixed populations with three ancestral populations.
#'
#' @param chrom chromosome number that you are analyzing
#' @param binsize size (in cM) of distance bins for calculating correlation; default = 0.5 cM
#' @param map map file; data frame with, at minimum, columns 'chr' and 'cM'
#' @param pop1.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 1 ancestry, 0 = pop 2 or 3
#' @param pop2.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 2 ancestry, 0 = pop 1 or 3
#' @param pop3.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 3 ancestry, 0 = pop 1 or 2
#' @param verbose do you want to print updates to screen; default = TRUE
#'
#'
#' @return A data table with the observed correlation in local ancestry vectors for a subset of loci on this chromosome.
#'
#' @import data.table gdsfmt SeqArray SNPRelate
#'
#' @importFrom utils combn
#'
#' @importFrom stats cor
#'
#' @seealso \code{\link[STEAM]{get_g}} and \code{\link[STEAM]{combine_corr_chr}}
#'
#' @export
get_corr_chr_customized <- function(chrom, binsize = 0.5, map="", pop1.gds = "", 
                                    pop2.gds  = "", pop3.gds = "", snps.dt, verbose = TRUE){
  ## restrict map to chrom of interest
  snps.dt <- readRDS(snps.dt)
  
  # ## open GDS files
  print("opengds")
  pop1 <- SNPRelate::snpgdsOpen(pop1.gds)
  pop2 <- SNPRelate::snpgdsOpen(pop2.gds)
  pop3 <- SNPRelate::snpgdsOpen(pop3.gds)
  
  # set up correlation columns to store correlations
  snps.dt$corr_11 <- NA
  snps.dt$corr_12 <- NA
  snps.dt$corr_13 <- NA
  snps.dt$corr_14 <- NA
  snps.dt$corr_21 <- NA
  snps.dt$corr_22 <- NA
  snps.dt$corr_23 <- NA
  snps.dt$corr_24 <- NA
  snps.dt$corr_31 <- NA
  snps.dt$corr_32 <- NA
  snps.dt$corr_33 <- NA
  snps.dt$corr_34 <- NA
  snps.dt$corr_41 <- NA
  snps.dt$corr_41 <- NA
  snps.dt$corr_43 <- NA
  snps.dt$corr_44 <- NA
  
  # loop through pairs to get correlation
  for(i in 1:nrow(snps.dt)){
    # get snp names
    snp1.i <- snps.dt$snp1[i]
    snp2.i <- snps.dt$snp2[i]
    # get correlations
    print(paste0("getting corr, row",i))
    corr.i.afr <- get_corr_customized(snp1.i, snp2.i, anc_cor = "afr",
                                      afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_11[i] <- corr.i.afr[1]
    snps.dt$corr_12[i] <- corr.i.afr[2]
    snps.dt$corr_13[i] <- corr.i.afr[3]
    snps.dt$corr_14[i] <- corr.i.afr[4]
    rm(corr.i.afr)
    gc()
    print("afr done")
    corr.i.eur <- get_corr_customized(snp1.i, snp2.i, anc_cor = "eur",
                                      afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_21[i] <- corr.i.eur[1]
    snps.dt$corr_22[i] <- corr.i.eur[2]
    snps.dt$corr_23[i] <- corr.i.eur[3]
    snps.dt$corr_24[i] <- corr.i.eur[4]
    rm(corr.i.eur)
    gc()
    print("eur done")
    corr.i.amer <- get_corr_customized(snp1.i, snp2.i, anc_cor = "amer",
                                       afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_31[i] <- corr.i.amer[1]
    snps.dt$corr_32[i] <- corr.i.amer[2]
    snps.dt$corr_33[i] <- corr.i.amer[3]
    snps.dt$corr_34[i] <- corr.i.amer[4]
    rm(corr.i.amer)
    gc()
    print("amer done")
    corr.i.other <- get_corr_customized(snp1.i, snp2.i, anc_cor = "other",
                                        afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_41[i] <- corr.i.other[1]
    snps.dt$corr_42[i] <- corr.i.other[2]
    snps.dt$corr_43[i] <- corr.i.other[3]
    snps.dt$corr_44[i] <- corr.i.other[4]
    rm(corr.i.other)
    gc()
    print("other done")
    
    # progress report
    #if(i %% 1000 == 0) cat('Done with rep',i,'of',nrow(snps.dt),'\n')
  }
  # 
  # close gds files
  SNPRelate::snpgdsClose(pop1); SNPRelate::snpgdsClose(pop2); SNPRelate::snpgdsClose(pop3)
  
  # return data table with correlation
  return(snps.dt)
  
}


#' @param anc_cor the focal ancestry to compute correlations from
# get correlation for one pair of SNPs
get_corr_customized <- function(s1, s2, anc_cor, afr, eur, nam){
  # load ancestries at SNPs s1 and s2
  afr.anc <- snpgdsGetGeno(afr, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  eur.anc <- snpgdsGetGeno(eur, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  nam.anc <- snpgdsGetGeno(nam, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  other.anc <- nam.anc
  other.anc$genotype <- 2- (nam.anc$genotype + afr.anc$genotype + eur.anc$genotype)
  
  # get indices for each SNP (make sure they're in order we expected)
  e1 <- which(eur.anc$snp.id==s1); e2 <- which(eur.anc$snp.id == s2)
  a1 <- which(afr.anc$snp.id==s1); a2 <- which(afr.anc$snp.id == s2)
  n1 <- which(nam.anc$snp.id==s1); n2 <- which(nam.anc$snp.id == s2)
  o1 <- which(other.anc$snp.id==s1); o2 <- which(other.anc$snp.id==s2)
  # get correlation between all possible pairs of ancestries
  if(anc_cor == "afr"){
    c11 <- cor(afr.anc$g[,a1], afr.anc$g[,a2])
    c12 <- cor(afr.anc$g[,a1], eur.anc$g[,e2])
    c13 <- cor(afr.anc$g[,a1], nam.anc$g[,n2])
    c14 <- cor(afr.anc$g[,a1], other.anc$g[,o2])
  } else if(anc_cor == "eur"){
    c11 <- cor(eur.anc$g[,e1], afr.anc$g[,a2])
    c12 <- cor(eur.anc$g[,e1], eur.anc$g[,e2])
    c13 <- cor(eur.anc$g[,e1], nam.anc$g[,n2])
    c14 <- cor(eur.anc$g[,e1], other.anc$g[,o2])
  } else if(anc_cor == "amer"){
    c11 <- cor(nam.anc$g[,n1], afr.anc$g[,a2])
    c12 <- cor(nam.anc$g[,n1], eur.anc$g[,e2])
    c13 <- cor(nam.anc$g[,n1], nam.anc$g[,n2])
    c14 <- cor(nam.anc$g[,n1], other.anc$g[,o2])
  } else{
    c11 <- cor(other.anc$g[,o1], afr.anc$g[,a2])
    c12 <- cor(other.anc$g[,o1], eur.anc$g[,e2])
    c13 <- cor(other.anc$g[,o1], nam.anc$g[,n2])
    c14 <- cor(other.anc$g[,o1], other.anc$g[,o2])  
  }
  rm(afr.anc, eur.anc, nam.anc, other.anc)
  gc()
  return(c(c11,c12,c13,c14))}


#' Calculate Correlation in Local Ancestry (K = 4) 
#' -- the 4th ancestry is a combination of all other 4 ancestries in FLARE7
#' including Middle Eastern, East Asian, Central/South Asian, and Oceanian
#'
#' Create single data frame with average correlation of local ancestry
#' for pairs of SNPs varying distances apart. Combines data for each chromosome
#' generated by get_corr_chr(). Code is currently applicable
#' only to admixed populations with three ancestral populations.
#'
#' @param cor.list list of data tables, one per chromosome, generated by get_corr_chr().
#'
#' @return A data table with the average correlation in local ancestry, to be used for estimating the number of generations since admixture using \code{\link[STEAM]{get_g}}.
#'
#' @import data.table
#'
#' @importFrom stats reshape
#'
#' @seealso \code{\link[STEAM]{get_g}} and \code{\link[STEAM]{get_corr_chr}}
#'
#' @export
combine_corr_chr <- function(cor.list){
  # combine into single data table
  snps.df <- rbindlist(cor.list)
  
  # get average within each bin
  avgs <- snps.df[,.(mean(corr_11),mean(corr_12),mean(corr_13),mean(corr_14),
                     mean(corr_21),mean(corr_22),mean(corr_23),mean(corr_24),
                     mean(corr_31),mean(corr_32),mean(corr_33),mean(corr_34),
                     mean(corr_41),mean(corr_42),mean(corr_43),mean(corr_44)), by = .(bin)]
  avgs[,('theta') := L_to_theta(bin)] # add recomb fraction back on
  
  # add col names
  names(avgs) <- c('cM',
                   'corr_11','corr_12','corr_13','corr_14',
                   'corr_21','corr_22','corr_23','corr_24',
                   'corr_31','corr_32','corr_33','corr_34',
                   'corr_41','corr_42','corr_43','corr_44',
                   'theta')
  
  # convert from wide to long
  avgs.long <- reshape(avgs, direction = 'long', varying = list(2:10),
                       v.names = 'corr', idvar = 'bin', timevar = 'anc',
                       times = c('1_1','1_2','1_3','1_4','2_1','2_2','2_3',
                                 '2_4','3_1','3_2','3_3','3_4','4_1','4_2','4_3','4_4'))
  
  avgs.long$anc2 <- avgs.long$anc
  avgs.long$anc2[avgs.long$anc2=='2_1'] <- '1_2'
  avgs.long$anc2[avgs.long$anc2=='3_1'] <- '1_3'
  avgs.long$anc2[avgs.long$anc2=='4_1'] <- '1_4'
  avgs.long$anc2[avgs.long$anc2=='3_2'] <- '2_3'
  avgs.long$anc2[avgs.long$anc2=='4_2'] <- '2_4'
  avgs.long$anc2[avgs.long$anc2=='4_3'] <- '3_4'
  
  # save final results as data frame
  lacorr.df <- as.data.frame(avgs.long)
  lacorr.df <- lacorr.df[,c('theta','corr','anc2')]
  names(lacorr.df) <- c('theta','corr','anc')
  return(lacorr.df)
}
