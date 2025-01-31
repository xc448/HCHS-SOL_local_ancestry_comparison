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
  snps.dt$corr_21 <- NA
  snps.dt$corr_22 <- NA
  snps.dt$corr_23 <- NA
  snps.dt$corr_31 <- NA
  snps.dt$corr_32 <- NA
  snps.dt$corr_33 <- NA
  
  # loop through pairs to get correlation
  for(i in 1:nrow(snps.dt)){
    # get snp names
    snp1.i <- snps.dt$snp1[i]
    snp2.i <- snps.dt$snp2[i]
    # get correlations
    print(paste0("getting corr, row",i))
    corr.i.afr <- get_corr_customized(snp1.i, snp2.i, anc_cor = "afr",
                                      afr = pop1, eur = pop2, nam = pop3)
    corr.i.eur <- get_corr_customized(snp1.i, snp2.i, anc_cor = "eur",
                                      afr = pop1, eur = pop2, nam = pop3)
    corr.i.amer <- get_corr_customized(snp1.i, snp2.i, anc_cor = "amer",
                                       afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_11[i] <- corr.i.afr[1]
    snps.dt$corr_12[i] <- corr.i.afr[2]
    snps.dt$corr_13[i] <- corr.i.afr[3]
    snps.dt$corr_21[i] <- corr.i.eur[1]
    snps.dt$corr_22[i] <- corr.i.eur[2]
    snps.dt$corr_23[i] <- corr.i.eur[3]
    snps.dt$corr_31[i] <- corr.i.amer[1]
    snps.dt$corr_32[i] <- corr.i.amer[2]
    snps.dt$corr_33[i] <- corr.i.amer[3]
    # progress report
    if(i %% 1000 == 0) cat('Done with rep',i,'of',nrow(snps.dt),'\n')
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
  # get indices for each SNP (make sure they're in order we expected)
  e1 <- which(eur.anc$snp.id==s1); e2 <- which(eur.anc$snp.id == s2)
  a1 <- which(afr.anc$snp.id==s1); a2 <- which(afr.anc$snp.id == s2)
  n1 <- which(nam.anc$snp.id==s1); n2 <- which(nam.anc$snp.id == s2)
  # get correlation between all possible pairs of ancestries
  if(anc_cor == "afr"){
    c11 <- cor(afr.anc$g[,a1], afr.anc$g[,a2])
    c12 <- cor(afr.anc$g[,a1], eur.anc$g[,e2])
    c13 <- cor(afr.anc$g[,a1], nam.anc$g[,n2])
  } else if(anc_cor == "eur"){
    c11 <- cor(eur.anc$g[,e1], afr.anc$g[,a2])
    c12 <- cor(eur.anc$g[,e1], eur.anc$g[,e2])
    c13 <- cor(eur.anc$g[,e1], nam.anc$g[,n2])
  } else {
    c11 <- cor(nam.anc$g[,n1], afr.anc$g[,a2])
    c12 <- cor(nam.anc$g[,n1], eur.anc$g[,e2])
    c13 <- cor(nam.anc$g[,n1], nam.anc$g[,n2])
  }
  gc()
  return(c(c11,c12,c13))
}