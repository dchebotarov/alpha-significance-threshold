# Usage:
# Rscript num-effective-tests.R -b mydata.bed -w 800000
# (this will use PLINK binary file mydata.bed/bim/fam and 800kb window size
# Output:
# file listing the estimates of number of effective SNPs according to 3 methods:
# Cheverud
# Li and Ji 2005
# Gao et al 2009




library(getopt)
suppressPackageStartupMessages( library(snpStats))

if(!interactive()){
  # 1 - required, 2 - optional
  opt = getopt(spec = matrix(c(
    "gwas_result", "f", 1, "character",
    "plink_root", "b", 1, "character",
    "window_size", "w", 2, "numeric"
   # "annot_file",  "a", 2, "character"
  ), byrow = TRUE, ncol=4))


  gwas_result_file = opt$gwas_result
  if(is.null(gwas_result_file)){
    gwas = NULL
    #stop("Please specify GWAS result file as returned by GEMMA (columns chr,rs,ps,p_wald) using -f argument")
  }

  plink_root_file = opt$plink_root

  window_size = opt$window_size
  if(is.null(window_size)){
    window_size = 5e5
  }

} else {
  # replace with your own values
  #gwas_result_file = "output/lmm-SeedlHeight.assoc.txt"

}


# cumsum_thresh = 0.995  # for Gao
cumsum_thresh = 0.99

library(readr)
suppressPackageStartupMessages(library(dplyr))

X = snpStats::read.plink(opt$plink_root)
names(X) = c("bed", "fam", "bim")

if(!is.null(gwas_result_file)){
  gw = read_table2(gwas_result_file)

  if(!("P" %in% names(gw) ) & (! "p_wald" %in% names(gw)) ){
    gw = read_table2(gwas_result_file)
  }
  if("P" %in% names(gw)){
    names(gw)[ names(gw)=="P"]="p_wald"
  }

  if("EMP1" %in% names(gw)){
    names(gw)[ names(gw)=="EMP1"]="p_wald"
  }

  gw$logp = -log10(gw$p_wald)
  if("CHR" %in% names(gw)){
    names(gw)[ names(gw)=="CHR"]="chr"
  }
  if("BP" %in% names(gw)){
    names(gw)[ names(gw)=="BP"]="ps"
  }
  if("SNP" %in% names(gw)){
    names(gw)[ names(gw)=="SNP"]="rs"
  }
  gwas=gw

  gwas = gwas[ !is.na(gwas$p_wald) ,]

  gwas$chr = as.integer(gwas$chr)
  gwas$chr[ is.na(gwas$chr) ]  = 0

  if(nrow(gwas)==0){
    stop("No SNPs remain!")
  }

  gwas$logP = -log10(gwas$p_wald)

  used_in_gwas = intersect( X$bim[,2], gwas$rs) 
  if(length(used_in_gwas)==0){
    stop("No SNPs common between GWAS file and plink file ! \n")
  } else {
    cat("SNPs both in GWAS and plink file: ", length(used_in_gwas), "\n")
  }

  X$bed = X$bed[, X$bim[,2] %in% used_in_gwas  ]
  X$bim = X$bim[ X$bim[,2] %in% used_in_gwas, ]

}

win_size = window_size
win_step = window_size

NEFF = data.frame(
  chr = 1:12,
  Meff  = numeric(12),
  Neff_LJ = numeric(12),
  Meff_Gao = integer(12)
  )
win_df_list = list()

for(focus_chr in 1:12){
  cat("###  Chromosome" , focus_chr, " ###\n")

  select_chr = X$bim[,1] == focus_chr

  Xchr = list(bed = X$bed[, select_chr],
              fam = X$fam,
              bim = X$bim[ select_chr, ]
  )

  win_df = data.frame(
    start = seq(1,
                max(Xchr$bim[,4]) -win_size + win_step,
                by=win_step )
  )
  win_df$end = win_df$start + win_size -1
  win_df$win_ix = 1:nrow(win_df)

  win_df$meanR2 = NA_real_
  win_df$medianR2 = NA_real_
  win_df$propLD = NA_real_

  win_df$nsnp = NA_integer_

  win_df$Meff_Cheverud = NA_real_
  win_df$Meff_LJ = NA_real_
  win_df$Meff_Gao = NA_integer_

  # stopifnot( all(gwas$rs == X$bim[,2]) ) # TRUE

  for( i  in 1:dim(win_df)[[1]]) {
    cur_win = win_df$win_ix[[i]]
    cur_win_start = win_df$start[[i]] # win_df$win_ix == cur_win]
    cur_win_end = win_df$end[[i]] # win_df$win_ix == cur_win]

    chr_snp_selection = Xchr$bim[,4] <= cur_win_end &
      Xchr$bim[,4] >= cur_win_start

    win_df$nsnp[[i]] = sum(chr_snp_selection)
    if(sum(chr_snp_selection) < 2) next
    cat("Window ", i, ": ", sum(chr_snp_selection), "snps.  ")

    regX = Xchr$bed[, chr_snp_selection  ]
    regX = as(regX, "numeric")
    means_regX = colMeans(regX, na.rm=T)
    which_na = which(is.na(regX))
    if(length(which_na)>0){
      na_replace = rep( means_regX, each=nrow(regX))[ which_na]
      regX[ which_na ] = na_replace
    }

    C = cor(regX, use="pairw")
    #C2 = cor(rX, use="pairw")

    if(any(is.na(C))){
      browser()
    }
    EC = eigen(C, symmetric = T)
    eigvalues = EC$values
    cat(" min eigv=", format(min(eigvalues), digits=2), "  ")


    nn = sum(eigvalues)
    #cat("sum(eigv)=",nn,"  " )
    eigvalues[ eigvalues <  0]  = 0

    nn2 = sum(eigvalues)
    # Gao et al
    cumsum_eig = cumsum(eigvalues)
    Meff_Gao = which( cumsum_eig >= cumsum_thresh * nn )[[1]]
    Meff_Gao_995 = which( cumsum_eig >= 0.995 * nn2 )[[1]]
    Meff_Gao_99 = which( cumsum_eig >= 0.99 * nn2 )[[1]]

    # For other methods, the small values don't matter
    eigvalues[ eigvalues < 1e-5] = 0
    cat("sum(eigv)=", format(nn2, digits = 1) ,"  " )

    M = dim(C)[[1]]  # same as nsnp[[i]]
    V_lambda =sum ( (eigvalues-1)^2 / (M-1))
    # Cheverud:
    Meff = 1 + (M-1)*(1- V_lambda/M)
    # Li and Ji:
    Meff_LJ = sum( eigvalues - floor(eigvalues)) + sum(eigvalues >= 1)

    cat(" LJ=", format(Meff_LJ, digits = 3), " ")

    cat(" Gao=", Meff_Gao, " ")
    cat(" Gao99=", Meff_Gao_995, " ")
    cat(" Gao995=", Meff_Gao_995, " ")
    cat("\n")

    win_df$Meff_Cheverud[[i]] = Meff
    win_df$Meff_LJ[[i]] = Meff_LJ
    win_df$Meff_Gao[[i]] =  Meff_Gao
  }
  NEFF$Meff[[focus_chr]] = sum(win_df$Meff_Cheverud, na.rm = T)
  NEFF$Neff_LJ[[focus_chr]] = sum(win_df$Meff_LJ, na.rm=T)
  NEFF$Meff_Gao[[focus_chr]] = sum(win_df$Meff_Gao, na.rm=T)

  win_df_list[[ focus_chr ]] = win_df
}

write_csv(NEFF, paste0("Neff_windsize", win_size/1e6,
                       "Mb_step", win_step/1e6,
                       "Mb.csv"))

ans = data.frame(
  tot_Neff_Cheverud =  round(sum(NEFF$Meff), 2),
  tot_Neff_LJ =  round(sum(NEFF$Neff_LJ), 3),
  tot_Neff_Gao = sum(NEFF$Meff_Gao)
)

write_tsv(ans, paste0("Num_eff_tests_estimates-window", win_size/1e6,  "Mb.csv" ))
cat(" The results are saved. Significance level thresholds (-log10(p) scale) according to (Li and Ji 2005) : \n")
cat(" - for alpha=0.1: \n")
print( -log10(0.1/ans$tot_Neff_LJ) )
cat(" - for alpha=0.05: \n")
print( -log10(0.05/ans$tot_Neff_LJ) )
cat("A possible suggestive threshold (e.g. alpha=0.5): \n")
print( -log10(0.5/ans$tot_Neff_LJ) )

# -log10( 0.05/tot_Neff_LJ)  # 5.53


