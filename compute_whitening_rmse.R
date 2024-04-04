#! /usr/bin/env Rscript
library(getopt)

spec = matrix(c(
    'super_pop', 's', 1, "character",
    'out',    'o', 1, "character"
  ), byrow=TRUE, ncol=4)
opt = getopt(spec)

# /hpc/users/hoffmg01/work/decorrelate_analysis/compute_whitening_rmse.R --super_pop EUR --out df_EUR.RDS

suppressPackageStartupMessages({
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(decorrelate)
library(Rfast)
library(ggplot2)
library(CovTools)
library(corpcor)
library(beam)
library(ShrinkCovMat)
library(Matrix)
library(survival)
library(tidyverse)
library(whitening)
})


# from: /Users/gabrielhoffman/workspace/repos/eval_methods/decorrelate
# norm of matrix compared to identity,
# This is the root mean squared error of the off diagonal entries
normCov = function(Sigma){

  if( length(Sigma) == 1)
  if( is.na(Sigma) ) return(NA)

  # Based on Frobenius norm  
  # same as mse
  # value = base::norm(Sigma - diag(1, nrow(Sigma)), "F")^2 / length(Sigma)
  # sqrt(value)
    # base::norm(Sigma - diag(1, nrow(Sigma)), "F") / nrow(Sigma)

  mse = mean((Sigma-diag(1, nrow(Sigma)))^2)
  sqrt(mse)
}


file = "/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/release_2013502/integrated_call_samples_v3.20130502.ALL.panel"
infoAll = read.table(file, header=TRUE)
infoAll$sample = paste(infoAll$sample, infoAll$sample, sep="_")


methods = c("eb",  "0", "0.01", "Schafer", "LW", "Touloumis") 

df_grid = expand.grid(chrom=22:22, super_pop=opt$super_pop, stringsAsFactors=FALSE)

maf = function(x){
    af = sum(x) / (2*length(x))
    min(af, 1-af)
}


df = lapply(1:nrow(df_grid), function(i){

    # subset info to this super pop
    idx = which(infoAll$super_pop == df_grid$super_pop[i])
    info = infoAll[idx,]

    # read in genome-blocks for this population
    file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/ldetect-data/", df_grid$super_pop[i], "/fourier_ls-all_mod.bed")
    gr = import(file, format="bed")
    seqlevelsStyle(gr) = "NCBI"

    # subset Granges for this chrom
    gr_chr = gr[seqnames(gr) == df_grid$chrom[i]]

    # get training set
    idx_train = sample(nrow(info), 0.5*nrow(info))
    idx_test = setdiff(seq(nrow(info)), idx_train)

    df = lapply(seq(length(gr_chr)), function(k){
      # Read data in range
      vcf.file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/filter/", df_grid$super_pop[i], ".chr",df_grid$chrom[i], ".vcf.gz")
      res = readVcf( vcf.file, genome = "GRCh37", param = gr_chr[k] )
      data = genotypeToSnpMatrix(res)

      # Convert to numeric
      X = as(data$genotypes, "numeric")

      identical(rownames(X), info$sample)

      # pass MAF filter in both sets
      min_af = 0
      keep = (apply(X[idx_train,], 2, maf) > min_af) & (apply(X[idx_test,], 2, maf) > min_af)
      table(keep)

      X = scale(X)
      # Get SNPs with non-zero variance in both training and testing
      X_train = X[idx_train,]
      X_test = X[idx_test,]

      # learn transformation
      ecl = eclairs( X_train[,keep], compute="corr")

      rMSE_baseline = normCov(cora(X_test[,keep]))

      df = lapply( methods, function(method){

          message("\r", df_grid$chrom[i], " ", k, " ", method, "   ", appendLF=FALSE)

          # select lambda based on method
          lambda = switch( method, 
            eb = ecl$lambda, 
            "beam" = beam(X_train[,keep], verbose=FALSE)@alphaOpt,
            "0" = 0,
            "0.01" = 0.01,
            "LW" = CovEst.2003LW( scale(X_train[,keep]) )$delta,
            "Touloumis" = shrinkcovmat.identity(scale(X_train[,keep]))$lambdahat,
            "Schafer" = estimate.lambda(scale(X_train[,keep]), verbose=FALSE)  )

          lambda = min(1, max(0, lambda))

          # transform testing data
          X_test_white = decorrelate(X_test[,keep], ecl, lambda=lambda)

          rMSE = normCov(cora(X_test_white))

          # get rMSE
          data.frame(method, df_grid[i,], gr_chr[k], nsnps = ncol(X), lambda, rMSE, rMSE_baseline, averageCorr = averageCorr(ecl))
      })
      do.call(rbind, df)
    })
    do.call(rbind, df)
})
df = do.call(rbind, df)

saveRDS(df, file=opt$out)

