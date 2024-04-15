#! /usr/bin/env Rscript
library(getopt)

spec = matrix(c(
    'super_pop', 's', 1, "character",
    'out',    'o', 1, "character"
  ), byrow=TRUE, ncol=4)
opt = getopt(spec)

# SRC=/hpc/users/hoffmg01/work/decorrelate_analysis/compute_whitening_rmse.R
#
# $SRC --super_pop EUR --out df_EUR.RDS
# $SRC --super_pop AFR --out df_AFR.RDS
# $SRC --super_pop ASN --out df_ASN.RDS

# printf "EUR\nAFR\nASN" | parallel -P3 "$SRC --super_pop {} --out df_{}.RDS"


suppressPackageStartupMessages({
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(decorrelate)
library(Rfast)
library(CovTools)
library(survival)
library(Matrix)
library(corpcor)
library(RhpcBLASctl)
library(ShrinkCovMat)
})

omp_set_num_threads(4)


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

# whiten with pseudoinverse with fixed rank
get_w_ginv = function(X, k){
  C = cov(X)
  dcmp = eigen(C)
  W = with(dcmp, vectors[,seq(k)] %*% diag(1/sqrt(values[seq(k)])) %*% t(vectors[,seq(k)]))
  W
}

file = "/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/release_2013502/integrated_call_samples_v3.20130502.ALL.panel"
infoAll = read.table(file, header=TRUE)
infoAll$sample = paste(infoAll$sample, infoAll$sample, sep="_")
infoAll$super_pop[infoAll$super_pop == "EAS"] = "ASN"

methods = c('GIW-EB (current work)',  "0", "0.01", "Schäfer-Strimmer", "Ledoit-Wolf", "Touloumis", "OAS", 'Pseudoinverse') 

df_grid = expand.grid(chrom=1:22, super_pop=opt$super_pop, stringsAsFactors=FALSE)

maf = function(x){
    af = sum(x) / (2*length(x))
    min(af, 1-af)
}


df = lapply(1:nrow(df_grid), function(i){

    # subset info to this super pop
    idx = which(infoAll$super_pop == df_grid$super_pop[i])
    info = infoAll[idx,]

    # read in genome-blocks for this population
    # file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/ldetect-data/", df_grid$super_pop[i], "/fourier_ls-all_mod.bed")
    file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/adjclust/", df_grid$super_pop[i], ".bed")
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
      X_train = X[idx_train,keep]
      X_test = X[idx_test,keep]

      # learn transformation
      ecl = eclairs( X_train, compute="corr")

      rMSE_baseline = normCov(cora(X_test))

      df = lapply( methods, function(method){

          message("\r", df_grid$chrom[i], " ", k, " ", method, "   ", appendLF=FALSE)

          if( method == 'Pseudoinverse'){
            rnk = min(dim(X_train)-1)
            W = get_w_ginv(scale(X_train), rnk)

            X_test_white = tcrossprod(X_test, W)

            rMSE = normCov(cora(X_test_white))

            }else{
                X_tr = scale(X_train)

                # select lambda based on method
                lambda = switch( method, 
                  'GIW-EB (current work)' = ecl$lambda, 
                  # "beam" = beam(X_tr, verbose=FALSE)@alphaOpt,
                  "0" = 0,
                  "0.01" = 0.01,
                  "Ledoit-Wolf" = CovEst.2003LW( X_tr )$delta,
                  "OAS" = CovEst.2010OAS(X_tr)$rho,
                  "Touloumis" = shrinkcovmat.identity(X_tr)$lambdahat,
                  "Schäfer-Strimmer" = estimate.lambda(X_tr, verbose=FALSE)  )

                lambda = min(1, max(0, lambda))

                # transform testing data
                X_test_white = decorrelate(X_test, ecl, lambda=lambda)

                rMSE = normCov(cora(X_test_white))
            }

          # get rMSE
          data.frame(method, df_grid[i,], gr_chr[k], nsnps = ncol(X), lambda, rMSE, rMSE_baseline, averageCorr = averageCorr(ecl))
      })
      do.call(rbind, df)
    })
    do.call(rbind, df)
})
df = do.call(rbind, df)

saveRDS(df, file=opt$out)














