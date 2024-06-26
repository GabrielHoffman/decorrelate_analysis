#! /usr/bin/env Rscript
library(getopt)

spec = matrix(c(
    'super_pop', 's', 1, "character",
    'out',    'o', 1, "character"
  ), byrow=TRUE, ncol=4)
opt = getopt(spec)

# cd /hpc/users/hoffmg01/work/decorrelate_analysis
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

# matrix square root
msqrt = function(S){
  dcmp = eigen(S)
  # with(dcmp, vectors %*% diag(sqrt(values)) %*% t(vectors))
  with(dcmp, vectors %*% (sqrt(values) * t(vectors)))
}

minvsqrt = function(S, tol = 1e-14){
  dcmp = eigen(S)  
  k = sum(dcmp$values > tol)
  # with(dcmp, vectors %*% diag(1/sqrt(values)) %*% t(vectors))
  with(dcmp, vectors[,seq(k)] %*% ((1/sqrt(values[seq(k)])) * t(vectors[,seq(k)])))
}


# whiten with pseudoinverse with fixed rank
get_w_ginv = function(X, k, tol = 1e-14){
  C = cov(X)
  dcmp = eigen(C)

  k = sum(dcmp$values[seq(k)] > tol)
  # W = with(dcmp, vectors[,seq(k)] %*% diag(1/sqrt(values[seq(k)])) %*% t(vectors[,seq(k)]))
  W = with(dcmp, vectors[,seq(k)] %*% ((1/sqrt(values[seq(k)])) * t(vectors[,seq(k)])))
  W
}


file = "/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/release_2013502/integrated_call_samples_v3.20130502.ALL.panel"
infoAll = read.table(file, header=TRUE)
infoAll$sample = paste(infoAll$sample, infoAll$sample, sep="_")
infoAll$super_pop[infoAll$super_pop == "EAS"] = "ASN"

df_grid = expand.grid(chrom=1:22, super_pop=opt$super_pop, stringsAsFactors=FALSE)

maf = function(x){
    af = sum(x) / (2*length(x))
    min(af, 1-af)
}


# subset info to this super pop
idx = which(infoAll$super_pop == opt$super_pop)
info = infoAll[idx,]

# get training set
set.seed(1)
idx_train = sample(nrow(info), 0.5*nrow(info))
idx_test = setdiff(seq(nrow(info)), idx_train)

df = lapply(1:nrow(df_grid), function(i){

    # read in genome-blocks for this population
    file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/ldetect-data/", opt$super_pop, "/fourier_ls-all_mod.bed")
    # file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/adjclust/", opt$super_pop, ".bed")
    gr = import(file, format="bed")
    seqlevelsStyle(gr) = "NCBI"

    # subset Granges for this chrom
    gr_chr = gr[seqnames(gr) == df_grid$chrom[i]]

    df = lapply(seq(length(gr_chr)), function(k){
      # Read data in range
      vcf.file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/filter/", opt$super_pop, ".chr",df_grid$chrom[i], ".vcf.gz")
      res = readVcf( vcf.file, genome = "GRCh37", param = gr_chr[k] )
      data = genotypeToSnpMatrix(res)

      # Convert to numeric
      X = as(data$genotypes, "numeric")

      identical(rownames(X), info$sample)

      # pass MAF filter in both sets
      min_af = 0
      keep = (apply(X[idx_train,], 2, maf) > min_af) & (apply(X[idx_test,], 2, maf) > min_af)
      table(keep)

      Y = scale(X)

      df = data.frame()

      # Eval methods
      #-------------

      # decorrelate
      tm = system.time({
      ecl <- eclairs( Y[idx_train,])
      y_white <- decorrelate(Y[-idx_train,], ecl)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = 'GIW-EB', 
                        t(c(tm)),
                        rmse = rmse))

      # decorrelate
      tm = system.time({
      ecl <- eclairs( Y[idx_train,], k=50)
      y_white <- decorrelate(Y[-idx_train,], ecl)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                  Method = 'GIW-EB (k=50)', 
                  t(c(tm)),
                  rmse = rmse))

      tm = system.time({
      ecl <- eclairs( Y[idx_train,])
      y_white <- decorrelate(Y[-idx_train,], ecl, lambda = 0)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "lambda = 0",  
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      ecl <- eclairs( Y[idx_train,])
      y_white <- decorrelate(Y[-idx_train,], ecl, lambda = 0.01)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "lambda = 0.01", 
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      ecl <- eclairs( Y[idx_train,])
      y_white <- decorrelate(Y[-idx_train,], ecl, lambda = 1e-4)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "lambda = 1e-4",  
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      fit <- CovTools::CovEst.2003LW( Y[idx_train,] )
      y_white <- Y[-idx_train,] %*% minvsqrt(fit$S)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "Ledoit-Wolf",  
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      fit <- CovTools::CovEst.2010OAS( Y[idx_train,] )
      y_white <- Y[-idx_train,] %*% minvsqrt(fit$S)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "OAS",  
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      fit <- shrinkcovmat.equal( t(Y[idx_train,]) )
      y_white <- Y[-idx_train,] %*% minvsqrt(fit$Sigmahat)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "Touloumis", 
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      C <- cor.shrink( Y[idx_train,], verbose=FALSE )
      y_white <- Y[-idx_train,] %*% minvsqrt(C)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "Schäfer-Strimmer", 
                        t(c(tm)),
                        rmse = rmse))

      tm = system.time({
      # learn transformation
      k <- min(dim(Y[idx_train,]))-1
      W <- get_w_ginv(Y[idx_train,], k)
      y_white <- tcrossprod(Y[-idx_train,], W)
      })   
      rmse = normCov(cora(y_white))
      df = rbind(df, data.frame(
                        Method = "Pseudoinverse", 
                        t(c(tm)),
                        rmse = rmse))    

      df$rMSE_baseline = normCov(cora(Y[-idx_train,]))
      df$chrom = df_grid$chrom[i]
      df$super_pop = df_grid$super_pop[i]

      cat(df_grid$chrom[i], k, "  \r")

      cbind(df, nsnps = ncol(Y))
    })
    do.call(rbind, df)
})
df = do.call(rbind, df)

saveRDS(df, file=opt$out)














