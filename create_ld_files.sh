# Gabriel Hoffman
# Oct 20, 2022
# 
# Compute LD for 1kg reference panels

cd /sc/arion/projects/CommonMind/hoffman/ldref

module load plink2

DIR=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/

VCF=$DIR/VCF/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Exclude SNPS with low MAF
# MAF=$DIR/MAF/EUR_AF_1000G_phase3
# cat $MAF | awk '{if($2 <0.05) print $1'

# Keep Eur samples that are not in exclude list
EXCL=$DIR/release_2013502/20140625_related_individuals.txt
INDIVS=$DIR/release_2013502/integrated_call_samples_v3.20130502.ALL.panel

# Get sample ids for each population
for POP in $(echo 'EUR' 'EAS' 'AFR');
do
	comm -23 <(cat $INDIVS | tail -n +2 | awk -vpop=$POP '{if($3 == pop) print $1, $1}' | sort) <(cut -f1 $EXCL | awk '{print $1,$1}' | sort) > ${POP}.indivs
done


# Compute LD for each population
for POP in $(echo 'CEU' 'GBR');
do
	seq 1 22 | parallel -P1 plink --threads 6 --vcf $DIR/VCF/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep ${POP}.indivs --maf 0.05 --r2 gz --ld-window-r2 0.01 --ld-window 1000 --ld-window-kb 1000 --out ld/${POP}.chr{}
done

# format genomw-wide bed
cd /sc/arion/projects/CommonMind/hoffman/ldref/ldetect-data
for POP in $(echo 'EUR' 'ASN' 'AFR');
do
tail -n +2 ${POP}/fourier_ls-all.bed | sed 's/ //g' | grep -v "None" > ${POP}/fourier_ls-all_mod.bed 
done



# Analysis based on correlation matrix
######################################

# r^2 is calculated by plink

info = data.frame(LD='ld/CEU.chr22.ld.gz', BIM='/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3//PLINK/chr22.bim')
res.CEU = readCorrMatrix(info)


info = data.frame(LD='ld/GBR.chr22.ld.gz', BIM='/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3//PLINK/chr22.bim')
res.GBR = readCorrMatrix(info)

snp.ids = intersect(rownames(res.CEU$C.ld), rownames(res.GBR$C.ld))

# check concordance in correlation
i = snp.ids[1:300]

C1 = res.CEU$C.ld[i,i]
C2 = res.GBR$C.ld[i,i]

a = C1[lower.tri(C1)]
b = C2[lower.tri(C2)]

cor(a,b)


library(decorrelate)



# from: /Users/gabrielhoffman/workspace/repos/eval_methods/decorrelate

# norm of matrix compared to identity,
# This is the root mean squared error of the off diagonal entries
normCov = function(Sigma){

  if( length(Sigma) == 1)
    if( is.na(Sigma) ) return(NA)
  # base::norm(Sigma - diag(1, nrow(Sigma)), "F")^2 / length(Sigma)

  mse = mean((Sigma-diag(1, nrow(Sigma)))^2)
  sqrt(mse)
}

ecl = eclairs_corMat( C1, n=90 )


A = decorrelate(C1, ecl)
normCov(A)
normCov(C1)


A = decorrelate(C2, ecl)
normCov(A)
normCov(C2)


library(rtracklayer)
gr = import("ldetect-data/EUR/fourier_ls-all_mod.bed", format="bed")

gr_chr = gr[seqnames(gr) == levels(seqnames(res.CEU$gr))]

# for each 1) population, 2) method, 3) chrom, method

# for each region from Pickrell BED

df_metrics = lapply(1:length(gr_chr), function(i){
	message(i)
	df_metrics = lapply(c("eb", "0.01"), function(method){
		# extract SNPs from pop1
		incl1 = res.CEU$gr %within% gr_chr[i]
		C1 = res.CEU$C.ld[incl1,incl1]

		# extract SNPs from pop2
		incl2 = res.GBR$gr %within% gr_chr[i]
		C2 = res.GBR$C.ld[incl2,incl2]

		# get intersection of SNPs
		snp.ids = intersect(rownames(C1), rownames(C2))
		C1 = C1[snp.ids,snp.ids]
		C2 = C2[snp.ids,snp.ids]

		# lambda.LW = CovTools::CovEst.2003LW( C1 )$delta
	 #    lambda.Touloumis = ShrinkCovMat::shrinkcovmat.identity(C1)$lambdahat
	 #    lambda.Schafer = corpcor::estimate.lambda(C1, verbose=FALSE)

		# perform decomposition
		ecl = eclairs_corMat( C1, n=90 )

		lambda = switch( method, eb = ecl$lambda, "0.01" = 0.01 ) 

		# whiten
		A = decorrelate(C2[snp.ids,snp.ids], ecl, lambda=lambda)

		# get rMSE
		data.frame( as.data.frame(gr_chr[i]), method, lambda, rMSE = normCov(A))
	})
	df_metrics = do.call(rbind, df_metrics)
})
df_metrics = do.call(rbind, df_metrics)

library(ggplot2)



# INSTEAD: do this directly on the data matrix read in from VCF
# since other methods start with data matrix


# Data analysis figure
# 1) show rMSE versus lambda in (0,1) and show each methods estimate of lambda
# show corr matrix between training, testing, then whitenhed

# We make the much more realistic assumption that the correlation matricies are drawn from the same distribution

# Data directly from VCF
########################

cd /sc/arion/projects/CommonMind/hoffman/ldref

module load plink2 bcftools tabix

DIR=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/

# Subset for each population, and filter by LD
for POP in $(echo 'EUR' 'ASN' 'AFR');
do
seq 1 22 | parallel -P1 plink --vcf $DIR/VCF/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep ${POP}.indivs --maf 0.05 --geno 0 --hwe 1e-5 --snps-only --recode vcf --out filter/${POP}.chr{} 
done

ls filter/*vcf | parallel bgzip -f
ls filter/*vcf.gz | parallel tabix -fp vcf


























