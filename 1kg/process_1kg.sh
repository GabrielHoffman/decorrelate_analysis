


# Create filtered VCFs
ml plink/1.90b6.5 parallel
cd /sc/arion/projects/CommonMind/hoffman/ldref
seq 1 22 | parallel -P8 "plink  --geno 0  --hwe 1e-5  --keep ASN.indivs --maf 0.05 --out filter/ASN.chr{} --recode vcf --snps-only --vcf /sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3//VCF/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"



# Compute LD
ml plink/1.90b6.5 parallel
cd /sc/arion/projects/CommonMind/hoffman/ldref/filter

# write LD files
ls *.vcf.gz | sed 's/.vcf.gz//g' | parallel -P24 "plink --threads 1 --vcf {}.vcf.gz -r2 gz --ld-window 300 --ld-window-r2 .01 --out {}"

# Write BIM files
ls *.vcf.gz | sed 's/.vcf.gz//g' | parallel -P24 "plink --threads 1 --vcf {}.vcf.gz --make-bed --out {}"


# Compute LD clusters
######################
cd /sc/arion/projects/CommonMind/hoffman/ldref/filter

source("/hpc/users/hoffmg01/work/decorrelate_analysis/1kg/readLDMatrix.R")
library(adjclust)
library(parallel)

pop_array = c("EUR", "AFR", "ASN")
hcl.lst = mclapply(pop_array, function(pop){

	# for each chromosome
	hcl.lst = mclapply(seq(22), function(chr){
		id = paste0(pop, '.chr', chr, '.')
		df_files = data.frame(chrom = chr, 
				BIM = paste0(id, "bim"), 
				LD = paste0(id, "ld.gz"))

		# read LD matrix		
		res.ld = readLDMatrix(df_files)

		# Compute adjacency clustering
		hcl = adjClust(res.ld[[1]]$C.ld, h=300-1, strictCheck = FALSE)

		# add info
		hcl$chr = chr
		hcl$pop = pop
		hcl
	}, mc.cores=22)
	names(hcl.lst) = paste0("chr", seq(22))

	hcl.lst
}, mc.cores=3)
names(hcl.lst) = pop_array

# saveRDS("hcl.lst", file="hcl.lst_1kg.RDS")

# Cut into clusters
mean_size = 2000
res.clusters = mclapply( names(hcl.lst), function(pop){
	res = mclapply(hcl.lst[[pop]], function(hcl){
		k = length(hcl$labels) / mean_size

		# cut tree into clusters
		clID = cutree_chac(hcl, as.integer(k))

		# rename clusters
		clID.return = paste(hcl$pop, hcl$chr, clID, sep="_")
		names(clID.return) = names(clID)

		data.frame(pop = hcl$pop, 
					chr = hcl$chr, 
					name = names(clID.return),
					cluster = clID.return)
	}, mc.cores=22)
	res = do.call(rbind, res)
	rownames(res) = c()
	res
}, mc.cores=3)
names(res.clusters) = names(hcl.lst)

# write clusters to bed file
res_bed = mclapply( names(hcl.lst), function(pop){
	files = paste0(pop, '.chr', seq(22), '.', "bim")
	df = mclapply(files, fread, mc.cores=4)
	df = do.call(rbind, df)
	colnames(df) = c("chr", "name", "V3", "position", "A1", "A2")
	df = df[! duplicated(df$name),]

	# get clusters for this pop
	df_clust = res.clusters[[pop]]

	# intersect names
	shared = intersect(df$name,df_clust$name)
	df = df[df$name %in% shared,]
	df_clust = df_clust[df_clust$name %in% shared,]

	# assign SNPs to clusters
	i = match(df_clust$name, df$name)
	df$cluster[i] = df_clust$cluster

	res = mclapply(unique(df$cluster), function(id){

		df_sub = df[df$cluster == id,]

		data.frame(chr = paste0("chr", df_sub$chr[1]), 
			start = min(df_sub$position), 
			stop = max(df_sub$position))
	}, mc.cores = 12)
	do.call(rbind, res)
}, mc.cores=3)
names(res_bed) = names(hcl.lst)


lapply( names(res_bed), function(pop){
	file = paste0("/sc/arion/projects/CommonMind/hoffman/ldref/adjclust/", pop, ".bed")
	write.table(res_bed[[pop]], file=file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
} )















