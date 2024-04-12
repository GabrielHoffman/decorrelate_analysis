


# Create filtered VCFs
ml plink/1.90b6.5 parallel
cd /sc/arion/projects/CommonMind/hoffman/ldref
seq 1 22 | parallel -P8 "plink  --geno 0  --hwe 1e-5  --keep ASN.indivs --maf 0.05 --out filter/ASN.chr{} --recode vcf --snps-only --vcf /sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3//VCF/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"



# Compute LD
ml plink/1.90b6.5 parallel
cd /sc/arion/projects/CommonMind/hoffman/ldref/filter

# write LD files
ls *.vcf.gz | head -n 4 | sed 's/.vcf.gz//g' | parallel -P24 "plink --threads 1 --vcf {}.vcf.gz -r2 gz --ld-window 200 --ld-window-r2 .01 --out {}"

# Write BIM files

ls *.vcf.gz | head -n 4 | sed 's/.vcf.gz//g' | parallel -P24 "plink --threads 1 --vcf {}.vcf.gz --make-bed --out {}"





df_files = data.frame(chrom = 10, BIM = "AFR.chr10.bim", LD = "AFR.chr10.ld.gz")

res = readLDMatrix(df_files)