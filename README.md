# Detect genomic selection outlier loci using pcadapt [Luu et al 2016]( https://doi.org/10.1111/1755-0998.12592)

The script 00_snpsCall_haplotypeCaller.sh produces sorting of mapped files and duplicates removal, splits reads with Ns in the CIGAR string and trims overhangs. Subsequently, HaplotypeCaller of gatk v.4 (Van der Auwera and O'Connor 2020) is used to call variants.

The script 01_snpsCall_genomicDB.slrm merges multiple samples in gvcf format using the GenomicsDBImport utility.

The script 02_snpsCall_genotypeGVCF.slrm performs joint genotyping and processing of the resulting vcf file.

Filtering of SNPs was performed using vcftools v.0.1.16 (Danecek et al. 2011) with the following options:
```
vcftools \
--vcf input.vcf \
--max-alleles 2 \
--min-alleles 2 \
--minDP 4 \
--minGQ 20 \
--minQ 30 \
--remove-indels \
--max-missing 0.5 \
--recode
 ```
 
(Purcell et al. 2007) was used to convert vcf files into plink binary bed format:
 ```
 plink2 \
 --vcf input.vcf \
 --make-bed \
 --out out \
 --allow-extra-chr
 ```
 pcadapt analyses were run using the script pcadapt.r.
