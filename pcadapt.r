#! /usr/bin/env Rscript

# Author: Aglaia Szukala, aglaia.szukala [at] univie.ac.at
# Title: selectOutlierSNPs.r

# Usage: pcadapt.r input.bed input.vcf outdir K pop_names.txt pch.txt
# input.bed is a .bed file that was produced using plink [plink2 --vcf input.vcf --make-bed --out out --allow-extra-chr]
# K=number of populations in dataset. The choice of K should be checked looking at the output plot "Variance_explained.pdf"
# pop_names.txt and pch.txt are tables with population and pch (plot) ecotype identifiers respectively

# Load required packages
library("devtools")
library("qvalue")
library("ggplot2")
library("pcadapt")
library(vcfR)

args<-commandArgs(trailingOnly=T)
outdir <- args[3] #output directory

# Read in population identifiers
pop <- read.table(args[5])
pop_names <- paste0(pop[,1])

# Read genotype matrix in .bed format (plink)  into pcadapt
file_pop <- read.pcadapt(args[1], type = "bed")

# Choose number of K (components), which is big enough
K <- 20

# Find number of components 
x <- pcadapt(file_pop, min.maf = 0.05, K = K) #threshhold of maf = 0.05 use opt min.maf to change it
pdf(paste0(outdir,"Variance_explained.pdf"))
plot(x, option = "screeplot") #plot Proportion of Variance explained by each K
dev.off()

# Plot PCA, Manhattan plot and qqplot
locs_pch <- read.table(args[6])
locs_pch <- locs_pch[,1]

pdf(paste0(outdir,"PC1_2.pdf")) 
plot(x, option = "scores", pop = pop_names)+
  labs(x = paste0("PC1 (", round((x$singular.values[1])^2*100,2),"%)"), 
       y = paste0("PC2 (", round((x$singular.values[2])^2*100,2),"%)"))+
  geom_point(aes(fill=pop_names),pch=locs_pch, cex=4,show.legend = F)+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(paste0(outdir,"PC3_4.pdf"))
plot(x, option = "scores", i = 3, j = 4, pop = pop_names)+ #use i and j options to plot other components
  labs(x = paste0("PC3 (", round((x$singular.values[3])^2*100,2),"%)"),
       y = paste0("PC4 (", round((x$singular.values[4])^2*100,2),"%)"))+
  geom_point(aes(fill=pop_names), pch=locs_pch, cex=4,show.legend = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# Define K as number of demes
K <- as.integer(args[4]) #use K defined by user to perform downstream analyses

# Check LD: display the loadings and evaluate if the loadings are clustered in a single or several genomic regions
pdf(paste0(outdir,"LD_check.pdf"))
par(mfrow = c(2, 1))
for (i in 1:K)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
dev.off()
par(mfrow = c(1, 1))

# Thinning in case you observe LD in certain regions
# If you see that one or more components are determined by a single genomic region, you should thin the SNPs
x <- pcadapt(file_pop, K = K, min.maf = 0.05, 
             LD.clumping = list(size = 200, thr = 0.1), method = "componentwise")

pdf(paste0(outdir,"Q-Qplot_PC1.pdf"))
plot(x, option = "qqplot", threshold = 0.1, K=1, main ="First component")
dev.off()
pdf(paste0(outdir,"Q-Qplot_PC2.pdf"))
plot(x, option = "qqplot", threshold = 0.1, K=2, main ="Second component")
dev.off()

pdf(paste0(outdir,"Stat.distribution_PC1.pdf")) #The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷
plot(x, option = "stat.distribution", K=1)
dev.off()

pdf(paste0(outdir,"Stat.distribution_PC2.pdf")) #The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷
plot(x, option = "stat.distribution", K=2)
dev.off()

# Detect Fst outliers using pcadapt
# Compute Bonferroni adjusted p-values and define ouliers (most conservative)
x <- pcadapt(file_pop, K = K, min.maf = 0.05, 
             LD.clumping = list(size = 200, thr = 0.1))
x$padj <- p.adjust(x$pvalues, method = "bonferroni")
alpha <- 0.05
outliers_pcadapt <- which(x$padj < alpha) #returns the positions of the SNPs
N_outliers <- length(outliers_pcadapt)

# Select outliers along the first PC component, which separates the ecotypes
snp_pc <- get.pc(x, outliers_pcadapt)
out_snps <- as.list(snp_pc[snp_pc$PC == 1,][1])
out_snps <- out_snps$SNP
N_outliers_pc1 <- length(out_snps)

# Output file with outlier SNPs 
write.table(print(out_snps),
            file = paste0(outdir,N_outliers_pc1,"outliers_pcadapt.txt"),
            sep = "\t", quote = F)

# Output file with outlier SNPs coordinates
vcf1<-read.vcfR(args[2]) # read in vcf file
vcfann<-as.data.frame(getFIX(vcf1)) # converting to data frame
write.table(vcfann[rownames(vcfann) %in% out_snps,], 
            file = paste0(outdir,"outliers_pcadapt_coordinates.txt"), 
            sep = "\t", quote = F, row.names = F)

# Plot Manhattan with outlier
x$chrom <- vcfann$CHROM
x$Tag = ifelse(as.numeric(x$chrom) %% 2 == 0, 0, 1)
x$Tag[x$padj < 0.05] <- "firebrick"
x$Tag[x$Tag == "0"] <- "black"
x$Tag[x$Tag == "1"] <- "black" #change to "grey" if you want to colour by chromosome
tag <- x$Tag[x$pass]

h = -(log10(max(x$pvalues[which(x$padj < alpha)])))

pdf(paste0(outdir,"Manhattan.pdf"))
plot(x, option = "manhattan") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(aes(fill=tag), 
             col=tag,
             pch=16, cex=0.95,
             show.legend = F) +
  geom_hline(aes(yintercept=h),
             linetype="dashed", color = "firebrick")
dev.off()
