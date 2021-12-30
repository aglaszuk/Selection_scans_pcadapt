#!/bin/bash

###########  Script to call variants using GATK V.4  ###########

# Define files and program variables needed for the analysis

wd=/path/to/the/starting/working/dir/where/bamfiles/are/located
ref=/path/to/genome.scf.fasta
ref_dir=/path/to/genome/directory
picard=/path/to/picard.jar

# Sort files mapped to reference genome using picard v 2.23.4:

cd $wd

mkdir sorted

for file in *.bam; \
do \
java -Xmx32g \
-jar $picard \
AddOrReplaceReadGroups \
I=$file \
O=./sorted/${file/.bam/.sorted.bam} \
SORT_ORDER=coordinate \
RGID=${file/.bam/} \
RGLB=${file/.bam/} \
RGPL=illumina \
RGPU=machine \
RGSM=${file/.bam/}; \
done

cd sorted

# Mark and remove duplicates in sorted files and create index using picard v. 2.23.4

mkdir deduplicated

for file in *.bam; \
do \
java -Xmx32g \
-jar $picard \
MarkDuplicates \
I=$file \
O=./deduplicated/${file/.sorted.bam/.dedupli.bam} \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=$file.metrics; \
done

cd ./deduplicated/

# split'N'trim using GATK v.4.1.8.1

java -Xmx32g -jar $picard \
CreateSequenceDictionary \
-R $ref 

samtools faidx $ref 

mkdir split

for file in *.bam; \
do \
$gatk \
SplitNCigarReads \
-R $ref \
-I $file \
-O ./split/${file/.dedupli.bam/.split.bam}; \
done

# Call variants using HaplotypeCaller

cd split

mkdir gvcf

for file in *.bam; \
do \
$gatk \
—java-options "-Xmx32g" \
HaplotypeCaller \
-R $ref \
-I ./$file \
-O ./gvcf/${file/.split.bam/.g.vcf} \
-ERC GVCF; \
done

# Create genome subsets to parallelize genomeDB jobs in the next step

cd $ref_dir/

/apps/freebayes/1.3.2/scripts/fasta_generate_regions.py \
$ref \
100000 > \
regions.txt

mkdir genome_chunks

split \
-l 1000 \
-d regions.txt \
genome_chunks/chunks

for file in genome_chunks/*; \
do \
mv $file genome_chunks/$file.bed; \
done
