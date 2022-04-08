#!/bin
# Data preparation before R analysis

# Download data from ENA - fastq
wget \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR066/SRR066787/SRR066787.fastq.gz \
-O SRR066787_wce.fastq.gz
wget \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR066/SRR066766/SRR066766.fastq.gz \
-O SRR066766_H3K27AC.fastq.gz
wget \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR066/SRR066767/SRR066767.fastq.gz \
-O SRR066767_H3K27AC.fastq.gz

#-----(R ShortRead package for QC)

# Download genome annotation from iGenomes
# https://support.illumina.com/sequencing/sequencing_software/igenome.html

# Download genomic data with Bowtie2 formatting
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
mkdir indexes
mkdir indexes/mm10
mv mm10.zip indexes/mm10
unzip indexes/mm10/mm10.zip

# Run Bowtie
#Unpack the fastq archives
gzip –d SRR*gz
#Run Bowtie to align to the genome
bowtie2 -p 20 \
-q /home/guest101/chip_seq/indexes/mm10/mm10SRR066787_wce.fastq \
-S ES_input.sam
bowtie2 -p 20 \
-q /home/guest101/chip_seq/indexes/mm10/mm10SRR066766_H3K27AC.fastq \
-S H3K27ac_rep1.sam
bowtie2 -p 20 \
-q /home/guest101/chip_seq/indexes/mm10/mm10SRR066767_H3K27AC.fastq \
-S H3K27ac_rep2.sam

# Samtools select hq scoring alignments
samtools view -bS -q 40 ES_input.sam > ES_input_bestAlignment.bam
samtools view -bS -q 40 H3K27ac_rep1.sam > H3K27ac_rep1_bestAlignment.bam
samtools view -bS -q 40 H3K27ac_rep2.sam > H3K27ac_rep2_bestAlignment.bam

# Remove PCR duplicates (overrepped fragments)
samtools rmdup -s ES_input_bestAlignment.bam ES_input_filtered.bam
samtools rmdup -s H3K27ac_rep1_bestAlignment.bam H3K27ac_rep1_filtered.bam
samtools rmdup -s H3K27ac_rep2_bestAlignment.bam H3K27ac_rep2_filtered.bam
# or use Picard function markduplicates

# Convert to BED format
bedtools bamtobed -i ES_input_filtered.bam > ES_input_filtered.bed
bedtools bamtobed -i H3K27ac_rep1_filtered.bam > H3K27ac_rep1_filtered.bed
bedtools bamtobed -i H3K27ac_rep2_filtered.bam > H3K27ac_rep2_filtered.bed

#*Analyse Chromosome 6 genes only
awk '{if($1=="chr6") print $0}' ES_input_filtered.bed > ES_input_filtered_ucsc_chr6.bed
awk '{if($1=="chr6") print $0}' H3K27ac_rep1_filtered.bed > H3K27ac_rep1_filtered_ucsc_chr6.bed
awk '{if($1=="chr6") print $0}' H3K27ac_rep2_filtered.bed > H3K27ac_rep2_filtered_ucsc_chr6.bed

#-----(R build annotations)

# Peak finding using MACS
# Two replicates and one control
macs14 -t H3K27ac_rep1_filtered.bed \
-c ES_input_filtered.bed \
-f BED \
-g mm \
--nomodel \
-n Rep1
macs14 -t H3K27ac_rep2_filtered.bed \
-c ES_input_filtered.bed \
-f BED \
-g mm \
--nomodel \
-n Rep2

#*Filter by Chr6 data only
awk '{if($1=="chr6") print $0}' Rep1_peaks.bed > Rep1_peaks_ucsc_chr6.bed
awk '{if($1=="chr6") print $0}' Rep2_peaks.bed > Rep2_peaks_ucsc_chr6.bed