# Bioinformatics
Bioinformatics - NGS, NCBI, BWA, SAMtools, Picard, GATK, & more

# Segatella (Prevotella) copri Genome Analysis Project

## Table of Contents
1. Introduction
2. Background and Concepts
3. Workflow Overview
4. Detailed Methodology
5. Code and Commands
6. Interpretation of Results
7. Future Applications
8. Conclusion

## 1. Introduction

This project aims to analyze the genome of Segatella (Prevotella) copri using next-generation sequencing (NGS) data. The workflow includes quality control, read trimming, alignment to a reference genome, and variant calling to identify genetic variations within the bacterial genome.

## 2. Background and Concepts

### DNA Structure and Sequencing
DNA (deoxyribonucleic acid) is composed of four types of nucleotides: Adenine (A), Thymine (T), Cytosine (C), and Guanine (G). These nucleotides form base pairs (A-T and C-G) in the double-helix structure of DNA. DNA sequencing determines the order of these nucleotides.

### Sequencing Reads
Modern sequencing technologies produce short fragments of DNA sequences called "reads". This project uses paired-end sequencing, generating two reads for each DNA fragment:
- Forward Read (R1): Sequenced from one end of the fragment
- Reverse Read (R2): Sequenced from the opposite end of the same fragment

### Quality Scores
Each base in a sequencing read is assigned a quality score (Phred score), representing the probability of an incorrect base call. Higher scores indicate higher confidence.

### Trimming and Quality Control
Raw sequencing data often contains low-quality bases and adapter sequences. Trimming removes these, resulting in:
- Paired reads: Both forward and reverse reads pass quality thresholds
- Unpaired reads: Only one read of a pair passes quality thresholds

### Alignment
Trimmed reads are aligned to a reference genome, mapping each read to its most likely origin in the genome.

### Variant Calling
After alignment, genetic variants (differences between sequenced DNA and the reference genome) are identified, including:
- Single Nucleotide Polymorphisms (SNPs): Single base changes
- Insertions and Deletions (Indels): Addition or removal of bases

## 3. Workflow Overview

1. Quality Control (FastQC)
2. Trimming (Trimmomatic)
3. Alignment (BWA)
4. Post-alignment Processing (Samtools, Picard)
5. Variant Calling (GATK HaplotypeCaller)

## 4. Detailed Methodology

### 4.1 Quality Control
FastQC is used to assess the quality of raw sequencing reads, providing insights into sequence quality, GC content, and overrepresented sequences.

### 4.2 Trimming
Trimmomatic removes low-quality bases and adapter sequences, improving the overall quality of the dataset.

### 4.3 Alignment
BWA (Burrows-Wheeler Aligner) maps the trimmed reads to the Segatella copri reference genome.

### 4.4 Post-alignment Processing
Samtools and Picard tools are used to sort, index, and mark duplicate reads in the aligned data.

### 4.5 Variant Calling
GATK HaplotypeCaller identifies genetic variants by comparing the aligned reads to the reference genome.

## 5. Code and Commands

### Environment Setup
conda create -n bioinfo_env python=3.8
conda activate bioinfo_env
conda install -c bioconda bwa samtools gatk4 fastqc trimmomatic picard

### Download Reference Genome and Sequence Reads
datasets download genome accession GCF_020735445.1 --filename GCF_020735445.1.zip
unzip GCF_020735445.1.zip
prefetch SRX23539992
fastq-dump --split-files SRR27878129

### Quality Control
fastqc /home/youngdeblynn/project/raw_data/SRR27878129_1.fastq /home/youngdeblynn/project/raw_data/SRR27878129_2.fastq -o /home/youngdeblynn/project/qc/

### Trimming
trimmomatic PE -phred33 \
  /home/youngdeblynn/project/raw_data/SRR27878129_1.fastq /home/youngdeblynn/project/raw_data/SRR27878129_2.fastq \
  /home/youngdeblynn/project/trimmed_data/SRR27878129_1_paired.fastq /home/youngdeblynn/project/trimmed_data/SRR27878129_1_unpaired.fastq \
  /home/youngdeblynn/project/trimmed_data/SRR27878129_2_paired.fastq /home/youngdeblynn/project/trimmed_data/SRR27878129_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### Alignment
bwa index /home/youngdeblynn/ncbi_dataset/data/GCF_020735445.1/GCF_020735445.1_ASM2073544v1_genomic.fna
bwa mem /home/youngdeblynn/ncbi_dataset/data/GCF_020735445.1/GCF_020735445.1_ASM2073544v1_genomic.fna \
  /home/youngdeblynn/project/trimmed_data/SRR27878129_1_paired.fastq \
  /home/youngdeblynn/project/trimmed_data/SRR27878129_2_paired.fastq \
  > /home/youngdeblynn/project/aligned_data/aligned_reads.sam

### Post-alignment Processing
samtools view -Sb /home/youngdeblynn/project/aligned_data/aligned_reads.sam | samtools sort -o /home/youngdeblynn/project/aligned_data/aligned_reads_sorted.bam
samtools index /home/youngdeblynn/project/aligned_data/aligned_reads_sorted.bam

picard MarkDuplicates \
  I=/home/youngdeblynn/project/aligned_data/aligned_reads_sorted.bam \
  O=/home/youngdeblynn/project/aligned_data/aligned_reads_marked.bam \
  M=/home/youngdeblynn/project/aligned_data/marked_dup_metrics.txt

### Variant Calling
gatk HaplotypeCaller \
  -R /home/youngdeblynn/ncbi_dataset/data/GCF_020735445.1/GCF_020735445.1_ASM2073544v1_genomic.fna \
  -I /home/youngdeblynn/project/aligned_data/aligned_reads_marked.bam \
  -O /home/youngdeblynn/project/variant_data/output_variants.vcf

## 6. Interpretation of Results

The final output is a VCF (Variant Call Format) file containing identified genetic variants. These variants represent differences between the sequenced Segatella copri genome and the reference genome. Analysis of these variants can provide insights into:

- Genetic diversity within the species
- Potential functional variations affecting bacterial behavior or interactions
- Evolutionary relationships between different strains

## 7. Future Applications

This genomic analysis of Segatella (Prevotella) copri opens up several avenues for future research and applications:

1. Microbiome Studies: Understanding genetic variations in S. copri can inform broader studies on gut microbiome composition and its impact on human health.

2. Strain Identification: The identified variants can be used to develop strain-specific markers for more accurate identification of S. copri in complex microbial communities.

3. Functional Genomics: Variants in specific genes can be further studied to understand their impact on bacterial metabolism, antibiotic resistance, or interactions with the host.

4. Comparative Genomics: This data can be used in comparative studies with other Prevotella species or strains to understand evolutionary relationships and adaptations.

5. Precision Medicine: Understanding genetic variations in gut bacteria like S. copri could contribute to personalized medicine approaches, particularly in gastrointestinal health.

6. Probiotic Development: Insights from this genomic analysis could inform the development of probiotic strains with specific beneficial characteristics.

7. Diagnostic Tools: The identified genetic markers could be used to develop diagnostic tools for detecting specific S. copri strains in clinical samples.

## 8. Conclusion

This project provides a comprehensive genomic analysis of Segatella (Prevotella) copri, from raw sequencing data to the identification of genetic variants. The workflow and results presented here serve as a foundation for further research into this important gut microbe, potentially contributing to our understanding of microbiome dynamics and human health.
