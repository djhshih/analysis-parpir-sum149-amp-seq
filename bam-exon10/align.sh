#!/bin/bash

#sample=DS01-M
sample=DS02-P
#sample=DS03-C2
#sample=DS04-C19

ref=../reference/brca1-enst00000471181-exon10.fasta
r1=../fastq/${sample}_R1_001.fastq.gz
r2=../fastq/${sample}_R2_001.fastq.gz

bwa mem $ref $r1 $r2 | samtools view -b > ${sample}.raw.bam
samtools sort ${sample}.raw.bam -o ${sample}.bam
rm ${sample}.raw.bam
samtools index ${sample}.bam
