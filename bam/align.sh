#!/bin/bash

#sample=DS01-P
#sample=DS02-M
#sample=DS03-C2
#sample=DS04-C19

#sample=DS01-P-2
sample=DS02-P-3

ref=../reference/brca1-variants.fasta
r1=../fastq/${sample}_R1_001.fastq.gz
r2=../fastq/${sample}_R2_001.fastq.gz

bwa mem $ref $r1 $r2 | samtools view -b > ${sample}.raw.bam
samtools sort ${sample}.raw.bam -o ${sample}.bam
rm ${sample}.raw.bam
samtools index ${sample}.bam
