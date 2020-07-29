#!/bin/bash

#sample=DS01-M
#sample=DS02-P
sample=DS03-C2
#sample=DS04-C19

ref=../reference/brca1-enst00000471181-exon10.fasta
#ref=../reference/brca1-enst00000471181-amplicon1.fasta
r1=../fastq/${sample}_R1_001.fastq.gz
r2=../fastq/${sample}_R2_001.fastq.gz

# critical parameter
library_length=357

razers3 -o ${sample}.raw.bam -rr 100 -i 90 -mr 10 --unique -ll $library_length $ref $r1 $r2 -vv
samtools sort ${sample}.raw.bam -o ${sample}.bam
rm ${sample}.raw.bam
samtools index ${sample}.bam
