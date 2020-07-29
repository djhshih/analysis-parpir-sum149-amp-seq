#!/bin/bash

bam=$1

# only look at read 1, since read 2 often cannot be mapped uniquely
# due to the nature of the amplicon
samtools view -f 0x0040 $1 -q 5 | cut -f 3 | uniq -c
