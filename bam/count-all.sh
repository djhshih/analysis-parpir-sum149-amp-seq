#!/bin/bash

set -euo pipefail

for f in *.bam; do
	echo $f
	./count.sh $f 
done
