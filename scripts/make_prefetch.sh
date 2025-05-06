#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

ARGUMENTS="$@"
NUM_OF_ARGUMENTS="$#"

if [ $NUM_OF_ARGUMENTS -gt 1 ]; then
	echo "Detected multiple SRA accessions"
	for ARGUMENT in "$ARGUMENTS"
	do
    		echo "Prefetching $ARGUMENT"
    		prefetch $ARGUMENT
	done
	exit 0
elif [ $NUM_OF_ARGUMENTS -eq 1 ]; then
	echo "SRA accession file detected"
	cat $ARGUMENTS | xargs -n 1 prefetch
	exit 0
else
	echo "No argument detected"
	exit 1
fi


