#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

ARGUMENTS="$@"
NUM_OF_ARGUMENTS="$#"

if [ $NUM_OF_ARGUMENTS -gt 1 ]; then
	echo "Detected multiple SRA accessions"
	for dir in "$ARGUMENTS"
	do
    		echo "Retrieving reads for $dir"
    		fasterq-dump --split-files "$dir/$dir.sra" -O $dir
		find $dir -type f \( -name "*.fastq" -o -name "*.sra" \) -exec pigz --fast {} +
	done
	echo "Done!"
	exit 0
elif [ $NUM_OF_ARGUMENTS -eq 1 ]; then
	echo "SRA accession file detected"
	
	while IFS= read -r accession; do
		echo "Retrieving reads for $accession"
		fasterq-dump --split-files "$accession/$accession.sra" -O "$accession";
		find "$accession" -type f \( -name "*.fastq" -o -name "*.sra" \) -exec pigz --fast {} +
	done < $ARGUMENTS

	echo "Done!"
	exit 0
else
	echo "No argument detected"
	exit 1
fi
