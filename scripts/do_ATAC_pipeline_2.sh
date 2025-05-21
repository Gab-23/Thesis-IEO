#!/bin/bash

NUM_THREADS=$1
FIXED_BAM=$2; ACC_BAM=$3; SHIFTED_BAM=$4
LOWER=$5; UPPER=$6

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo -e "\n"
  echo "
This script is to be run after do_ATAC_pipeline.sh
It uses the clean BAM provided as output and filters for a user specified IS
After doing that it computes default ATAC read shifting to account for Tn5 insertion
  "
  echo "HOW TO USE: "
  echo -e "\n"
  echo "./do_ATAC_pipeline_2.sh 
  [NUM_THREADS] 
  [FIXED_BAM] [ACC_BAM] [SHIFTED_BAM] 
  [LOWER] [UPPER] 				NOTA: these parameters determine the [LOWER] and [UPPER] boundaries of selected fragments
  "
  echo -e "\n"
  exit 0
fi

# define mamba environment
ENV_NAME=atacseq_preproc_2

 # no locks have to be set
export MAMBA_NO_LOCK=1

# initialize mamba
eval "$(mamba shell hook --shell bash)"

# activate mamba environment
mamba activate $ENV_NAME

if [ ! -f $SHIFTED_BAM ]; then
	echo "Filtering BAM file to retain NFR"
	# filter BAM to retain only the specified IS
	samtools view -@ $NUM_THREADS -h $FIXED_BAM | \
		awk -v lower=$LOWER -v upper=$UPPER 'substr($0,1,1)=="@" || ($9 >= lower && $9 <= upper) || ($9 <= -lower && $9 >= -upper)' | \
		samtools view -@ $NUM_THREADS -b > $ACC_BAM

else
	echo "Skipping BAM file filtering"
fi

if [ ! -f $SHIFTED_BAM ]; then
	echo "Indexing filtered BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $ACC_BAM
else
	echo "Skipping BAM file indexing"
fi

if [ ! -f $SHIFTED_BAM ]; then
	echo "Shifting reads "
	# perform BAM shifting
	alignmentSieve --ATACshift -p $NUM_THREADS -b $ACC_BAM -o $SHIFTED_BAM
else
	echo "Skipping shifting step"
fi

# compute the number of reads in the final BAM file
NUM_READS=$(samtools view -@ $NUM_THREADS -c $SHIFTED_BAM)

echo "ENCODE suggests more than 50M reads (25M fragments) for PE assays"
echo "This run yielded $(($NUM_READS / 1000000))M reads"

echo "Done"
exit 0

