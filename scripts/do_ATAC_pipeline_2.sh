#!/bin/bash

NUM_THREADS=$1
NOBLACK_BAM=$2; ACC_BAM=$3; SHIFTED_BAM=$4
LOWER=$5; UPPER=$6

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo -e "\n"
  echo "HOW TO USE: "
  echo "./do_ATAC_pipeline_2.sh 
  [NUM_THREADS] 
  [NOBLACK_BAM] [ACC_BAM] [SHIFTED_BAM] 
  [LOWER] [UPPER]
  "
  echo -e "\n"
  exit 0
fi

ENV_NAME=atacseq_preproc_2

if [ ! -f $SHIFTED_BAM ]; then
	echo "Filtering BAM file to retain NFR"
	mamba run -n $ENV_NAME samtools view -@ $NUM_THREADS -h $NOBLACK_BAM | \
		awk -v lower=$LOWER -v upper=$UPPER 'substr($0,1,1)=="@" || ($9 >= lower && $9 <= upper) || ($9 <= -lower && $9 >= -upper)' | \
		mamba run -n $ENV_NAME samtools view -@ $NUM_THREADS -b > $ACC_BAM
else
	echo "Skipping BAM file filtering"
fi

if [ ! -f $SHIFTED_BAM ]; then
	echo "Indexing filtered BAM file"
	mamba run -n $ENV_NAME samtools index -@ $NUM_THREADS $ACC_BAM
else
	echo "Skipping BAM file indexing"
fi

if [ ! -f $SHIFTED_BAM ]; then
	echo "Shifting reads "
	mamba run -n $ENV_NAME alignmentSieve --ATACshift -p $NUM_THREADS -b $ACC_BAM -o $SHIFTED_BAM
else
	echo "Skipping shifting step"
fi

NUM_READS=$(samtools view -@ $NUM_THREADS -c $SHIFTED_BAM)

echo "ENCODE suggests more than 50M reads (25M fragments) for PE assays"
echo "This run yielded $(($NUM_READS / 1000000))M reads"

echo "Done"
exit 0




