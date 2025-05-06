#!/bin/bash

HIC_PROCESSING_ENVIRONMENT=$1
MERGED_OUTFILE=$2; #  MERGED_OUTFILE_SYMM=$3
NUM_THREADS=$3

shift 3

COOLER_PATHS=${@}

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
        echo -e "\n"
        echo "HOW TO USE: "
        echo "./do_hic_pipeline.sh 
 	[HIC_PROCESSING_ENVIRONMENT] 
	[MERGED_OUTFILE] [MERGED_OUTFILE_SYMM]
       	[NUM_THREADS] 
	[COOLER_PATHS]
       	"
        echo -e "\n"
        exit 0
fi

source ~/miniconda3/etc/profile.d/conda.sh

echo "Activating $HIC_PROCESSING_ENVIRONMENT environment"
conda activate $HIC_PROCESSING_ENVIRONMENT

echo "Merging coolers..."
cooler merge $MERGED_OUTFILE $COOLER_PATHS

echo "Balancing merged cooler..."
cooler balance -p $NUM_THREADS $MERGED_OUTFILE

# echo "Symmetrizing matrix..."
# cooler dump --balanced --join --fill-lower -o $MERGED_OUTFILE_SYMM $MERGED_OUTFILE

echo "Done!"
exit 0




