#!/bin/bash


SUMMITS_BED=$1
INPUT_GENOME=$2
HOW_MUCH=$3
OUTFILE_NAME=$4


if [ $1 = "-h" ] || [ $1 = "--help" ]; then
        echo -e "\n"
        echo "
This function takes as input summit coordinates and expands them to 10kbp using bedtools slop
        "
        echo "HOW TO USE: "
        echo -e "\n"
        echo "./expand_summits.sh [SUMMITS_BED] [INPUT_GENOME] [HOW_MUCH] [OUTFILE_NAME]"
        echo -e "\n"
        exit 0
fi

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate environment
conda activate peak_calling
echo "Activated peak calling environment"

bedtools slop -i $SUMMITS_BED -g $INPUT_GENOME -b $HOW_MUCH > $OUTFILE_NAME









