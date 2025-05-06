#!/bin/bash

CANCER_TYPE="$1"
ORIGINAL_PEAKS="$2"
COMPLEMENTARY="$3"

source ~/miniconda3/etc/profile.d/conda.sh

if [ $1 = "-h" ] | [ $1 = "--help" ]; then
        echo -e "\n"
        echo "HOW TO USE: "
        echo "./filter_fragments.sh 
  [CANCER_TYPE]
  [ORIGINAL_PEAKS]
  [COMPLEMENTARY]"
        echo -e "\n"
        exit 0
fi

conda activate fragments_filtering_env

fragments_array=$(find /home/ieo7429/Desktop/THESIS_GAB/samples/$CANCER_TYPE/HT* -type f -name "*fragments_lifted.tsv.gz")

for fragment_file in $fragments_array;
do
	sample_name=$(echo $fragment_file | cut -d '/' -f 8)
	echo "Processing $sample_name"
	if [ $COMPLEMENTARY == "-v" ]
	then
		echo "Looking for fragments NOT overlapping peaks" 
		bedtools intersect -a $fragment_file -b $ORIGINAL_PEAKS -wa -v | gzip > /home/ieo7429/Desktop/THESIS_GAB/samples/$CANCER_TYPE/$sample_name/"$sample_name"_fragments_lifted_complementary.tsv.gz
	else
		echo "Looking for fragments overlapping peaks uniquely"
		bedtools intersect -a $fragment_file -b $ORIGINAL_PEAKS -wa -u | gzip > /home/ieo7429/Desktop/THESIS_GAB/samples/$CANCER_TYPE/$sample_name/"$sample_name"_fragments_lifted_filtered.tsv.gz
	fi
done
echo "Done!"
exit 0
