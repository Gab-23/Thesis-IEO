#!/bin/bash

METADATA_TAB=$1
TYPE=$2; STUDY=$3
NUM_THREADS=$4
MIN_LENGTH=$5; QVALUE=$6
PQ_VALS_OUTFILE=$7; PILEUPS_OUTFILE=$8; 
ANALYZED_READS_OUTFILE=$9
NARROWPEAK_OUTFILE=${10}

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo -e "\n"
	echo "
This function uses Genrich framework to call peaks given a set of .bam files.
The function expects a .tsv table as input containing the sampleIDs and the condition to filter them by
Genrich automatically aggregates peaks over samples using Fischer's method, instead of IDR
	"
	echo "HOW TO USE: "
	echo -e "\n"
	echo "./do_peak_calling.sh
	[METADATA_TAB] 				NOTA: this assumes that 1st column contains sampleID and the 6th column contains the sample type	
	[TYPE] [STUDY]				NOTA: type is assumed to be sample type of interest
        [NUM_THREADS]	
	[MIN_LENGTH] [QVALUE]			NOTA: minimum peak length / FDR threshold 
	[PQ_VALS_OUTFILE] [PILEUPS_OUTFILE] 
	[ANALYZED_READS_OUTFILE]
	[NARROWPEAK_OUTFILE]
	"
  	echo -e "\n"
	exit 0
fi

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate environment
conda activate peak_calling
echo "Activated peak calling environment" 

# load *_shifted.bam array
readarray -t BAM_ARRAY< <(sed -e '1d' $METADATA_TAB | \
	awk -v type=$TYPE -v study=$STUDY '$6 == type && $9 == study' | \
	cut -f 1 | \
	xargs -I{} echo "$STUDY/bam/{}/{}_shifted.bam")

# load *_name_sorted.bam array
readarray -t BAM_ARRAY_NAME_SORTED< <(sed -e '1d' $METADATA_TAB | \
	awk -v type=$TYPE -v study=$STUDY '$6 == type && $9 == study' | \
	cut -f 1 | \
	xargs -I{} echo "$STUDY/bam/{}/{}_name_sorted.bam")

# initialize Genrich comma separated path variable
GENRICH_INPUT=$(printf "%s\n" "${BAM_ARRAY_NAME_SORTED[@]}" | paste -sd ,)

# iterate over indexes
for idx in ${!BAM_ARRAY[@]}; do
	# take current input_bam
	input_bam=${BAM_ARRAY[$idx]}
	# take current output_bam
	output_bam=${BAM_ARRAY_NAME_SORTED[$idx]}
	
	# if output_bam does not exist
	if [ ! -f $output_bam ]; then
		echo "Sorting BAM file according to read name"
		# sort according to read name
		samtools sort -@ $NUM_THREADS -n $input_bam -o $output_bam
	else
		echo "Skipping BAM file sorting [QNAME]"
	fi

done

# if output narrowPeak does not exist
if [ ! -f $NARROWPEAK_OUTFILE ]; then
	echo "Calling peaks (hope they answer)"
	
	# -j is for ATAC setting
	# -D skips Tn5 correction

	Genrich -t $GENRICH_INPUT -j -D \
		-l $MIN_LENGTH \
		-q $QVALUE \
		-f $PQ_VALS_OUTFILE \
		-k $PILEUPS_OUTFILE \
		-b $ANALYZED_READS_OUTFILE \
		-o $NARROWPEAK_OUTFILE
else
	echo "Skipping peak calling step"
fi

echo "Done!"
exit 0
