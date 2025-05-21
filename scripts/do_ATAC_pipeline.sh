#!/bin/bash

OUTDIR_FASTQ=$1; READ_FWD=$2; READ_REV=$3
NUM_THREADS=$4
INDEX=$5
READ_VAL_FWD=$6; READ_VAL_REV=$7
INITIAL_BAM=$8; INITIAL_SORTED_BAM=$9
CANON_BAM=${10}; CANON_SORTED_BAM=${11}
BAMQC_OUTDIR=${12}
DEDUP_BAM=${13}; DEDUP_METRICS=${14}
FILTERED_BAM=${15}
BLACKLISTED_REGIONS=${16}
NO_BLACKLISTED_BAM=${17}; FIXED_BAM=${18}
IS_METRICS=${19}; IS_HIST=${20}

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo -e "\n"
	echo "
This function takes several arguments and performs standard ATAC-seq pipeline processing.
First it auto-detects adapter sequences and trims reads, then it computes QC.
Reads are then aligned to the reference genome and BAM files are filtered in order to obtain clean results
The script removes duplicates, discordant reads and also filters for MAPQ > 20. 
Finally blacklisted regions are removed.

The output is a clean BAM file, ready for IS filtering and Tn5 correction.
"
	echo "HOW TO USE: "
	echo -e "\n"
	echo "./do_ATAC_pipeline.sh 
  [OUTDIR_FASTQ] [READ_FWD] [READ_REV]
  [NUM_THREADS]
  [INDEX] 
  [READ_VAL_FWD] [READ_VAL_REV]
  [INITIAL_BAM] [INITIAL_SORTED_BAM] 
  [CANON_BAM] [CANON_SORTED_BAM]
  [BAMQC_OUTDIR]
  [DEDUP_BAM] [DEDUP_METRICS]
  [FILTERED_BAM]
  [BLACKLISTED_REGIONS]
  [NO_BLACKLISTED_BAM] [FIXED_BAM]
  [IS_METRICS] [IS_HIST]"
  	echo -e "\n"
	exit 0
fi

source ~/miniconda3/etc/profile.d/conda.sh # initialize conda
conda activate atacseq_preproc # activate environment
echo "Activated ATAC-seq preprocessing environment"

# check for presence of output file
if [ ! -f $FIXED_BAM ]; then
	echo "Trimming reads..."
	# trim reads using Trim Galore!
	trim_galore --fastqc --gzip --paired -o $OUTDIR_FASTQ $READ_FWD $READ_REV
else
	echo "Skipping read trimming step"
fi


if [ ! -f $FIXED_BAM ]; then
	echo "Aligning reads..."
	 # align reads on reference genome using bowtie2
	bowtie2 -q -p $NUM_THREADS --local --very-sensitive-local --no-discordant --no-mixed --dovetail \
		-x $INDEX -1 $READ_VAL_FWD -2 $READ_VAL_REV | samtools view -@ $NUM_THREADS -b -o $INITIAL_BAM

# NOTA: -q reads for .fastq
#	--local performs local alignment
#	--no-discordant avoids rescuing discordant pairs
# 	--no-mixed avoids having unpaired reads
#	--dovetail allows for reads to overlap

else
	echo "Skipping alignment step"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Sorting initial BAM file"
	# sort BAM
	samtools sort -@ $NUM_THREADS -o $INITIAL_SORTED_BAM $INITIAL_BAM
else
	echo "Skipping initial BAM file sorting"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Indexing initial BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $INITIAL_SORTED_BAM
else
	echo "Skipping initial BAM file indexing"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Filtering for canonical chromosomes..."
	# retain canonical chromosomes
	samtools view -@ $NUM_THREADS -o $CANON_BAM $INITIAL_SORTED_BAM $(echo chr{1..22} chrX chrY)
else
	echo "Skipping filtering for canonical chromosomes"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Sorting canon BAM file"
	# sort BAM
	samtools sort -@ $NUM_THREADS -o $CANON_SORTED_BAM $CANON_BAM
else
	echo "Skipping canon BAM file sorting"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Indexing canon BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $CANON_SORTED_BAM
else
	echo "Skipping canon BAM file indexing"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Performing BAMQC"
	# perform BAMQC
	qualimap bamqc --java-mem-size=20G -bam $CANON_SORTED_BAM -outdir $BAMQC_OUTDIR -outformat HTML -c 
else
	echo "Skipping BAMQC execution"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Spotting and removing duplicates"

	# spot and remove PCR and optical duplicates
	# 2500 works best for patterned sequencers, in my case they used Illumina NovaSeq 6000
	# which is patterned

	picard MarkDuplicates I=$CANON_SORTED_BAM O=$DEDUP_BAM M=$DEDUP_METRICS \
		REMOVE_DUPLICATES=true \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		ASSUME_SORTED=true
else
	echo "Skipping duplicates identification"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Filtering deduplicated BAM file"
	# remove reads based on bitwise flags
	samtools view -@ $NUM_THREADS -q 20 -f 0x2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800 -b $DEDUP_BAM > $FILTERED_BAM
else
	echo "Skipping filtering of deduplicated BAM file"
fi

if [ ! -f $FIXED_BAM ]; then
	echo "Removing blacklisted regions from BAM file"
	# remove blacklisted regions
	bedtools intersect -v -abam $FILTERED_BAM -b $BLACKLISTED_REGIONS > $NO_BLACKLISTED_BAM
else
	echo "Skipping blacklisted regions removal"
fi

if [ ! -f $FIXED_BAM ]; then
        echo "Fixing mate information from BAM file"
        # fixed read orientation
        picard FixMateInformation \
	I=$NO_BLACKLISTED_BAM \
	O=$FIXED_BAM \
	ADD_MATE_CIGAR=true

else
        echo "Skipping BAM fixing"
fi

if [ ! -f $IS_METRICS ] || [ ! -f $IS_HIST ]; then
	echo "Computation of Insert Size metrics"
	# obtain IS metrics
	picard CollectInsertSizeMetrics I=$FIXED_BAM O=$IS_METRICS H=$IS_HIST
else
	echo "Skipping Insert Size metrics computation"
	echo "Done!"
	exit 0
fi

echo "Removing intermediate files"
rm $INITIAL_BAM $INITIAL_SORTED_BAM "$INITIAL_SORTED_BAM.bai"
rm $CANON_BAM $CANON_SORTED_BAM "$CANON_SORTED_BAM.bai"
rm $DEDUP_BAM $FILTERED_BAM $NO_BLACKLISTED_BAM

echo "Done!"
exit 0



