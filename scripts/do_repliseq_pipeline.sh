#!/bin/bash

OUTDIR_FASTQ=$1; READ_FWD=$2; READ_REV=$3
NUM_THREADS=$4
INDEX=$5
READ_VAL_FWD=$6; READ_VAL_REV=$7
INITIAL_BAM=$8; INITIAL_SORTED_BAM=$9
CANON_BAM=${10}; CANON_SORTED_BAM=${11}
BAMQC_OUTDIR=${12}
IS_PATTERNED=${13}; DEDUP_BAM=${14}; DEDUP_METRICS=${15}
FILTERED_BAM=${16}
BLACKLISTED_REGIONS=${17}
NO_BLACKLISTED_BAM=${18}
BINFILE=${19}

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo -e "\n"
	echo "placeholder text"
	echo "HOW TO USE: "
	echo -e "\n"
	echo "./do_hic_pipeline.sh 
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
  [NO_BLACKLISTED_BAM]
  [BINFILE]	
	"
  	echo -e "\n"
	exit 0
fi

source ~/miniconda3/etc/profile.d/conda.sh # initialize conda
conda activate repliseq_env # activate environment
echo "Activated ATAC-seq preprocessing environment"

if [ ! -f $cicciopasticcio ]; then
	echo "Trimming reads..."
	# trim reads using TrimGalore!
	trim_galore --fastqc --gzip --paired -o $OUTDIR_FASTQ $READ_FWD $READ_REV
else
	echo "Skipping read trimming step"
fi

if [ ! -f $cicciopasticcio ]; then
	echo "Aligning reads..."
	# align reads on reference genome using bwa-mem
	bwa mem -t $NUM_THREADS $INDEX $READS_VAL_FWD $READS_VAL_REV |  samtools view -S -h -b > $INITIAL_BAM
else
	echo "Skipping alignment step"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Sorting initial BAM file"
	# sort BAM
	samtools sort -@ $NUM_THREADS -o $INITIAL_SORTED_BAM $INITIAL_BAM
else
	echo "Skipping initial BAM file sorting"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Indexing initial BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $INITIAL_SORTED_BAM
else
	echo "Skipping initial BAM file indexing"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Filtering for canonical chromosomes..."
	# retain canonical chromosomes
	samtools view -@ $NUM_THREADS -o $CANON_BAM $INITIAL_SORTED_BAM $(echo chr{1..22} chrX chrY)
else
	echo "Skipping filtering for canonical chromosomes"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Sorting canon BAM file"
	# sort BAM
	samtools sort -@ $NUM_THREADS -o $CANON_SORTED_BAM $CANON_BAM
else
	echo "Skipping canon BAM file sorting"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Indexing canon BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $CANON_SORTED_BAM
else
	echo "Skipping canon BAM file indexing"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Performing BAMQC"
	# perform BAMQC
	qualimap bamqc --java-mem-size=20G -bam $CANON_SORTED_BAM -outdir $BAMQC_OUTDIR -outformat HTML -c 
else
	echo "Skipping BAMQC execution"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Spotting and removing duplicates"

	# spot and remove PCR and optical duplicates
	# 2500 works best for patterned sequencers, in my case they used Illumina NovaSeq 6000
	# which is patterned
	# NOTA: This works only with original sequencer headers. If you download data from SRA the header will be different
	# 	and thus the optical duplicates won't be spotted. This means that to download reads one should use the deprecated
	#	fastq-dump with --origfmt option, which is not supported in fasterq-dump (wtf)
	
	if [ $IS_PATTERNED == true ] ; then
		PIXEL_DISTANCE=2500   		 
	else
		PIXEL_DISTANCE=100
	fi

	echo "Is the flowcell patterned? $IS_PATTERNED"
        echo "I will use $PIXEL_DISTANCE as pixel distance"

	picard MarkDuplicates I=$CANON_SORTED_BAM O=$DEDUP_BAM M=$DEDUP_METRICS \
		REMOVE_DUPLICATES=true \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=$PIXEL_DISTANCE\
		ASSUME_SORTED=true
else
	echo "Skipping duplicates identification"
fi

if [ ! -f $NO_BLACKLISTED_BAM ]; then
	echo "Removing blacklisted regions from BAM file"
	# remove blacklisted regions
	bedtools intersect -v -abam $DEDUP_BAM -b $BLACKLISTED_REGIONS > $NO_BLACKLISTED_BAM
else
	echo "Skipping blacklisted regions removal"

if [ ! -f $OUTFILE_COVERAGE ]; then
	# compute coverage
	bedtools coverage -counts -sorted -a $BINFILE -b $NO_BLACKLISTED_BAM > $OUTFILE_COVERAGE
else
	echo "Skipping coverage computation"
fi

echo "Done"
exit 0



