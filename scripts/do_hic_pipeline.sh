#!/bin/bash

HIC_ENVIRONMENT=$1; HIC_QC_ENVIRONMENT=$2; HIC_PROCESSING_ENVIRONMENT=$3
HIC_QC_SCRIPT=$4
INDEX=$5; ASSEMBLY_CODE=$6
READS_FWD=$7; READS_REV=$8
NUM_THREADS=$9
OUTFILE_BAM=${10}; OUTFILE_SORTED_COORDS_BAM=${11}; OUTFILE_FILTERED_BAM=${12}; OUTFILE_SORTED_QNAME_BAM=${13}
OUTFILE_QC=${14}
OUTFILE_PAIRS=${15}; OUTFILE_SORTED_PAIRS=${16}; OUTFILE_STATS=${17}
CHROM_SIZES_CANON=${18}
OUTFILE_NODUPS_PAIRS=${19}; OUTFILE_DEDUP_STATS=${20}
OUTFILE_NODUPS_UU_PAIRS=${21}; OUTFILE_CLEAN_STATS=${22}
OUTDIR_MULTIQC=${23}
BACKBONE_CLEAN_NO_HEADER=${24}; OUTFILE_COOL_MATRIX=${25}

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo -e "\n"
	echo "HOW TO USE: "
	echo "./do_hic_pipeline.sh 
  [HIC_ENVIRONMENT] [HIC_QC_ENVIRONMENT] [HIC_PROCESSING_ENVIRONMENT] 
  [HIC_QC_SCRIPT]
  [INDEX] [ASSEMBLY_CODE]
  [READS_FWD] [READS_REV]
  [NUM_THREADS] 
  [OUTFILE_BAM] [OUTFILE_SORTED_COORDS_BAM] [OUTFILE_FILTERED_BAM] [OUTFILE_SORTED_QNAME_BAM]
  [OUTFILE_QC]
  [OUTFILE_PAIRS] [OUTFILE_SORTED_PAIRS] [OUTFILE_STATS]
  [CHROM_SIZES_CANON]
  [OUTFILE_NODUPS_PAIRS] [OUTFILE_DEDUP_STATS] 
  [OUTFILE_NODUPS_UU_PAIRS] [OUTFILE_CLEAN_STATS]
  [OUTDIR_MULTIQC]
  [BACKBONE_CLEAN_NO_HEADER] [OUTFILE_COOL_MATRIX]"
  	echo -e "\n"
	exit 0
fi

source ~/miniconda3/etc/profile.d/conda.sh

echo "Activating $HIC_ENVIRONMENT environment"
conda activate $HIC_ENVIRONMENT

echo "Executing Pipeline according to PhaseGenomics"
echo "For more details, look at: "
echo "https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html"

if [ ! -f $OUTFILE_SORTED_QNAME_BAM ]; then
	echo "Aligning reads to the genome using BWA-MEM"
	bwa mem -t $NUM_THREADS -5SP $INDEX $READS_FWD $READS_REV | samblaster | samtools view -S -h -b -F 2316 > $OUTFILE_BAM
else
	echo "Skipping alignment step"
fi

if [ ! -f $OUTFILE_SORTED_QNAME_BAM ]; then
	echo "Sorting BAM file according to coordinates"
	samtools sort -@ $NUM_THREADS $OUTFILE_BAM -o $OUTFILE_SORTED_COORDS_BAM
else
	echo "Skipping BAM file sorting [COORDS]"
fi

if [ ! -f $OUTFILE_SORTED_QNAME_BAM ]; then
	echo "Index BAM file"
	samtools index -@ $NUM_THREADS $OUTFILE_SORTED_COORDS_BAM
else
	echo "Skipping BAM file indexing"
fi

if [ ! -f $OUTFILE_SORTED_QNAME_BAM ]; then
	echo "Filter only canonical chromosomes" 
	samtools view -@ $NUM_THREADS -o $OUTFILE_FILTERED_BAM $OUTFILE_SORTED_COORDS_BAM $(echo chr{1..22} chrX chrY)
else
	echo "Skipping BAM file chromosome filtering"
fi

if [ ! -f $OUTFILE_SORTED_QNAME_BAM ]; then
	echo "Sorting BAM file according to read name"
	samtools sort -@ $NUM_THREADS -n $OUTFILE_FILTERED_BAM -o $OUTFILE_SORTED_QNAME_BAM
	
	echo "Removing all intermediate BAM files"
	rm $OUTFILE_BAM $OUTFILE_SORTED_COORDS_BAM "$OUTFILE_SORTED_COORDS_BAM.bai" $OUTFILE_FILTERED_BAM
else
	echo "Skipping BAM file sorting [QNAME]"
fi

echo "Activating $HIC_QC_ENVIRONMENT environment"
conda activate $HIC_QC_ENVIRONMENT

if ls ${OUTFILE_QC}* 1> /dev/null 2>&1; then
	echo "Skipping HiC-QC execution"	
else
	echo "Running HiC-QC script"
	python $HIC_QC_SCRIPT -b $OUTFILE_SORTED_QNAME_BAM -o $OUTFILE_QC
fi

echo "Activating $HIC_PROCESSING_ENVIRONMENT environment"
conda activate $HIC_PROCESSING_ENVIRONMENT

if [ ! -f $OUTFILE_PAIRS ]; then
	echo "Running pairtools parse"
	pairtools parse -o $OUTFILE_PAIRS -c $CHROM_SIZES_CANON \
  	--drop-sam --drop-seq --output-stats $OUTFILE_STATS \
  	--assembly $ASSEMBLY_CODE \
  	--flip \
  	--add-columns mapq \
  	--walks-policy mask \
  	--nproc-in $NUM_THREADS \
  	--nproc-out $NUM_THREADS \
  	$OUTFILE_SORTED_QNAME_BAM
else
	echo "Skipping pairtools parse execution"
fi

if [ ! -f $OUTFILE_SORTED_PAIRS ]; then
	echo "Running pairtools sort"
	pairtools sort --nproc $NUM_THREADS \
		-o $OUTFILE_SORTED_PAIRS \
		$OUTFILE_PAIRS
else
	echo "Skipping pairtools sort execution"
fi	

if [ ! -f $OUTFILE_NODUPS_PAIRS ]; then
    echo "Running pairtools dedup"

     pairtools dedup \
         --nproc-in    $NUM_THREADS \
         --nproc-out   $NUM_THREADS \
         --max-mismatch 3 \
         --mark-dups \
         --output       $OUTFILE_NODUPS_PAIRS \
         --output-stats $OUTFILE_DEDUP_STATS \
         $OUTFILE_SORTED_PAIRS &

    DEDUP_PID=$!
    while kill -0 $DEDUP_PID ; do
    	echo "Process is still active..."
   	sleep 1
    done

else
    echo "Skipping pairtools dedup and execution"
fi

if [ ! -f $OUTFILE_NODUPS_UU_PAIRS ]; then
	echo "Running pairtools select"
	pairtools select "mapq1>=30 and mapq2>=30" \
		--nproc-in $NUM_THREADS \
		--nproc-out $NUM_THREADS \
		$OUTFILE_NODUPS_PAIRS \
		-o $OUTFILE_NODUPS_UU_PAIRS
else
	echo "Skipping pairtools select execution"
fi

if [ ! -f $OUTFILE_CLEAN_STATS ]; then
	echo "Running pairtools stats"
	pairtools stats $OUTFILE_NODUPS_UU_PAIRS -o $OUTFILE_CLEAN_STATS
else
	echo "Skipping pairtools stats execution"
fi

if ls ${OUTDIR_MULTIQC}* 1> /dev/null 2>&1; then
	echo "Skipping multiQC excution"
else
	echo "Running multiQC"
	multiqc --outdir $OUTDIR_MULTIQC $OUTFILE_CLEAN_STATS
fi

if [ ! -f $OUTFILE_COOL_MATRIX ]; then
	echo "Cooling the data..."
	cooler cload pairs -c1 2 -p1 3 \
		--no-symmetric-upper \
		-c2 4 -p2 5 \
       		--assembly $ASSEMBLY_CODE \
      		$BACKBONE_CLEAN_NO_HEADER \
      		$OUTFILE_NODUPS_UU_PAIRS \
      		$OUTFILE_COOL_MATRIX
else
	echo "Skipping .cool matrix generation"
fi

echo "What's Next?"
echo "[OPTIONAL] Merge the .cool matrices"
echo "[MANDATORY] Apply ICE normalization"
echo "[MANDATORY] Apply fixed-distance Z-score normalization"
echo "Done!"
exit 0
