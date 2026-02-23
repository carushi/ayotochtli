#!/bin/bash
set +eu
if [[ -z "$1" ]]; then
	TOP_DIR=../data/bam_relaxed/
else
	TOP_DIR=$1
fi
echo $TOP_DIR
SCRIPT=../scripts/alignment
GENOME=../data/genome/quad16-90
GENOME_SIZE=../data/genome/genome.chrom.sizes
STAR=../reference/STAR-2.7.7a/bin/Linux_x86_64/STAR

activate_genome () {
	set +eu
	source activate genome
	set -exuo pipefail
}

activate_snap () {
	set +eu
	source activate snapatac
	set -exuo pipefail
}

for file in ` ls $TOP_DIR/*.out.bam `
do
	bfname="${file%__relaxedAligned.sortedByCoord.out.bam}"
	echo ${bfname}.bam
	samtools view ${file} -H > ${bfname}.sam
	# create a bam file with the barcode embedded into the read name
	cat <( cat ${bfname}.sam ) \
		<( samtools view ${file} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CR:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CR"], $0 }' ) \
		| samtools view -bS - | samtools sort -n  - > ${bfname}.bam
	snaptools snap-pre \
		--input-file=${bfname}.bam \
		--output-snap=${bfname}.snap \
		--genome-size=${GENOME_SIZE}/data/rkawaguc/data_armadillo/genome_data/quad16-90/quad16-90.chrom.sizes \
		--genome-name=dasNov3 \
		--min-map=30 \
		--min-flen=0 \
		--max-flen=1000 \
		--keep-chrm=TRUE \
		--keep-single=TRUE \
		--keep-secondary=False \
		--overwrite=True \
		--max-num=1000000 \
		--min-cov=100 \
		--verbose=True
 	snaptools snap-add-bmat  \
		--snap-file=${bfname}.snap  \
		--bin-size-list 1000 5000  \
		--verbose=True		
done

