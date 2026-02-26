#!/bin/bash
set -exuo pipefail
TOP_DIR=../data/bam/scATAC/
OUT_DIR=$TOP_DIR

STAR=../reference/STAR-2.7.9a/bin/Linux_x86_64/STAR
SCRIPT=../scripts/alignment

# personal
# GENOME=../data/genome/quad16-90/
# OUTDIR=../data/bam_pub_genome_cd4/
# public
GENOME=../data/genome/genome/
OUTDIR=../data/scATAC/bam_pub_genome/
# public refseq
# GENOME=../data/genome/genome_refseq/
# OUTDIR=../data/bam_pub_refseq/

if false ; then
    READFILE=""
    EXT=""
else
    READFILE=" --readFilesCommand zcat "
    EXT=".gz"
fi

echo ${EXT}
# for DIR in $1
for DIR in 307369 307370 307371 307372
do
find ${TOP_DIR}/$DIR/*_R1_*.fastq${EXT} | while read file
do
bfname="${file%R1_001.fastq${EXT}}"
echo $bfname
mkdir ${OUTDIR}/$DIR
obfname="${OUTDIR}/$DIR/$(basename -- $bfname)"
echo $obfname
echo $bfname
echo $file
if test -f "${obfname}_Log.final.out"; then
	continue
fi

echo ${bfname}R3_001.fastq.gz
fastqc ${bfname}R1_001.fastq.gz
fastqc ${bfname}R3_001.fastq.gz
cutadapt --cores=6 --minimum-length 5 --nextseq-trim=20 -a AGATCGGAAGAG -A AGATCGGAAGAG -o ${bfname}_trimmed_R1.fastq -p ${bfname}_trimmed_R2.fastq ${bfname}R1_001.fastq.gz ${bfname}R2_001.fastq.gz
bunzip2 ${bfname}_trimmed_R*fastq.bz2
gunzip ${bfname}R1_001.fastq.gz
gunzip ${bfname}R3_001.fastq.gz
gunzip ${bfname}R2_001.fastq.gz
if test -d "../startmp_${DIR}" ; then
    rm -r  ../startmp_${DIR}
fi


${STAR} --outSAMtype BAM Unsorted \
--outSAMattributes NH HI CR CB \
--soloType CB_samTagOut \
--soloCBmatchWLtype 1MM \
--soloBarcodeReadLength 16 \
--twopassMode Basic \
--soloFeatures Gene GeneFull Velocyto \
--soloCBwhitelist ${SCRIPT}/whitelist/737K-cratac-v1.txt \
--runThreadN 5 --genomeDir $GENOME $READFILE \
--readFilesIn ${bfname}R1_001.fastq${EXT} ${bfname}R3_001.fastq${EXT} ${bfname}R2_001.fastq${EXT} --outFileNamePrefix ${obfname} \
--outTmpDir ../startmp_${DIR} \
--alignIntronMax 1 --alignMatesGapMax 1000
bzip2 ${bfname}_trimmed_R*.fastq
done
done

find $TOP_DIR/ -name "*.bam" | xargs -I{} echo "samtools index {}" | bash
find $TOP_DIR/ -name "*.sam" -exec rm {} \;

