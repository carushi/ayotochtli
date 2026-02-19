#!/bin/bash
set -exuo pipefail 
TOP_DIR=../data/scATAC-seq/
OUT_DIR=$TOP_DIR
# personal
# GENOME=../data/genome/quad16-90/
# public
# GENOME=../data/genome/genome/
# public refseq
GENOME=../data/genome/genome_refseq/


SCRIPT=../scripts/alignment/
STAR=../reference/STAR-2.7.7a/bin/Linux_x86_64/STAR

for DIR in 307369 307370 307371 307372
do
find ${TOP_DIR}/$DIR/*R1_*.fastq.gz | while read file
do
bfname="${file%R1_001.fastq.gz}"
echo $bfname
obfname="${OUT_DIR}/$DIR/$(basename -- $bfname)"
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
${STAR} --outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI CR CB \
--soloType CB_samTagOut \
--soloCBmatchWLtype 1MM \
--soloBarcodeReadLength 16 \
--twopassMode Basic \
--soloFeatures Gene GeneFull SJ Velocyto \
--soloCBwhitelist ${SCRIPT}/whitelist/737K-cratac-v1.txt \
--sjdbFileChrStartEnd $GENOME/sjdbList.fromGTF.out.tab \
--readFilesCommand zcat \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 15 \
--outTmpDir ./startmp --runThreadN 5 --genomeDir $GENOME \
--readFilesIn ${bfname}R1_001.fastq.gz ${bfname}R3_001.fastq.gz ${bfname}R2_001.fastq.gz --outFileNamePrefix ${obfname}_relaxed
bzip2 ${bfname}_trimmed_R*.fastq
done
done
find $TOP_DIR/ -name "*.bam" | xargs -I{} echo "samtools index {}" | bash
find $TOP_DIR/ -name "*.sam" -exec rm {} \;


