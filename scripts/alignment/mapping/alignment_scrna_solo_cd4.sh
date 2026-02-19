#!/bin/bash
set -exuo pipefail
TOP_DIR=../data/scRNA-seq/
SCRIPT=../scripts/alignment/
GENOME=../data/genome/genome_p_cd4


cd $TOP_DIR
for dir in 305264 305265 305266 305267
do
break
cd fastq/$dir
ls  *.gz | xargs -I{} fastqc -o ${TOP_DIR}fastqc/${dir} {}
cd ../../
done


for dir in 305264 305265 305266 305267
do
cd ${TOP_DIR}/fastq/${dir}
ls
find ./ -iname "Ballouz*R2_001.fastq.gz" 2>&1 | grep -v Permission | xargs -I{} echo -n {} | sed 's/\.\//,/g' | sed 's/^,//g' | xargs -I{} echo {} | while read -r line; do
STAR --outSAMtype BAM Unsorted \
 --soloType Droplet \
 --soloCBwhitelist $SCRIPTS/3M-february-2018.txt \
 --soloCellFilter CellRanger2.2 8000 0.99 10 --soloFeatures Gene Velocyto \
 --quantMode GeneCounts \
 --readFilesCommand zcat \
 --soloBarcodeReadLength 28 \
 --runThreadN 5 \
 --genomeDir $GENOME \
 --readFilesIn ${line} ${line//R2/R1} --outFileNamePrefix ${TOP_DIR}/bam_p_cd4/${dir}
done
