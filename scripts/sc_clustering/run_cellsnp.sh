#!/usr/bin/bash

for sample in 1 2 3 4
do
    for lane in 1 2 3 4
    do
        BARCODE=barcode_atac_S${sample}.tsv
        FILE=../data/bam/1690${sample}_S${sample}_L00${lane}_sorted.bam
        cellsnp-lite -s $FILE -b $BARCODE -O ./${sample}_${lane} -p 10 --minCOUNT 10 --gzip --genotype --chrom chrM --cellTAG CR:Z --UMItag None
        if [ $lane -eq 1 ];
        then
            BAM=$FILE
            echo
        else 
            BAM=${BAM},${FILE}
        fi
    done
done
# GT=genotype
# DP=combined depth across samples
# OTH=
# AD=Allelic depth
# PL=PHred scaled