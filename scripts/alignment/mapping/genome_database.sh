#!/bin/bash

GENOME=dasNov3_spike.fa
GENE=Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf

STAR --genomeSAindexNbases 14 --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $GENOME --sjdbGTFfile $GENE --limitGenomeGenerateRAM 40115748224

GENOME=16-90.trans.fa
GENE=16-90.merged.gtf

STAR --genomeSAindexNbases 14 --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $GENOME --sjdbGTFfile $GENE --limitGenomeGenerateRAM 40115748224


GENOME=dasNov3_spike.fa
GENE=Dasypus_novemcinctus.Dasnov3.0.95_mod.spike_cd4.gtf

STAR --genomeSAindexNbases 14 --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $GENOME --sjdbGTFfile $GENE --limitGenomeGenerateRAM 40115748224

GENOME=dasNov3_spike.fa
GENE=Dasypus_novemcinctus.Dasnov3.0.95_mod.spike_cd4.gtf

STAR --genomeSAindexNbases 14 --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $GENOME --sjdbGTFfile $GENE --limitGenomeGenerateRAM 40115748224
