## Data

- Reference genome (dasNov3)
  - https://www.ensembl.org/Dasypus_novemcinctus/Info/Index
  - ../data/genome/
- Personal genome
  - 16-90.trans.fa
  - 16-90.merged.gtf
- Spikein
  - dasNov3_spike.fa
  - Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf
- scATAC-seq raw reads
  - ../data/GSE273197/
- 10X barcodes whitelist
  - 737K-cratac-v1.txt
- see rna_seq_processing.md for upstream analyses
- see scRNA_scanpy.ipynb for scRNA-seq read preprocessing

## Scripts

* bash genome_database.sh
* alignment_scatac_by_star.sh
  * --alignIntronMax 1 --alignMatesGapMax 1000
  * ref: [https://github.com/alexdobin/STAR/issues/558](https://github.com/alexdobin/STAR/issues/558)
  * default alignMatesGapMax: 0
* alignment_scatac_by_star_relaxed.sh (not applied in this study)
  * --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 15 \
  * default outFilterMismatchNmax: 10