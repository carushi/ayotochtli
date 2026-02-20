---
title: "R Notebook"
output: html_notebook
---

# Download data:

- Project here: [PRJNA595587](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA595587)
- DNA: [PRJNA591897](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA591897)
- RNA: [PRJNA595370](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA595370) [SRP237365](https://www.ncbi.nlm.nih.gov/sra/?term=SRP237365) [GSE141951](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141951)
- Virus: [AL450380](https://www.ncbi.nlm.nih.gov/nuccore/AL450380)


# Reference genomes  
## Setup genome index
Create genome file with spike-ins 
```
cat dasNov3.fa spike.fa > dasNov3_spike.fa 
```

```


```
~/STAR/STAR_2.7  \
--runThreadN 4       \
--runMode genomeGenerate  \
--genomeDir /dasNov3.0.95/  \
--genomeFastaFiles dasNov3_spike.fa  \
--sjdbGTFfile Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf    \
--sjdbOverhang 100                 \
--limitGenomeGenerateRAM 40115748224
```

## Run STAR 
```
~/STAR/STAR_2.7  \
--genomeDir dasNov3.0.95/ \
--runThreadN 12 \
--readFilesIn *1.fastq *2.fastq  \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts 
```

## Run STAR with 


## Ignoring spike-ins
```{r}
files = as.character(unlist(read.table("runs") ))
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"
load("gene_annotations_v0_95_mod.Rdata")

dir = "outs/"

                Ns = list()
                i = 1

                for( n in files ){
                        N = list()
                        filedir = paste(dir, n, sep="/")
                        countfile = paste(filedir, genecounts, sep=".")
                        logfile = paste(filedir, logname, sep=".")

                        if( file.exists(countfile) ) {
                                print(countfile)
                                counts =  read.table(countfile)


                                log1 =read.table(logfile, sep="\t", nrows=6)
                                log2 =read.table(logfile, sep="\t", skip=8, nrows=14)
                                log3 =read.table(logfile, sep="\t", skip=23, nrows=4)
                                log4 =read.table(logfile, sep="\t", skip=28, nrows=3)

                                N$mapinfo = rbind(log1,log2,log3,log4)
                                N$unmapped =  counts[1,]
                                N$multimapping = counts[2,]
                                N$noFeature =   counts[3,]
                                N$ambiguous = counts[4,]
                                N$length = dim(counts)[1]-4
                                N$genes = counts[ (1:N$length)+4,1]
                                N$counts1 = counts[ (1:N$length)+4,2]
                                N$counts2 = counts[ (1:N$length)+4,3]
                                N$counts3 = counts[ (1:N$length)+4,4]


                        } else {

                                 N$counts3 = rep(0, length(attr$ensemblID ) )

                        }
                        if( i > 1  ){
                                counts_exp = cbind(counts_exp, N$counts3)
                        } else {
                                counts_exp = N$counts3
                        }
                        Ns[[i]] = N
                        print(i)
                        i = i + 1

                }

                rownames(counts_exp) = attr$ensemblID
                colnames(counts_exp) = files
                save(Ns, counts_exp, file=paste(dir, "/counts.Rdata", sep=""))

```



## Including spike ins
```{r}
files = as.character(unlist(read.table("runs") ))
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"

load("gene_annotations_v0.95_spikeins.Rdata")
load("ercc.conc.Rdata")

dir = "outs/"

                Ns = list()
                i = 1

                for( n in files ){
                        N = list()
                        filedir = paste(dir, n, sep="/")
                        countfile = paste(filedir, genecounts, sep=".")
                        logfile = paste(filedir, logname, sep=".")

                        if( file.exists(countfile) ) {
                                print(countfile)
                                counts =  read.table(countfile)


                                log1 =read.table(logfile, sep="\t", nrows=6)
                                log2 =read.table(logfile, sep="\t", skip=8, nrows=14)
                                log3 =read.table(logfile, sep="\t", skip=23, nrows=4)
                                log4 =read.table(logfile, sep="\t", skip=28, nrows=3)

                                N$mapinfo = rbind(log1,log2,log3,log4)
                                N$unmapped =  counts[1,]
                                N$multimapping = counts[2,]
                                N$noFeature =   counts[3,]
                                N$ambiguous = counts[4,]
                                N$length = dim(counts)[1]-4
                                N$genes = counts[ (1:N$length)+4,1]
                                N$counts1 = counts[ (1:N$length)+4,2]
                                N$counts2 = counts[ (1:N$length)+4,3]
                                N$counts3 = counts[ (1:N$length)+4,4]


                        } else {

                                 N$counts3 = rep(0, length(attr$ensemblID ) )

                        }
                        if( i > 1  ){
                                counts_exp = cbind(counts_exp, N$counts3)
                        } else {
                                counts_exp = N$counts3
                        }
                        Ns[[i]] = N
                        print(i)
                        i = i + 1

                }

                rownames(counts_exp) = attr$ensemblID
                colnames(counts_exp) = files
                save(Ns, counts_exp, file=paste(dir, "/counts.Rdata", sep=""))
                save(counts_exp, file=paste(dir, "/armadillo_ref_counts.Rdata", sep=""))
                

```

### Plotting spike-ins

```{r}
load('armadillo_ref_counts.Rdata") 
X.cpm = calc_cpm(counts_exp)
f.a = (attr$assembly=="spikein")
m = match(attr[f.a,1], concentrations[,2] ) 
f.o = !is.na(m)
f.con = m[f.o]
o = order(concentrations[f.con ,4] ) 

boxplot( t(log2(1+X.cpm[f.a,][f.o,][o,]) ) , pch=19, col=magma(10)[6],xlab="Spikeins", ylab="Expression (log2 1 + CPM)") 
text( 1:length(o), -2, attr[f.a,][f.o,][o,1], xpd=-1, srt = 90 )
boxplot( t(log10(counts_exp[f.a,][f.o,][o,]) ),pch=19, col=makeTransparent(1), xlab="Spikeins" , ylab="Expression (log10 counts)") 

```
 
 
### Stability analysis  
```{r}
# i # gene
# j # dataset
# k # value/cell
# g # gene set
# c # expression data
c = X.cpm
e = f.a

sercc = colSums(c[e,] )
rERCC = sapply(which(!e), function(i) cor( c[i,] , sercc , m="s") )

hist(rERCC[!f.a][f.zz], border=NA, col=viridis(10)[3], main="", xlab="Absolute stability - correlations")


f.old = (attr$assembly=="ensembl")


load("U:/armadillo/functionalsets.Rdata")
# Proportional stability
g = f.old
g2 = list() 
sg = list() 
rSG = list() 

for(i in 1:length(functionalsetnames)) { 
  g2[[i]] = functionalsets[,i]== 1
  sg[[i]] = colSums(c[g,][g2[[i]],] )
  rSG[[i]] = sapply(which(g), function(ii) cor( c[ii,] , sg[[i]] , m="s") )
} 

i = 22 
hist(rSG[[i]][f.zz & !g2[[i]]], col=3, freq=F, main=functionalsetnames[i], xlab="Proportional stability - correlations")
hist(rSG[[i]][f.zz & g2[[i]]], add=T, col=2, freq=F )

```


# Personal genomes 
## Setup
- For each quad, download their personal genomes. 
- Make sure you have samtools, GATK (v3). 
```
picard='~/GenomeAnalysisTK/picard.jar'
GATK='~/GenomeAnalysisTK/GenomeAnalysisTK.jar'
igvtools='~/IGVTools/igvtools.jar'
```



## Call and filter for variants
Using quad15-50 as an example. Replace with different quads. 
### Call variants 
```
~/gatk-4.1.0.0/gatk HaplotypeCaller \
   -R dasNov3.fa \
   -I 15F501-15F502-15F503-15F504.sort.dedup.realign.bam \
   -O 15F501-15F502-15F503-15F504.raw_variants_nong.vcf

~/gatk-4.1.0.0/gatk SelectVariants \
-R dasNov3.fa  \
-V 15F501-15F502-15F503-15F504.raw_variants_nong.vcf  \
-O raw_snps.vcf \
--select-type-to-include SNP

~/gatk-4.1.0.0/gatk SelectVariants \
-R dasNov3.fa  \
-V 15F501-15F502-15F503-15F504.raw_variants_nong.vcf \
-O raw_indels.vcf \
--select-type-to-include INDEL

```
### Filter 
```
~/gatk-4.1.0.0/gatk VariantFiltration \
-R dasNov3.fa  \
-V raw_indels.vcf \
-O filtered_indels.vcf \
--filter-name "basic_snp_filter1" \
--filter-expression "QD < 2.0" \
--filter-name "basic_snp_filter2" \
--filter-expression "FS > 200.0" \
--filter-name "basic_snp_filter3" \
--filter-expression "ReadPosRankSum < -20.0" \
--filter-name "basic_snp_filter4" \
--filter-expression "SOR > 10.0" 
```

```
~/gatk-4.1.0.0/gatk VariantFiltration \
-R dasNov3.fa  \
-V raw_snps.vcf \
-O filtered_snps.vcf \
--filter-name "basic_snp_filter1" \
--filter-expression "QD < 2.0" \
--filter-name "basic_snp_filter2" \
--filter-expression "FS > 60.0" \
--filter-name "basic_snp_filter3" \
--filter-expression "MQ < 40.0" \
--filter-name "basic_snp_filter4" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "basic_snp_filter5" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "basic_snp_filter6" \
--filter-expression "SOR > 4.0" 
```
 


### Select homozygous alts 
```
zgrep '^#'  filtered_snps.vcf.gz > header.vcf
zgrep PASS filtered_snps.vcf.gz | grep '1/1' > filtered_snps.pass.vcf
 cat header.vcf filtered_snps.pass.vcf > filtered_snps.pass2.vcf
mv filtered_snps.pass2.vcf filtered_snps.pass.vcf

zgrep '^#'  filtered_indels.vcf.gz > header2.vcf
zgrep PASS filtered_indels.vcf.gz | grep '1/1' > filtered_indels.pass.vcf
 cat header2.vcf filtered_indels.pass.vcf > filtered_indels.pass2.vcf
mv filtered_indels.pass2.vcf filtered_indels.pass.vcf

```

```zgrep filtered_snps.vcf.gz | awk '{OFS="\t"; if (!/^#/){print $1,$2,$2,$4"/"$5,".","+","."}}' > 15-50.bed
```

```~/g2gtools convert -i 15-50.bed -c 15-50.vci.gz -o 15-50.conv.bed
```


### Count
```
grep '1/1' raw_indels.vcf -c
grep '1/1' raw_snps.vcf -c
grep '0/1' raw_indels.vcf -c
grep '0/1' raw_snps.vcf -c
zgrep '1/1' filtered_indels.vcf.gz  | wc -l 
zgrep '1/1' filtered_snps.vcf.gz  | wc -l 
zgrep '0/1' filtered_indels.vcf.gz | wc -l 
zgrep '0/1' filtered_snps.vcf.gz | wc -l 
```

 

## Convert reference to personal genome 
If you do not have enough memory for these tasks, you can split the genome into ~40 or so chunks and repeat on each part. It is a pain. 
```
##  Create VCI file of the sample
~/g2gtools-0.2.7/g2gtools vcf2vci -o 15-50.vci -p 16 -i  filtered_snps.pass.vcf.gz -i filtered_indels.pass.vcf.gz -s '15F501-15F502-15F503-15F504'  -f dasNov3.fa --pass --quality 


## Incorporate SNPs into the reference genome
~/g2gtools-0.2.7/bin/g2gtools patch   -p 16 -i dasNov3.fa  -c 15-50.vci.gz -o 15-50.patched.dasNov3.fa 

## Incorporate indels into the reference genome and get custom diploid genome for the sample.
~/g2gtools-0.2.7/bin/g2gtools transform -i 15-50.patched.dasNov3.fa  -c 15-50.vci.gz -o 15-50.dasNov3.fa -p 16

## Liftover gene annotation onto sample coordinates.
~/g2gtools-0.2.7/bin/g2gtools convert -i Dasypus_novemcinctus.Dasnov3.0.95_test.gtf -c 15-50.vci.gz -o 15-50.gtf              
cat 15-50.gtf 15-50.gtf.unmapped > 15-50.merged.gtf
```


In R, generate the gene annotation file for the quad: 
```{r}
id = "quad15-50"
id2 = "15-50"

setwd("../")
setwd(id)
 
temp = read.table( paste0(id2,".merged.gtf"),  header=F, sep="\t")
filt = temp[,3] == "gene"
mat = strsplit( as.character(temp[filt,9]), ";")

ids = sapply(1:length(mat), function(i) mat[[i]][1] )
ids = gsub( "gene_id ", "", ids)

checks = grep('gene_name', mat)
genename = rep ("", length(ids))
for ( i in checks) {
	grep('gene_name', mat[[i]])
	genename[i] = mat[[i]][grep('gene_name', mat[[i]])]
}
genename = gsub( " gene_name ", "", genename)


checks = grep('gene_biotype', mat)
genetype = rep ("", length(ids))
for ( i in checks) {
	genetype[i] = mat[[i]][grep('gene_biotype', mat[[i]])]
}
genetype = gsub( " gene_biotype ", "", genetype)


part1 = cbind( temp[filt,1:8], ids , genetype, genename)
part1[,2] = id 
part1 = part1[,-3]
part1 = part1[,-5]
part1 = part1[,-6]

colnames(part1) = c("scaffold", "quad", "start","end",  "strand", "ensemblID", "type",  "name" )
attr = part1
save(attr, file=paste0("gene_annotations_v0.95.", id, ".Rdata")) 
```



## Map to personal genomes 
First, generate personal genome index: 
```
~/STAR/STAR_2.7  \
--runThreadN 4       \
--runMode genomeGenerate  \
--genomeDir  quad15-50/   \
--genomeFastaFiles 15-50.trans.fa  \
--sjdbGTFfile 15-50.gtf  \
--sjdbOverhang 100                 \
--limitGenomeGenerateRAM 40048639360
```
Then map: 
```
~/STAR/STAR_2.7  \
--genomeDir dasNov3.0.95/quad15-50/ \
--runThreadN 12 \
--readFilesIn *1.fastq *2.fastq  \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts 
```

## Allele-specific expression with personal genome alignments 

```
genome=~/quad15-10/15-10.trans.fa

samtools faidx 15-50.trans.fa
java -jar ~/GenomeAnalysisTK/picard.jar CreateSequenceDictionary R=15-50.trans.fa O=15-50.trans.dict
```

```
echo "Adding read groups to bam file"
 java -jar  $picard AddOrReplaceReadGroups \
    I=Aligned.sorted.filt.bam  \
    O=Aligned.sorted.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

 echo "Marking duplicates"
java -jar $picard MarkDuplicates \
    I=Aligned.sorted.rg.bam  \
    O=Aligned.sorted.dedupped.bam  \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "Splitting and trimming"
java -jar $GATK \
   -T SplitNCigarReads         \
   -R $genome \
   -I Aligned.sorted.dedupped.bam \
   -o Aligned.sorted.split.bam  \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

echo "Counting SNPs"
java -jar $GATK \
   -T ASEReadCounter \
   -R quad15-50/15-50.trans.fa \
   -o test.counts.csv \
   -I Aligned.sorted.split.bam \
   -sites 15-50.bed.vcf \
   -U ALLOW_N_CIGAR_READS
```



