---
title: "R Notebook"
output: html_notebook
---
 
#Reference genomes  
# Ignoring spike-ins
```{r}
files = as.character(unlist(read.table("runs") ))
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"
load("gene_annotations_v0.95.Rdata")

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



# Spike in version 
```{r}
files = as.character(unlist(read.table("runs") ))
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"

load("Y:/armadillo/ensembl/dasNov3.0.95/gene_annotations_v0.95_spikeins.Rdata")
load("Y:/armadillo/ensembl/ercc.conc.Rdata")

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

# Plot spike-ins

```{r}
load("U:/armadillo/reference_genomes/outs/armadillo_ref_counts.Rdata") 
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
 
 
# Stability 
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
 
  