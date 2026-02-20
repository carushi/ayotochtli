---
title: "Armadillo - human comparisons"
output: html_notebook
---
 
 
```{r}


conv_smoother2 <- function (X, Y, window, xlab = "", ylab = "", raw = FALSE, xlim, ylim, se_flag = FALSE)
{
  filt <- is.finite(X) & is.finite(Y) & !is.na(X) & !is.na(Y)
  X <- X[filt]
  Y <- Y[filt]
  n <- order(X)
  m <- X[n]
  i <- length(m)/window
  ymax <- max(Y)
  ymin <- min(Y)
  xmax <- max(X)
  xmin <- min(X)
  X_c <- m[(((1:i) * window) - window/2)]
  Y_c <- convolve(Y[n], rep(1, window), type = "filter")
  Y_c <- Y_c[(1:i) * window - (window - 1)]/window
  var_Y_c <- abs(convolve(Y[n]^2, rep(1, window)/window, type = "filter") -
                   (convolve(Y[n], rep(1, window)/window, type = "filter"))^2)
  std_Y_c <- var_Y_c^(1/2)
  std_Y_c <- std_Y_c[(1:i) * window - (window - 1)]
  se_Y_C <- se(std_Y_c)
    

  plot(X_c, Y_c, ylim = ylim, xlim = xlim,
       lwd = 4, type = "l", col = 0, bty = "n", xlab = xlab,
       ylab = ylab, cex.lab = 1.5, cex.axis = 1.2)
  if(se_flag == TRUE) { 
    polygon(c(X_c, rev(X_c)), c(Y_c - se_Y_C, rev(Y_c + se_Y_C)),
            col = "lightgrey", border = NA)  
  } else { 
    polygon(c(X_c, rev(X_c)), c(Y_c - std_Y_c, rev(Y_c + std_Y_c)),
          col = "lightgrey", border = NA)
  }
  if (raw == TRUE) {
    points(X, Y, col = rgb(75, 0, 150, 10, maxColorValue = 255),
           pch = 19, cex = 0.1)
  }
  lines(X_c, Y_c, col = 1, lwd = 2)
  smoothed <- cbind(X_c, Y_c, std_Y_c)
  return(smoothed)
}

```
 
 
# Armadillo fold changes
```{r}
 
fcs.arm = list()
for(j in 1:5){
  X.ms = sapply(1:4, function(i) rowMeans(X.cpm.all[ ,pData$Quad==j & pData$Ind ==i] ))
  X.fcs = sapply(1:4, function(i) sapply(i:4, function(ii) abs(log2(X.ms[,i]/X.ms[,ii] ) )))
  X.fcs = do.call(cbind, X.fcs)
  X.fcs = X.fcs[, -c(1,5,8,10)]
  fcs.arm[[j]] =  X.fcs 
  #fcs.arm[[j]] = rowMeans(X.fcs, na.rm=T)
  
}
fcs.arm = do.call(cbind, fcs.arm)
fcs.avg = rowMeans(fcs.arm[f.zz,], na.rm=T)
hist(fcs.avg)


fcs.arm.cross = list()
for(j in 1:5){
 X.ms.j = sapply(1:4, function(i) rowMeans(X.cpm.all[ ,pData$Quad==j & pData$Ind ==i] ))
 X.ms = list()
 for(ji in (1:5)[-j]){
   X.ms[[ji]] = sapply(1:4, function(i) rowMeans(X.cpm.all[ ,pData$Quad==ji & pData$Ind ==i] ))
 }
 X.ms = do.call(cbind, X.ms)
 nnn = dim(X.ms)[2]
 X.fcs = lapply(1:4, function(i) sapply(1:nnn, function(ii) abs(log2(X.ms.j[,i]/X.ms[,ii] ) )))
 X.fcs = do.call(cbind, X.fcs)
 fcs.arm.cross[[j]] =  X.fcs
 }
fcs.arm.cross = do.call(cbind, fcs.arm.cross)
fcs.arm.cross[!is.finite(fcs.arm.cross)] = NA


hist(fcs.arm.cross[f.zz,][f.arm,])
hist(fcs.arm[f.zz,][f.arm,])
AA = fcs.arm[f.zz,]
BB = fcs.arm.cross[f.zz,]
AAh = fcs.arm[f.zz,][f.arm,][f.na,]
BBh =  fcs.arm.cross[f.zz,][f.arm,][f.na,]
hb = hist(BB, breaks=100 )
ha = hist(AA, breaks=hb$breaks )
hc = hist(AAh, breaks=hb$breaks )
hd = hist(BBh, breaks=100 )

ha$counts = ha$counts/dim(AA)[2]
hb$counts = hb$counts/dim(BB)[2]
hc$counts = hc$counts/dim(AAh)[2]
hd$counts = hd$counts/dim(BBh)[2]
xrange = c(0,4)
yrange = c(0,2000)
xlab="|log2 FCs|"


par(mfrow=c(2,1))
plot(hc, xlim=xrange, ylim=yrange, border=NA, col=magma(10)[6], main="Within", xlab=xlab)
abline(v=mean(AAh,na.rm=T))
plot(hd, xlim=xrange, ylim=yrange, border=NA, col=magma(10)[5], main="Cross", xlab=xlab)
abline(v=mean(BBh,na.rm=T))


```
 
 
# Human fold changes
```{r}
# quantify RNA expression for TwinsUK data https://www.omicsdi.org/dataset/ega/EGAD00001001088
load("countsQC.Rdata")
load("cpm.QC.Rdata")

load("../data/genome/gene_annotations_v19.Rdata")
attr.human = attr
samples.cor = cor(counts, m="s")

consTree = hclust( dist(samples.cor) )
consDend = as.dendrogram(consTree)
bad = which(cutree(consTree,3) != 1  )
X.cpm = calc_cpm(counts)


dens.list = lapply(1:391, function(i) get_density(hist( log2(X.cpm[,i], plot=F ) ) )  ) 

plot( dens.list[[1]], type="l", lwd=2, col=makeTransparent(1), xlim=c(-5,20), ylim=c(0,0.3))
tmp = sapply((1:391)[-bad], function(i) lines(dens.list[[i]], col=makeTransparent(1)) )
tmp = sapply(63, function(i) lines(dens.list[[i]], col=makeTransparent(3), lwd=2) )


plot( dens.list[[110]], type="l", lwd=2, col=makeTransparent(1))
tmp = sapply(bad, function(i) lines(dens.list[[i]], col=makeTransparent(2), lwd=2) )


load("../data/genome/gene_annotations_v0.95_mod.Rdata")
attr.arm = attr

m =  match( rownames(X.cpm) , attr.human$ensemblID )
f.x = !is.na(m)
f.a = m[f.x]

m = match (attr.arm$name[f.zz ], attr.human$name[f.a] )
f.arm = !is.na(m)
f.hum = m[f.arm]
X.m.arm = rowMeans(X.cpm.all[f.zz,][f.arm,], na.rm=T)
X.m.hum  = rowMeans(X.cpm[f.x,][f.hum,], na.rm=T)

plot(log2(1+X.m.arm), log2(1+X.m.hum), pch=19, xlab="Armadillo average expression (log2 1 + CPM)", ylab="Human (MZ) average expression (log2 1 + CPM)")


ids1 = unlist(twins.pair.uniq[,1])
ids2 = unlist(twins.pair.uniq[,2])
types = unlist(twins.pair.uniq[,3]) 
f.mz2 = unlist(twins.pair.uniq[,3]) == 2 
f.dz2 = unlist(twins.pair.uniq[,3]) == 1 

f.mz2[!is.na(match(ids2, bad)) ] = F 
f.mz2[!is.na(match(ids1, bad)) ] = F 
f.dz2[!is.na(match(ids1, bad)) ] = F 
f.dz2[!is.na(match(ids2, bad)) ] = F



fcs.hum = list()
for(i in 1:length(ids1)){
  X.fcs = abs(log2(X.cpm[,ids1[i]]/X.cpm[,ids2[i]] )) 
  fcs.hum[[i]] = X.fcs
}
fcs.hum = do.call(cbind, fcs.hum)
fcs.avg.hum = rowMeans(fcs.hum[,f.mz2], na.rm=T)
hist(fcs.avg.hum)
fcs.avg.hum[!is.finite(fcs.avg.hum)] = NA
f.na = !is.na(fcs.avg.hum[f.x][f.hum])
fcs.hum[!is.finite(fcs.hum)] = NA



fcs.arm[!is.finite(fcs.arm)] = NA
A = fcs.arm[f.zz,][f.arm,][f.na,]
B = fcs.hum[f.x,][f.hum,][f.na,f.mz2]
C = fcs.hum[f.x,][f.hum,][f.na,f.dz2]

 
 f.good = (f.mz2 | f.dz2)
 fcs.hum.rand = list()
 s1 = sample( c(ids1[f.good],ids2[f.good]), 1000, replace=T)
 s2 = sample( c(ids2[f.good],ids1[f.good]), 1000, replace=T)
 redo = which(s1==s2)
 s1[redo] = sample( ids1[f.good], length(redo))
 
 for(i in 1:1000){
    X.fcs = abs(log2(X.cpm[,s1[i]]/X.cpm[,s2[i]] )) 
   fcs.hum.rand[[i]] = X.fcs
 }
 fcs.hum.rand = do.call(cbind, fcs.hum.rand)
 fcs.hum.rand[!is.finite(fcs.hum.rand)] = NA
 E = fcs.hum.rand[f.x,][f.hum,][f.na,]
  
 
yrange = c(0,2000)
xrange = c(0,4)

 ha = hist(A, breaks=100 )
 hab = hist(BB, breaks=100 )
 
 hb = hist(B, breaks=100 )
 hc = hist(C, breaks=100 )
 he = hist(E, breaks=100 )
 ha$counts = ha$counts/dim(A)[2]
 hb$counts = hb$counts/dim(B)[2]
 hc$counts = hc$counts/dim(C)[2]
 he$counts = he$counts/dim(E)[2]


 par(mfrow=c(3,1))
 plot(hb, xlim=xrange,  ylim=yrange,  border=NA, col=magma(10)[7], main="Human (MZ) ", xlab="|log2 FCs|")
 abline(v=mean(B,na.rm=T))
 abline(v=median(B,na.rm=T))
 
 plot(hc, xlim=xrange,  ylim=yrange,  border=NA,  col=magma(10)[8], main="Human (DZ)", xlab="|log2 FCs|") 
 abline(v=mean(C,na.rm=T))
 abline(v=median(C,na.rm=T))

 plot(he, xlim=xrange, ylim=yrange,  border=NA, col=magma(10)[9], main="Human (shuffled) ", xlab="|log2 FCs|")
abline(v=mean(E,na.rm=T))
abline(v=median(E,na.rm=T),lty=2)


```
 
 
# DE prior comparisons 
```{r}
load("../data/genome/gene_annotations_v25.Rdata")
m = match(rownames(B), attr.human$ensemblID)
f.b = !is.na(m)
f.h = m[f.b]
delist = read.table("../data/DE_prior/DE_prior.txt", header=T)
m = match(delist[,5], attr.human$name[f.h])
f.d = !is.na(m)
f.dd = m[f.d] 

rerank = rank(delist[f.d,4])
delist.r=round(rerank/length(rerank),1)


beanplot(A[f.b,][f.dd,] ~ delist.r, log="", what=c(1,1,1,0)) 
beanplot(B[f.b,][f.dd,] ~ delist.r, log="", what=c(1,1,1,0))  
beanplot(C[f.b,][f.dd,] ~ delist.r, log="", what=c(1,1,1,0))  

xlab= "DE prior rank (binned)"
ylab = "|log2 FCs|"
par(mfrow=c(1,3))
boxplot(A[f.b,][f.dd,] ~ delist.r, ylim=xrange, main="Armadillo", ylab=ylab, xlab=xlab,  col= magma(10)[6] ) 
boxplot(B[f.b,][f.dd,] ~ delist.r, ylim=xrange, main="Human - MZ" , ylab=ylab, xlab=xlab ,  col= magma(10)[7] )     
boxplot(C[f.b,][f.dd,] ~ delist.r, ylim=xrange, main="Human - DZ" , ylab=ylab , xlab=xlab ,  col= magma(10)[8] )    



m = match(delist[,5], attr.human$name[f.h])
f.d = !is.na(m)
f.dd = m[f.d]
rerank = rank(delist[f.d,4])
delist.r=round(rerank/length(rerank),1)
rerank = rerank/length(rerank)

repA = array(A[f.b,][f.dd,]  )
repDE = rep(rerank, dim(A)[2] )
f.na = !is.na(repA)

EGAD::conv_smoother( repDE[f.na], repA[f.na], 1000, xlab, ylab)
# SE 
conv_smoother2( repDE[f.na], repA[f.na], 1000, xlab, ylab, se_flag= TRUE)






plot ( rowMeans(AAh), rowMeans(BBh), pch=19)
abline(0,1,col="lightgrey", lwd=2)


pdf("decomp.pdf")
for(ii in c(0.9, 0.9, 0.95) ){
  filt =  rerank > ii
x =  rowMeans(AAh)[f.b][f.dd] 
y =  rowMeans(BBh)[f.b][f.dd]

zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))

xlab= "Mean |log2 FC| - within"
ylab= "Mean |log2 FC| - across"
xhist0 = hist(x, plot=FALSE)
yhist0 = hist( y, plot=FALSE)
xhist1 = hist(x[filt], plot=FALSE, breaks=xhist0$breaks)
yhist1 = hist( y[filt], plot=FALSE, breaks=yhist0$breaks)
yhist1$counts = yhist1$counts/sum(yhist1$counts)
yhist0$counts = yhist0$counts/sum(yhist0$counts)
xhist1$counts = xhist1$counts/sum(xhist1$counts)
xhist0$counts = xhist0$counts/sum(xhist0$counts)

top = max(c(xhist1$counts, xhist0$counts, yhist0$counts, yhist1$counts))
par(mar=c(3,3,1,1))
plot(x,y, pch=19, main="", col=makeTransparent(1,50))
points( x[filt] , y[filt], pch=19 , col=4, cex=1.2) 

par(mar=c(0,3,1,1))
barplot(xhist0$counts, axes=FALSE, ylim=c(0, top), space=0)
axis(2)
barplot(xhist1$counts, axes=FALSE, ylim=c(0, top), space=0, add=T, col=makeTransparent(4)) 

par(mar=c(3,0,1,1))
barplot(yhist0$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
axis(1)
barplot(yhist1$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, add=T, col=makeTransparent(4))

par(oma=c(3,3,0,0))
mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
      at=.5 * (mean(x) - min(x))/(max(x)-min(x)))
mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
      at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

dev.off() 
```

# eGene recurrence comparisons 
```{r}
load("eqtl.mat.Rdata")

m = match(rownames(B), rownames(eqtl.mat) )
f.eb = !is.na(m)
f.eq = m[f.eb]
recur.eq = rowSums(eqtl.mat[f.eq,c(1:4,8,9,10)])

hist(recur.eq[recur.eq>0], col=magma(5)[3], main="", xlab="Recurrence",  breaks=10 )
fdrs = calc_fdrs_recur(eqtl.mat[f.eq,c(1:4,8,9,10)])

recur.eq2 = rowSums(eqtl.mat[,c(1:4,8,9,10)])

fdrs2 = calc_fdrs_recur(eqtl.mat[,c(1:4,8,9,10)])
hist( recur.eq2[recur.eq2>0], main="", xlab="Recurrence", col=4, breaks=10 )
abline(v=fdrs2$Pt)

pdf("eqtl.recurrence.pdf")
hist( recur.eq2[recur.eq2>0], main="", xlab="Recurrence", col=4, breaks=10 )
abline(v=fdrs2$Pt)
hist(recur.eq[recur.eq>0], col=magma(5)[3], main="", xlab="Recurrence",  breaks=10 )
abline(v=fdrs$Pt)
dev.off() 

 
xlab="eGene recurrence"
par(mfrow=c(1,3))
boxplot(A[f.eb,] ~ recur.eq , ylim=xrange, main="Armadillo", ylab=ylab, xlab=xlab,  col= magma(10)[6] ) 
boxplot(B[f.eb,] ~ recur.eq , ylim=xrange, main="Human - MZ", ylab=ylab, xlab=xlab,  col= magma(10)[7] ) 
boxplot(C[f.eb,] ~ recur.eq , ylim=xrange, main="Human - DZ", ylab=ylab, xlab=xlab,  col= magma(10)[8] ) 
 
#save(recur.eq, f.eb, f.eq, file="eGene.recur.Rdata")
# 

repA = array(A[f.eb,]  )
repDE = rep(recur.eq, dim(A)[2] )
f.na = !is.na(repA)

# SD 
EGAD::conv_smoother( repDE[f.na], repA[f.na], 1000, xlab, ylab)
# SE 
conv_smoother2( repDE[f.na], repA[f.na], 1000, xlab, ylab, se_flag= TRUE)





plot ( rowMeans(AAh), rowMeans(BBh), pch=19)
abline(0,1,col="lightgrey", lwd=2)

pdf("egenescomp.pdf")


for( ii in c(2,2, 4,7)){  
  filt = recur.eq >= ii
x =  rowMeans(AAh)[f.eb] 
y =  rowMeans(BBh)[f.eb]
zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))

xlab= "Mean |log2 FC| - within"
ylab= "Mean |log2 FC| - across"
xhist0 = hist(x, plot=FALSE)
yhist0 = hist( y, plot=FALSE)
xhist1 = hist(x[filt], plot=FALSE, breaks=xhist0$breaks)
yhist1 = hist( y[filt], plot=FALSE, breaks=yhist0$breaks)
yhist1$counts = yhist1$counts/sum(yhist1$counts)
yhist0$counts = yhist0$counts/sum(yhist0$counts)
xhist1$counts = xhist1$counts/sum(xhist1$counts)
xhist0$counts = xhist0$counts/sum(xhist0$counts)

top = max(c(xhist1$counts, xhist0$counts, yhist0$counts, yhist1$counts))
par(mar=c(3,3,1,1))
plot(x,y, pch=19, main="", col=makeTransparent(1,50))
points( x[filt] , y[filt], pch=19 , col=2, cex=1.2) 
abline(0,1, col="grey", lwd=2)

par(mar=c(0,3,1,1))
barplot(xhist0$counts, axes=FALSE, ylim=c(0, top), space=0)
axis(2)
barplot(xhist1$counts, axes=FALSE, ylim=c(0, top), space=0, add=T, col=makeTransparent(2)) 

par(mar=c(3,0,1,1))
barplot(yhist0$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
axis(1)
barplot(yhist1$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, add=T, col=makeTransparent(2))

par(oma=c(3,3,0,0))
mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
      at=.5 * (mean(x) - min(x))/(max(x)-min(x)))
mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
      at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
} 

dev.off() 

```
 
