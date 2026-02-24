---
title: "Armadillo transcriptional identity"
output: html_notebook
---


# Load data
```{r}
source("../scripts/helper.r")
DIR <<- "../bulk_transcriptome/data/"
load(paste0(DIR, "cpm_strand_comb.Rdata"))
load(paste0(DIR, "counts_strand_comb.Rdata"))
load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
load(paste0(DIR, "armadillo_helper.Rdata"))
load(paste0(DIR, "armadillo_hem.Rdata"))
load(paste0(DIR, "ref.strand.test_train.Rdata"))
# Variables 
N = dim(X.cpm.all)[1]
n_quads = 5
n_times = 3 
n_qt = n_quads * n_times
n_q = 4 
n_samp = n_quads * n_q
r_samp = 1:n_samp 

# Labels 
quad = as.numeric(pData$Quad)
sex = as.numeric(pData$Sex)
lane = as.numeric(pData$Batches)
tlab = c("t1", "t2", "t3")
quads = unique(substr(pData$ID, 1,4))

labels = paste(sapply(1:n_quads, function(i) rep(quads[i],n_times)), tlab  ) 
sexlabels = unlist(lapply(1:n_quads, function(i) rep(unique(cbind(pData$Quad, pData$Sex) )[i,2],n_times)) ) 
timelabels = rep(1:n_times, n_quads)
quadlabels = unlist(lapply(1:n_quads, function(i) rep(unique(cbind(pData$Quad, pData$Sex) )[i,1],n_times))   ) 

# Z-scores per quad per timepoint 
X.zscores = lapply(1:n_qt, function(i)  t(scale( t(X.cpm.all[, (1:4) + 4*(i-1)  ]))))   
X.zscores = do.call(cbind, X.zscores)


# Ranks per quad per timepoint (across quads)
X.rank = lapply(1:n_qt, function(i)   t(apply(X.cpm.all[,(1:4) + 4*(i-1) ],  1, rank )))   
X.rank2 = do.call(cbind, X.rank)

# Gene filters 
f.z = sapply(1:n_qt, function(i) rowSums(X.cpm.all[,1:4+ 4*(i-1)]) > 0  )
f.zz = rowSums(f.z) == n_qt 
f.zzz = (rowSums(X.cpm.all == 0 )  == 0  )
```
 

# Transciptional analysis 
```{r}
sample.cors = cor(X.cpm.all, method ="spearman", use="pair")
heatmap.3(sample.cors[r_samp,r_samp], 
          col=cols5, 
          ColSideColors = matrix(candy_colors[pData$Quad[r_samp]], ncol=1),
          RowSideColors = matrix(cols15[pData$Sex[r_samp]], ncol=20))

heatmap.3(sample.cors[r_samp,r_samp], 
          col=cols5, 
          ColSideColors = matrix(candy_colors[pData$Quad[r_samp]], ncol=1),
          RowSideColors = matrix(candy_colors[-(1:3)][pData$Lanes[r_samp]], ncol=20))
```

# Between siblings
```{r}
i = 1 
j = 2 
plot( log2(1+X.cpm.all[,i]),  log2(1+X.cpm.all[,j]) , 
      pch=19,
      xlab=paste("Armadillo", pData$altID[i], "log2 (1+CPM)"),
      ylab=paste("Armadillo", pData$altID[j], "log2 (1+CPM)"))

```

# Across litters
```{r}
i = 1 
j = 18 
plot( log2(1+X.cpm.all[,i]),  log2(1+X.cpm.all[,j]) , 
      pch=19,
      xlab=paste("Armadillo", pData$altID[i], "log2 (1+CPM)"),
      ylab=paste("Armadillo", pData$altID[j], "log2 (1+CPM)"))

```


# Cross quad comparisons 
```{r}
# Filtered correlations
samples.cor2  = cor(X.cpm.all[f.zz ,], m="s") 

bottom = row(samples.cor2[r_samp,r_samp]) > col(samples.cor2[r_samp,r_samp])
bottom2 = row(samples.cor2[1:n_q,1:n_q]) > col(samples.cor2[1:n_q,1:n_q])

# per timepoint 
t1 = samples.cor2[r_samp,r_samp][bottom]
t2 = samples.cor2[r_samp+n_samp,r_samp+n_samp][bottom]
t3 = samples.cor2[r_samp+(n_samp*2),r_samp+(n_samp*2)][bottom]

# Within
sim.within = lapply(1:n_times, function(i) lapply(1:n_quads, function(j) samples.cor2[r_samp+n_samp*(i-1),r_samp+n_samp*(i-1)][ (1:n_q)+n_q*(j-1),(1:n_q)+n_q*(j-1)][bottom2]) ) 
# across   
sim.across = lapply(1:n_times, function(i) lapply(1:n_quads, function(j) samples.cor2[r_samp+n_samp*(i-1),r_samp+n_samp*(i-1)][ (1:n_q)+n_q*(j-1),-((1:n_q)+n_q*(j-1))]) ) 

sim.comb = lapply(1:n_times, function(i) 
  sapply(1:n_quads, function(j)  list(sim.within[[i]][[j]], sim.across[[i]][[j]])) )
sim.comb2 = lapply(1:n_times, function(i)   list(unlist(sim.within[[i]]), unlist(sim.across[[i]]) ) ) 
sim.comb3 = list(sim.comb2[[1]][[1]], sim.comb2[[1]][[2]], sim.comb2[[2]][[1]], sim.comb2[[2]][[2]], sim.comb2[[3]][[1]], sim.comb2[[3]][[2]] ) 

colslist = lapply(1:n_quads, function(j) original_colors[j])
colslist2 = lapply( sort(rep(1:n_quads,2)), function(j) original_colors[j] ) 
colslist3 = lapply( sort(rep(1:n_times,2)), function(j) mm_colors[j] ) 

# beanplot(sim.comb[[1]], sim.comb[[2]], sim.comb[[3]], col=colslist2)
# abline(v=c(10.5), lwd=3)
# abline(v=c(20.5), lwd=3)
# beanplot(sim.comb2[[1]], sim.comb2[[2]], sim.comb2[[3]], col=colslist3)
# vioplot.2(sim.comb3, colin=unlist(colslist3) ) 


sim.within.all = unlist(lapply(1:n_times, function(i) sim.comb2[[i]][[1]]) )
sim.across.all = unlist(lapply(1:n_times, function(i) sim.comb2[[i]][[2]]) )


h1 = hist(sim.across.all, breaks=30, plot=F) 
h2 = hist(sim.within.all, breaks=20, plot=F) 
h2$density = h2$counts/sum(h2$counts)
h1$density = h1$counts/sum(h1$counts)
ylim= range(c( h1$density, h2$density) )
xlim= range(c( h1$breaks, h2$breaks) )

plot(h1, freq=F, col=makeTransparent("blue", alpha=0.5), border=NA, xlim=xlim, ylim=ylim, xlab="Sample-sample correlations", main="")
plot(h2, freq=F,add=T, border=NA, col=makeTransparent("purple", alpha=0.5) )
abline(v= mean(sim.across.all), lwd=2, col="blue" )
abline(v= mean(sim.within.all), lwd=2, col="purple" )
wilcox.test( sim.across.all, sim.within.all )
```

# Highly variable genes
```{r}
require(DESeq2)
require(statmod)


# Per quad per timepoint
var.fits = list()
for(i in 1:n_qt){
  ed <- X.cpm.all[f.zz,(1:n_q)+n_q*(i-1)] 
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  
  minMeanForFit <- unname( quantile( means[ which( cv2 > 1 ) ], .95 ) )
  useForFit <- means >= minMeanForFit # & spikeins
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  df <- ncol(ed) - 1
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio,decreasing=T)
  oed <- ed[varorder,]

  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  
  var.fits[[i]] = list(varorder,means, vars, cv2, adj.pval) 
  
}

plot_cv2 <- function(varorder, means, cv2, ...){
  minMeanForFit <- unname( quantile( means[ which( cv2 > 1 ) ], .95 ) )
  useForFit <- means >= minMeanForFit # & spikeins
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  df <- ncol(ed) - 1
  
  smoothScatter(log(means),log(cv2),colramp=colorRampPalette(c("white","purple","orange")), ... ) ; 
  lines( log(xg), log(vfit), col="black", lwd=3 ); 
  lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); 
  lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
  points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2, pch=19)
  points(log(means[rev(varorder)[1:100]]),log(cv2[rev(varorder)[1:100]]),col=4, pch=19)
  
}

for(i in 1:n_qt){
   plot_cv2( var.fits[[i]][[1]], var.fits[[i]][[2]], var.fits[[i]][[4]] , main="")
}
  # save(var.fits, var.list, f.zz, file=paste0("hvgs.", dataset, ".Rdata") ) 

```

```{r}
n_hvgs = 100 

hvgs = matrix(NA, ncol=15, nrow=N)
# HVGS 
for(i in 1:n_qt){
  hvgs[f.zz,i][var.fits[[i]][[1]]] = (sum(f.zz):1)/sum(f.zz) 
}
rownames(hvgs) = attr$name


frac = 1 - n_hvgs/ sum(f.zz)
genesets  =  1*(hvgs[ ,1:n_qt]> frac)
rownames(genesets ) = attr$ensemblID
colnames(genesets ) = labels
recur.hvgs    = rowSums(1*(hvgs[ ,1:n_qt]> frac ) ) 
fdrs.hvgs     = calc_fdrs_recur(1*(hvgs[ ,1:n_qt]> frac) )

recur.hvgs.q  = sapply(1:n_quads, function(j) rowSums(genesets[f.zz,(1:n_times) + n_times*(j-1)]))
fdrs.hvgs.q   = calc_fdrs_recur( 1*(recur.hvgs.q > 0 ) )


recur.q  = rowSums(recur.hvgs.q  > 0)  

hist(recur.q[recur.q >0]+0.01, 
     border=NA, 
     col=magma(10)[3],
     main="", 
     xlab=paste("Recurrence of top", n_hvgs, "HVGs" ))
abline(v=fdrs.hvgs.q$Pt, lwd=2)
```

## HVGs cross-quad
```{r} 
frac = 1 - n_hvgs/ sum(f.zz)
genesets  =  1*(hvgs[ ,1:n_qt]> frac)
rownames(genesets) = attr$ensemblID
colnames(genesets) = labels

mf_cor = calculate_multifunc(genesets[f.zz,])
mf_cor.q = lapply(1:n_quads, function(j) calculate_multifunc(genesets[f.zz, -c((1:n_times) + n_times*(j-1))]))

aurocs.q =  sapply(1:n_qt, function(i) auc_multifunc( genesets[f.zz,], rank(hvgs[f.zz,i])) )
aurocs.q2 = aurocs.q
for(j in 1:n_quads){ 
  aurocs.q2[ (1:n_times) + n_times*(j-1) , (1:n_times) + n_times*(j-1) ] = NA 
} 

hist(aurocs.q2[row(aurocs.q2) > col(aurocs.q2)], xlim = c(0.5,1), freq=F, col=magma(10)[4], border=NA, xlab="AUROC HVGs", main="")
abline(v=0.5, col=2, lty=2, lwd=2)
abline(v= mean(aurocs.q2[row(aurocs.q2) > col(aurocs.q2)], na.rm=T),   lwd=2)

```

# Phenotypic comparisons 
```{r}
colquad <-  pData$Quad[match(hemotology[,1] , pData$altID) ]

quads2 = paste0(quads, "0")
o2 = c(2,1,3:5)
is = c(4,15,8,12,13,14)
par(mfrow=c(2,3))
for(i in is) { 
  beeswarm(  as.numeric(as.character(hemotology[,i])) ~  quads2[colquad] ,  
             pch=19 , xlab="Quadruplet", 
             ylab=colnames(hemotology)[i], 
             main=colnames(hemotology)[i], 
             col=candy_colors[o2]) 
 bxplot(as.numeric(as.character(hemotology[,i])) ~  quads2[colquad] , add=T)
} 

```


# ANOVA 
```{r}
aof <- function(x, indv, time) {
  m<-data.frame(indv,time, x);
  anova(aov(x ~ indv * time , m))
}

for(j in 1:n_quads) {
time = sort(rep(1:n_times, n_samp) )[pData$Quad == j]
indv = pData$Ind[pData$Quad == j]

anovaresults <- apply(X.cpm.all[,pData$Quad == j], 1, aof, indv, time)

#load(file=paste0("anova.quad",j,".Rdata")) 
temp.1 <- unlist(lapply(anovaresults, function(x) { x["Pr(>F)"][1,] }))
temp.2 <- matrix(temp.1, length(anovaresults),1, byrow=T)
dimnames(temp.2) <- list(names(anovaresults), dimnames(anovaresults[[1]])[[1]][1])
pvalues <- data.frame( t(temp.2))
pvalues = t(pvalues) 
p.adjust 
par(mfrow=c(2,2) )
hist( -log10(pvalues[f.zz,1]),  main="ANOVA p-values" ,xlab="- log10 p-value" , col=candy_colors[j], border=NA ) 
hist(-log10(p.adjust(pvalues[f.zz,1]) ), main="ANOVA q-values", xlab="-log10 q-value"  , col=candy_colors[j], border=NA) 

 
hist(  (pvalues[f.zz,1]),  main="ANOVA p-values" ,xlab=" p-value" , col=candy_colors[j], border=NA )
hist( (p.adjust(pvalues[f.zz,1]) ), main="ANOVA q-values", xlab="  q-value"  , col=candy_colors[j], border=NA) 


save(anovaresults, pvalues, file=paste0("anova.quad",j,".Rdata")) 

}
```

```{r}
pdf("anova.plots.pdf") 
for(j in 1:n_quads){ 
  nnj = min(as.numeric(pData$altID[pData$Quad==j]))  - 1 
  
  load(paste0("anova.quad",j,".Rdata"))
  
  quad.hi.p = data.frame(pvalues[f.zz,1][ which(pvalues[f.zz] < 0.05)])
  quad.hi.pdata <- X.cpm.all[ row.names(quad.hi.p), ]
  quad.hi.pzsco <- X.zscores[ row.names(quad.hi.p), ]

 
 heatmap.3( log2(1+quad.hi.pdata[,pData$Quad==j] ), 
            col=plasma(100),    
            ColSideColors=tropical_colors[pData$Ind[pData$Quad==j]],
            main=quads[j])


boxplot( t( log2(1+quad.hi.pdata[,pData$Quad==j]))  ~ as.character(pData$ID[pData$Quad==j] ),
         col=tropical_colors[sort(as.numeric(pData$altID[pData$Quad==j] ))-nnj])

bxplot( t( log2(1+quad.hi.pdata[,pData$Quad==j]))  ~ as.character(pData$ID[pData$Quad==j] ) )
beeswarm( t( log2(1+quad.hi.pdata[,pData$Quad==j]))  ~ as.character(pData$ID[pData$Quad==j] ), 
         col=tropical_colors[sort(as.numeric(pData$altID[pData$Quad==j] ))-nnj], 
         pch=19, corral="random", 
         cex=0.5, 
         add=T)


heatmap.3( (quad.hi.pzsco[,pData$Quad==j] ), 
           col=plasma(100), 
           ColSideColors=tropical_colors[pData$Ind[pData$Quad==j]])

beeswarm( t( (quad.hi.pzsco[,pData$Quad==j]))  ~ as.character(pData$ID[pData$Quad==j] ), 
          col=tropical_colors[sort(as.numeric(pData$altID[pData$Quad==j] ))-nnj], 
          pch=19, 
          corral="random", 
          cex=0.5,
          xlab="Individual per time point", 
          ylab="Z-score expression")
bxplot( t( (quad.hi.pzsco[,pData$Quad==j]))  ~ as.character(pData$ID[pData$Quad==j] ), add=T )


beeswarm( t( (quad.hi.pzsco[,pData$Quad==j]))  ~ as.character(pData$altID[pData$Quad==j] ), 
          col=tropical_colors[1:4], 
          pch=19, corral="random", 
          cex=0.5,
          xlab="Individual", 
          ylab="Z-score expression")
bxplot( t( (quad.hi.pzsco[,pData$Quad==j]))  ~ as.character(pData$altID[pData$Quad==j] ), add=T )

quad.hi.pzsco.list = lapply(1:4, function(i)  quad.hi.pzsco[,pData$Quad==j][ ,pData$Ind[pData$Quad==j]==i ]  ) 
vioplot.2( quad.hi.pzsco.list, 
           colin = tropical_colors[1:4],
          xlab="Individual", 
          ylab="Z-score expression")

} 
dev.off() 

```


# Identity analysis 
```{r}
# Testing 
X.tt1 = X.rank2[,r_samp]            # time point 1
X.tt2 = X.rank2[,r_samp+n_samp]     # time point 2
X.tt3 = X.rank2[,r_samp+(n_samp*2)] # time point 3

# Training
X.t1 = X.cpm.all[,-(r_samp)]           # time point 2 & 3
X.t2 = X.cpm.all[,-(r_samp+n_samp)]    # time point 1 & 3
X.t3 = X.cpm.all[,-(r_samp+(n_samp*2))] # time point 1 & 2

X.t1.q = lapply(1:n_quads, function(j) t(sapply(1:N, function(i) rank(rowMeans(cbind(X.t1[i,((1:n_q)+n_q*(j-1)) ], X.t1[i,(((1:n_q)+n_q*(j-1))+n_samp)]) )  ) ))) 
X.t2.q = lapply(1:n_quads, function(j) t(sapply(1:N, function(i) rank(rowMeans(cbind(X.t2[i,((1:n_q)+n_q*(j-1)) ], X.t2[i,(((1:n_q)+n_q*(j-1))+n_samp)]) )  ) ))) 
X.t3.q = lapply(1:n_quads, function(j) t(sapply(1:N, function(i) rank(rowMeans(cbind(X.t3[i,((1:n_q)+n_q*(j-1)) ], X.t3[i,(((1:n_q)+n_q*(j-1))+n_samp)]) )  ) ))) 

X.tr1 = do.call(cbind,X.t1.q )
X.tr2 = do.call(cbind,X.t2.q )
X.tr3 = do.call(cbind,X.t3.q )

colnames(X.tr1) = pData$altID[r_samp]
colnames(X.tr2) = pData$altID[r_samp]
colnames(X.tr3) = pData$altID[r_samp]

# Correlations between timepoints 
## Correlations between t1 and t2
cor.t1t2 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.cpm.all[i,(1:n_q)+j], X.cpm.all[i,((1:n_q)+n_samp )+j], m="s") )  )
## Correlations between t1 and t3
cor.t1t3 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.cpm.all[i,(1:n_q)+j], X.cpm.all[i,((1:n_q)+(n_samp*2))+j], m="s") )  )
## Correlations between t2 and t3 
cor.t2t3 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.cpm.all[i,((1:n_q)+(n_samp*2))+j], X.cpm.all[i,((1:n_q)+n_samp )+j], m="s") )  )


cor.all = do.call(cbind, lapply(1:n_quads, function(j) cbind( cor.t2t3[[j]], cor.t1t3[[j]], cor.t1t2[[j]]))) 

# Correlations between training and testing 
## Correlations between t1+t2 and t3
pred.t1t2 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.tr3[i,(1:n_q)+j], X.tt3[i,(1:n_q)+j], m="s") )  )
## Correlations between t1+t3 and t2
pred.t1t3 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.tr2[i,(1:n_q)+j], X.tt2[i,(1:n_q)+j], m="s") )  )
## Correlations between t2+t3 and t1 
pred.t2t3 = lapply( (0:n_q)*n_q, function(j) sapply(1:N, function(i) cor( X.tr1[i,(1:n_q)+j], X.tt1[i,(1:n_q)+j], m="s") )  )

predcor.all = do.call(cbind, lapply(1:n_quads, function(j) cbind( pred.t2t3[[j]], pred.t1t3[[j]], pred.t1t2[[j]]))) 

# Gene predictability scores
pred1 = sapply(1:n_quads, function(j) rowSums(t(sapply(1:N, function(k) (diag(sapply(((1:n_q)+n_q*(j-1)), function(i) (X.tr1[k,i] ==  X.tt1[k,((1:n_q)+n_q*(j-1))] ) * 1) ) ) ) )))  
pred2 = sapply(1:n_quads, function(j)  rowSums(t(sapply(1:N, function(k) (diag(sapply(((1:n_q)+n_q*(j-1)), function(i) (X.tr2[k,i] ==  X.tt2[k,((1:n_q)+n_q*(j-1))] ) * 1) ) ) ) )) )
pred3 = sapply(1:n_quads, function(j)  rowSums(t(sapply(1:N, function(k) (diag(sapply(((1:n_q)+n_q*(j-1)), function(i) (X.tr3[k,i] ==  X.tt3[k,((1:n_q)+n_q*(j-1))] ) * 1) ) ) ) )))
pred.all = do.call(cbind, lapply(1:n_quads, function(i) cbind(pred1[,i], pred2[,i], pred3[,i])) ) 


rownames(cor.all)     = rownames(X.cpm.all)
rownames(predcor.all) = rownames(X.cpm.all)
rownames(pred.all)    = rownames(X.cpm.all)


colnames(cor.all)     = labels
colnames(predcor.all) = labels
colnames(pred.all)    = labels

```

 
```{r}
for(j in 1:n_quads){
  x = rank(rowMeans(X.cpm.all[f.zz,pData$Quad == j]) ) 
  y = rowMeans(cor.all[f.zz,quadlabels==j])  
  filt = y > 0 
  EGAD::conv_smoother(x[filt], y[filt], 100)
  title(quads[j])
}

```


# Task 1 
## A: correlations
```{r}
nr = 1:20
nr2 = nr/20
xlab="Correlated genes threshold"


pred.scores.list = list() 
pval.scores.list = list() 
pred.indscores.list = list()
setsizes.scores.list  = list() 
genesets.scores.list = list() 
scorematrix.list = list() 

i = 1 
for(k in nr2){ 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  genesets.temp = list() 
  scorematrix.temp = list() 
  
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    q1 = which(cor.t2t3[[j]][f.zz]>=k & !is.na(cor.t2t3[[j]][f.zz]))
    q2 = which(cor.t1t3[[j]][f.zz]>=k & !is.na(cor.t1t3[[j]][f.zz])) 
    q3 = which(cor.t1t2[[j]][f.zz]>=k & !is.na(cor.t1t2[[j]][f.zz]))
    
    t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
    t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
    t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
    genesets.temp[[j]] = list( (q1),  (q2),  (q3))   
    scorematrix.temp[[j]] = list( t1[[1]], t2[[1]], t3[[1]] )   
    
  }
  pred.scores.list[[i]] = pred.temp
  pval.scores.list[[i]] = pval.temp
  pred.indscores.list[[i]] = pred.ind.temp
  setsizes.scores.list[[i]] = setsizes.temp
  genesets.scores.list[[i]] = genesets.temp
  scorematrix.list[[i]] = scorematrix.temp
  
  
  i = i + 1
}




cd = sapply(nr, function(i) mean(pred.scores.list[[i]]))
bc = sapply(nr, function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(nr, function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(nr, function(i) (pred.scores.list[[i]]))
de.pval = sapply(nr, function(i) (pval.scores.list[[i]]))

bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]

filename = "task1_results_corr_ref.strand.Rdata"; dataset = "ref.stranded"

save(
  pred.scores.list,
  pval.scores.list,
  pred.indscores.list,
  setsizes.scores.list, 
  genesets.scores.list,nr, nr2, 
  cd,ab,de,bc, cd.pval,ab.pval,de.pval,bc.pval, 
  file=filename)

k = 20
barplot(colMeans(pred.scores.list[[k]]), col=mm_colors); abline(h=1, col=4, lty=2)
barplot(rowMeans(pred.scores.list[[k]]), col=candy_colors); abline(h=1, col=4, lty=2)

heatmap.3(pred.scores.list[[k]], 
          Rowv=F, 
          Colv=F, 
          col=magma(5), 
          ColSideColors=mm_colors[1:3], 
          RowSideColors=candy_colors[1:5], 
          labRow=quads, 
          labCol=tlab)

 
```

```{r}
k = 20 
plot( pval.sum[,1], x15r/sum(x15r), pch=19, xlab="Average score", ylab="Density", type="o")
abline(v=cd[k], col=4, lwd=2)
abline(v=1, col=2, lty=2 )

# Fig 2D

plot( pval.sum, pch=19, xlab="Average score", ylab="Cummulative Density", type="o")
abline(v=cd[k], col=4, lty=2)
abline(h=cd.pval[k], col=4, lty=2)
abline(v=1, col=2, lty=2)
```

# Null - this takes time, run only once 
```{r, eval=FALSE }
filename="cpm_test_train.Rdata"; dataset ="cpm"

nr3 = unique(round(10^((1:430)/100)  ))
nr3 = nr3[nr3 < sum(f.zz)]

rand.q = list()  

for(rri in 1:length(nr3) )    { 
  rand.q2 = list()
  for(j in 1:5){
    ki = ((1:4)+4*(j-1))
    
    t1.rand = sapply(1:1000, function(ri) prediction_gene_scores(X.tr1[f.zz,ki], X.tt1[f.zz,ki], sample(sum(f.zz), nr3[rri]))[[2]]) 
    t2.rand = sapply(1:1000, function(ri) prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], sample(sum(f.zz), nr3[rri]))[[2]])  
    t3.rand = sapply(1:1000, function(ri) prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], sample(sum(f.zz), nr3[rri]))[[2]]) 
    rand.q2[[j]] = cbind(t1.rand, t2.rand, t3.rand)
  }
  rand.q[[rri]]= rowMeans(do.call(cbind, rand.q2))
  
} 

cr = sapply(1:length(nr3), function(i) mean(rand.q[[i]])) 
cr.sd = sapply(1:length(nr3), function(i) sd(rand.q[[i]]))
cr.se = sapply(1:length(nr3), function(i) se(rand.q[[i]]))

save(cr, cr.sd, cr.se, rand.q, nr3, file=paste0("random.sets.", dataset, ".Rdata")) 

```

# Plot

```{r}
# Fig 2E
load(paste0("random.sets.", dataset, ".Rdata") )

Y_c = cr[-1]
X_c = log10(nr3)[-1]
f.plot = (X_c>=1)

std_Y_c = cr.sd[-1][f.plot]
X_c = X_c[f.plot]
Y_c = Y_c[f.plot]
#X_c = (nr3)[-1]
# std_Y_c = cr.se[-1]

k = 20 
ymax = max( cd[k])
ymin  = min( Y_c - std_Y_c)
ylab = "Mean gene set average scores"
xlab = "Gene set size"
 
plot( X_c, Y_c, ylim = c(ymin, ymax), xlim = c(1,4.3), lwd = 4, type = "l", col = 0, bty = "n", xlab = xlab, ylab = ylab,  axes=F)
polygon(c(X_c, rev(X_c)), c(Y_c - std_Y_c, rev(Y_c + std_Y_c)), col = "lightgrey", border = NA)
#lines(X_c, Y_c, ylim = c(ymin, ymax),  lwd = 4)
points(X_c, Y_c, ylim = c(ymin, ymax) , pch=19, cex=0.5)
axis(2)
axis(1, at = (0:4), lab = 10^(0:4))
abline(h=cd[k], col=2)
abline(v=log10(mean(setsizes.scores.list[[k]])), col=2)

```



# Perfect predictors
```{r}
require(venn)

hist(pred.all[f.zz,], breaks=c(-1:4), col=magma(5), main="", xlab="Gene scores", sub=dataset)
hist(predcor.all[f.zz,], col=viridis(20), main="Testing", sub=dataset, xlab="Correlations")
hist(cor.all[f.zz,], col=viridis(20), main="Training", sub=dataset, xlab="Correlations")


perfect.pred.list = lapply( c(1,4,7,10,13), function(i) which(cor.all[f.zz,i]==1 & predcor.all[f.zz,i]==1 ) )
names(perfect.pred.list) = quads


venn( perfect.pred.list, zcol=candy_colors , cexil = 1, cexsn = 1.5, ellipse=T)
text(200,900,"Perfect predictors", cex=1.5)

perfect.pred = sapply( c(1,4,7,10,13), function(i) (cor.all[ ,i]==1 & predcor.all[ ,i]==1 ) )
colnames(perfect.pred) = quads 

for(j in 1:n_quads){ 
 temp4 =  X.rank2[f.zz,][which(perfect.pred[f.zz,j]==1),pData$Quad==j]
 heatmap.3(temp4, main=quads[j], col=cividis(5))
} 

pdf("heatmaps.modules.pdf")
for(j in 1:n_quads){ 
 temp4 =  cor(t(X.rank2[f.zz,][which(perfect.pred[f.zz,j]==1),pData$Quad==j]), m="s")
 heatmap.3(temp4, main=quads[j], col=magma(10))
} 
dev.off() 

```

```{r}
# Fig S2
pdf("perfect.pred.prop.pdf")
hist(log2(1+X.cpm.all[f.zz,]), 
     freq=F, 
     col=viridis(10)[5], 
     main="All expression", 
     xlab="log2 (1 + CPM)")
abline( v = mean(log2(1+X.cpm.all[f.zz,])) , col=2 , lwd=2   )

hist(log2(1+X.cpm.all[f.zz,][unique(unlist(perfect.pred.list)) ,]), 
     freq=F, 
     col=viridis(10)[3], 
     main="Perfect predictor expression", 
     xlab="log2 (1 + CPM)")
abline( v = mean(log2(1+X.cpm.all[f.zz,][unique(unlist(perfect.pred.list)) ,])) , col=2, lwd=2   )



par(mfrow=c(2,3))
for(i in 1:n_quads){
  
plot(log2(1+rowMeans(X.cpm.all[f.zz,pData$Quad==i][which(perfect.pred[f.zz,i]==1),])),
log2(1+rowMeans(X.cpm.all[f.zz,pData$Quad!=i][which(perfect.pred[f.zz,i]==1),])), pch=19, main=quads[i], 
xlab="Average log2 CPM + 1", ylab="Average log2 CPM+1 all other quads")
abline(0,1)  
} 



par(mfrow=c(2,3))
for(i in 1:n_quads){
  
  plot(log2(1+rowSD(X.cpm.all[f.zz,pData$Quad==i][which(perfect.pred[f.zz,i]==1),])),
       log2(1+rowSD(X.cpm.all[f.zz,pData$Quad!=i][which(perfect.pred[f.zz,i]==1),])), pch=19, main=quads[i], 
       xlab="SDs (log2+1)", ylab="SDs (og2+1) all other quads")
  abline(0,1)  
} 
  
  
  res.sds= sapply(1:n_quads, function(i) residuals( 
    rowSD(log2(1+X.cpm.all[f.zz,pData$Quad==i][which(perfect.pred[f.zz,i]==1),])),  
    rowSD(log2(1+X.cpm.all[f.zz,pData$Quad!=i][which(perfect.pred[f.zz,i]==1),])),-1,1,0) )
  
  res.sds2= sapply(1:n_quads, function(i) residuals( 
    rowSD(log2(1+X.cpm.all[f.zz,pData$Quad==i] )),  
    rowSD(log2(1+X.cpm.all[f.zz,pData$Quad!=i])),-1,1,0) )
  
  res.sds2= sapply(1:n_quads, function(i) residuals( 
    rowSD((X.cpm.all[f.zz,pData$Quad==i] )),  
    rowSD((X.cpm.all[f.zz,pData$Quad!=i])),-1,1,0) )
  


  
par(mfrow=c(2,3))
for(i in 1:n_quads) { 
hist(res.sds[[i]], main=quads[i], col=makeTransparent(candy_colors[i]), border=NA, xlab="Residuals of the SDs (log2 1+CPM)" )  
abline(v= mean( res.sds[[i]]), col=candy_colors[i], lwd=2) 
}

hist(unlist(res.sds) , main="", xlab="Residuals of the SDs (log2 + 1 CPM)", col=1, border=NA) 
abline(v= mean( unlist(res.sds)), col=2, lwd=2) 


dev.off()

```


```{r}
temp.clust = list() 
clust = list() 
for( j in 1:5){
  temp = X.rank2[f.zz,pData$Quad==j][which(perfect.pred[f.zz,j]),1:4]  
  h = hclust( dist(temp) )
  clust[[j]] = cutree ( h, h=1)
  temp.clust[[j]] = temp
}



perfect.pred.clust = list() 
for( j in 1:n_quads){
  perfect.pred.temp = matrix(0, ncol=length(unique(clust[[j]])), nrow=N ) 
  rownames(perfect.pred.temp) = rownames(X.rank2)
  for(i in unique(clust[[j]]) ){ 
    f.t = clust[[j]] == i  
    perfect.pred.temp[f.zz, i][which(perfect.pred[f.zz,j])][f.t] = 1
  }
  colnames(perfect.pred.temp) =  paste( quads[j] , unique(clust[[j]]) )
   perfect.pred.clust[[j]] = perfect.pred.temp
}

perfect.pred.clust = do.call(cbind, perfect.pred.clust)


perfect.pred.clust = cbind(perfect.pred, 1*(rowSums(perfect.pred, na.rm=T)>0), perfect.pred.clust)

colnames(perfect.pred.clust)[6] = "union"


```


```{r}
for(ji in which(colSums(perfect.pred.clust) > 5)){ 
  j = which( substr( colnames(perfect.pred.clust)[ji],1,4 ) == quads) 
  temp  = X.zscores[f.zz,pData$Quad==j][which(perfect.pred.clust[f.zz,ji]==1),]
  temp4 = X.rank2[f.zz,][which(perfect.pred.clust[f.zz,ji]==1),pData$Quad==j]

  
  X_c = 1:3
  
  Y_c = colMeans(temp)
  std_Y_c = colSD(temp)
  
  xlab="Time point"
  ylab = "Z-scores"
  ymin = -2
  ymax = 2
  xmin = 0.5
  xmax = 3.5
  ft = 1:3 * 4 - 3
  
  plot(X_c, Y_c[ft], ylim = c(ymin, ymax), xlim = c(xmin, xmax),
       lwd = 4, type = "l", col = 0, bty = "n", xlab = xlab,
       ylab = ylab, cex.lab = 1.5, cex.axis = 1.2, main= colnames(perfect.pred.clust)[ji])
  for( ti in 1:4){
    ft = 1:3 * 4 - (4-ti)
    polygon(c(X_c, rev(X_c)), c(Y_c[ft] -std_Y_c[ft], rev(Y_c[ft] + std_Y_c[ft])),
            col = makeTransparent(tropical_colors[ti]), border = NA)
    points(X_c, Y_c[ft], pch=19, col = tropical_colors[ti])
  }
  
}
 
```

# Predictive gene models 
## Perfect predictor model 1 
```{r, eval=FALSE}

N = dim(X.cpm.all)[1]

rand.all= list()
for (r in 1:1000){ 
  print(r)
  X.temp = lapply(1:n_qt, function(i) shuffle_rows(X.cpm.all[,(1:n_q) + n_q*(i-1) ]))
  X.sim = do.call(cbind, X.temp)
  X.rank = lapply(1:n_qt, function(i)   t(apply(X.sim[,(1:n_q) + n_q*(i-1) ],  1, rank )))   
  X.rank2 = do.call(cbind, X.rank)
  X.tt1 = X.rank2[,(r_samp)]
  X.tt2 = X.rank2[,(r_samp+n_samp)]
  X.tt3 = X.rank2[,(r_samp+(n_samp*2))]
  
  X.t1 = X.sim[,-(r_samp)]
  X.t2 = X.sim[,-(r_samp+n_samp)]
  X.t3 = X.sim[,-(r_samp+(n_samp*2))]
  
  X.t1.q = (X.t1[,r_samp] + X.t1[,r_samp+n_samp] )/2
  X.t2.q = (X.t2[,r_samp] + X.t2[,r_samp+n_samp] )/2
  X.t3.q = (X.t3[,r_samp] + X.t3[,r_samp+n_samp] )/2
  
  X.t1.r = lapply(1:n_quads, function(i) t(apply(X.t1.q[,(1:n_q) + n_q*(i-1) ],  1, rank) ))
  X.t2.r = lapply(1:n_quads, function(i) t(apply(X.t2.q[,(1:n_q) + n_q*(i-1) ],  1, rank) ))
  X.t3.r = lapply(1:n_quads, function(i) t(apply(X.t3.q[,(1:n_q) + n_q*(i-1) ],  1, rank) ))
  
  X.tr1 = do.call(cbind,X.t1.r )
  X.tr2 = do.call(cbind,X.t2.r )
  X.tr3 = do.call(cbind,X.t3.r )
  #toc() 
  
  rand.all[[r]] = sapply(1:n_quads, function(j) 
    sum(rowSums(  X.tr1[f.zz,((1:n_q)+n_q*(j-1))] == X.tt1[f.zz,((1:n_q)+n_q*(j-1))]  
                & X.tr2[f.zz,((1:n_q)+n_q*(j-1))] == X.tt2[f.zz,((1:n_q)+n_q*(j-1))]  
                & X.tr3[f.zz,((1:n_q)+n_q*(j-1))] == X.tt3[f.zz,((1:n_q)+n_q*(j-1))]) == n_q )
  )
} 
save(rand.all , file="null_perfect_pred.Rdata") 

hist(unlist(rand.all ) , freq=F, col="grey", border=NA, xlim=c(0,100))
# abline(v=mean(unlist(rand.all ) ), lwd=2, col=1)
abline(v= colSums(perfect.pred[f.zz,], na.rm=T), col=original_colors[1:5], lwd=2)

perfect.pred.ns = colSums(perfect.pred[f.zz,]  )
rand.pred  = t(do.call (cbind, rand.all ) )

perfect.pred.ns.mean =  mean(perfect.pred.ns)
rand.pred.mean  = rowMeans(rand.pred) 

hist( rand.pred.mean , freq=F, col="grey", border=NA, xlim=c(0,100))
abline(v=mean(rand.pred.mean  ), lwd=2, col=1)
abline(v= perfect.pred.ns.mean , col=4, lwd=2)


```
## Perfect predictors model 2 (null)

```{r}

preds.cds = list()
maxns = sort(c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,193,200,300,400,500,600,700,800,900,1000)) 
 
for( maxn in maxns ){
  print(maxn)
  # Adding predictors
  pred.scores.list = list()
  n = floor(maxn/5) 
  for(r in 1:100){
    X.tt1 = X.tt1.orig
    X.tr1 = X.tr1.orig
    X.tt2 = X.tt2.orig
    X.tr2 = X.tr2.orig
    X.tt3 = X.tt3.orig
    X.tr3 = X.tr3.orig
    # Testing 
    for(j in 1:5){
      # j = 1   
      ki = 1:4 + 4 *(j-1)
      s1 = sample(sum(f.zz),n)
      s2 = s1
      s3 = s1
      # s2 = sample(sum(f.zz),n)
      #  s3 = sample(sum(f.zz),n)
      
      for(i in ki){
        X.tt1[f.zz,][s1,i] <- X.tr1[f.zz,][s1,i] 
        X.tt2[f.zz,][s2,i] <- X.tr2[f.zz,][s2,i]  
        X.tt3[f.zz,][s3,i] <- X.tr3[f.zz,][s3,i] 
      }
    }
    
    
    pred.temp = matrix(0, ncol=3, nrow=5)
    for(j in 1:5){
      inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
      ki = ((1:4)+4*(j-1))
      q1 = 1:sum(f.zz)
      q2 = q1
      q3 = q1 
      t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
      t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
      t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
      pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    }
    pred.scores.list[[r]] = pred.temp 
    
  }
  preds.cds[[maxn]] = sapply(1:100, function(i) mean(pred.scores.list[[i]]))
  
} 
save(preds.cds, maxns, file="model.preds.cds.Rdata" )

#load("model.preds.cds.Rdata")
xlab="Number of genes"
ylab="Average score"

Y_c = sapply(maxns, function(i) mean(preds.cds[[i]]) )
std_Y_c = sapply(maxns, function(i) sd(preds.cds[[i]]) )
X_c = maxns
#/5
ymin  = min( Y_c - std_Y_c)
ymax = max( Y_c + std_Y_c)


plot( X_c, Y_c, ylim = c(ymin, ymax),   lwd = 4, type = "l", col = 0, bty = "n", xlab = xlab, ylab = ylab )
polygon(c(X_c, rev(X_c)), c(Y_c - std_Y_c, rev(Y_c + std_Y_c)), col = "lightgrey", border = NA)
points(X_c, Y_c, ylim = c(ymin, ymax) , pch=19, cex=0.5)


```

## Perfect predictors model 3 (distributed)
```{r, eval = FALSE}
# Made by models.R
load("../data/bulk/distr.model.comb.Rdata")
load("../data/bulk/model.means.shuff.Rdata")


means.propv = sapply(1:length(fracs), function(fi) mean( unlist(rand.propv[[fi]]),na.rm=T))
means.fcs = sapply(1:length(fracs), function(fi) mean( unlist(rand.dist[[fi]]),na.rm=T))

 
i = 1 
propve.fcs = list()
for( propi in propve){ 
  mini = which.min(means.propv < propi ) 
 
  propve.fcs[[i]] = means.fcs[mini]
  i =i +1 
} 

  
propve2.fcs = get_value(propve2[3], propve[19], propve.fcs[[3]] , propve.fcs[[19]], propve) 


library(akima)
x = log10(maxns)
y = propve2.fcs[1:51]
xa = sort(rep(x, length(y)))
ya = rep(y, length(x))
ya[is.na(ya)] = 0
za = array(cd.dist.mat )
temp2 = cbind(xa,ya,za)[!is.na(za),]
zi = interp(temp2[,1],temp2[,2], temp2[,3])



conts = c(1.5, 1.8,2) 
x = list()
yo = list() 
for(i in 1:length(conts) ) {
  cont = which(round(zi$z, 1) == conts[i], arr.ind=T)
  x[[i]] = zi$x[cont[,1]]
  y = zi$y[cont[,2]]
  o = order(x[[i]])
  x[[i]]= x[[i]][o]
  y = y[o]
  lo = loess(y~x[[i]])
  yo[[i]] = predict(lo)
}  
 

filled.contour( log10(maxns), propve2.fcs[1:51], ylim=c(0.16,0.5), t( cd.dist.mat), col= magma(400), levels=(0:400)/100, axes=F,
                xlab="Number of genes",
                ylab="|log2 FCs|", frame.plot=F,
                plot.axes = { axis(1, at=0:4, lab=10^(0:4)) ; axis(2);
                sapply(1:3, function(i) lines(x[[i]],yo[[i]], col="white", lwd=2));
                #contour(log10(maxns),  propve2.fcs[1:51] ,   t( cd.dist.mat), levels = c(1.5,1.8,2 ), add=T, col="white");    
                }) 


```


## B: X scaffold genes 
```{r}
load("../data/ase/Xscaffolds.Rdata")
 
 
pred.temp = matrix(0, ncol=3, nrow=5)
pval.temp = matrix(0, ncol=3, nrow=5)
pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
setsizes.temp = matrix(0, ncol=3, nrow=5) 
genesets.temp = list() 
scorematrix.temp = list() 

for(j in 1:5){
  inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
  ki = ((1:4)+4*(j-1))
  m = match( attr$scaffold, scaffoldsX.sub)
  m = match( attr$scaffold, scaffoldsX.prop3)
  f.nx = !is.na(m)
  
  q1 = which(f.nx[f.zz] & !is.na(cor.t2t3[[j]][f.zz]))
  q2 = which(f.nx[f.zz] & !is.na(cor.t1t3[[j]][f.zz])) 
  q3 = which(f.nx[f.zz] & !is.na(cor.t1t2[[j]][f.zz]))
  
  t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
  t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
  t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
  
  pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
  pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
  pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
  
  setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
  genesets.temp[[j]] = list( (q1),  (q2),  (q3))   
  scorematrix.temp[[j]] = list( t1[[1]], t2[[1]], t3[[1]] )   
  
}

 
barplot(colMeans(pred.temp), col=mm_colors); abline(h=1, col=4, lty=2)
barplot(rowMeans(pred.temp), col=candy_colors); abline(h=1, col=4, lty=2)

heatmap.3(pred.temp, 
          Rowv=F, 
          Colv=F, 
          col=magma(5), 
          ColSideColors=mm_colors[1:3], 
          RowSideColors=candy_colors[1:5], 
          labRow=quads, 
          labCol=tlab)

plot( pval.sum2[,1], x9r/sum(x9r), pch=19, xlab="Average score", ylab="Density", type="o")
abline(v=mean(pred.temp[2:4,]), col=4, lwd=2)
abline(v=1, col=2, lty=2 )

pval.temp = pval.sum2[match(mean(pred.temp[2:4,]), pval.sum2[,1]),2]

plot( pval.sum2, pch=19, xlab="Average score", ylab="Cummulative Density", type="o")
abline(v=mean(pred.temp[2:4,]), col=4, lty=2)
abline(h= pval.temp, col=4, lty=2)
abline(v=1, col=2, lty=2)




plot( pval.sum[,1], x15r/sum(x15r), pch=19, xlab="Average score", ylab="Density", type="o")
abline(v=mean(pred.temp), col=4, lwd=2)
abline(v=1, col=2, lty=2 )

pval.temp = pval.sum[match(mean(pred.temp ), pval.sum [,1]),2]

plot( pval.sum , pch=19, xlab="Average score", ylab="Cummulative Density", type="o")
abline(v=mean(pred.temp), col=4, lty=2)
abline(h= pval.temp, col=4, lty=2)
abline(v=1, col=2, lty=2)

 
```


## C: ASE feature genes 
```{r}
load("../data/ase/ase.feature.genes.Rdata")


pred.temp = matrix(0, ncol=3, nrow=5)
pval.temp = matrix(0, ncol=3, nrow=5)
pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
setsizes.temp = matrix(0, ncol=3, nrow=5) 
genesets.temp = list() 
scorematrix.temp = list() 

for(j in 1:5){
  inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
  ki = ((1:4)+4*(j-1))
  q1 = which(!is.na(cor.t2t3[[j]][f.zz]) & !is.na(match( attr$ensemblID[f.zz], feature.genes[[j]][[1]])))
  q2 = which(!is.na(cor.t1t3[[j]][f.zz]) & !is.na(match( attr$ensemblID[f.zz], feature.genes[[j]][[2]]))) 
  q3 = which(!is.na(cor.t1t2[[j]][f.zz]) & !is.na(match( attr$ensemblID[f.zz], feature.genes[[j]][[3]]))) 

  t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
  t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
  t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
  
  pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
  pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
  pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
  
  setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
  genesets.temp[[j]] = list( (q1),  (q2),  (q3))   
  scorematrix.temp[[j]] = list( t1[[1]], t2[[1]], t3[[1]] )   
  
}

 
barplot(colMeans(pred.temp), col=mm_colors); abline(h=1, col=4, lty=2)
barplot(rowMeans(pred.temp), col=candy_colors); abline(h=1, col=4, lty=2)

heatmap.3(pred.temp, 
          Rowv=F, 
          Colv=F, 
          col=magma(5), 
          ColSideColors=mm_colors[1:3], 
          RowSideColors=candy_colors[1:5], 
          labRow=quads, 
          labCol=tlab)
 


plot( pval.sum[,1], x15r/sum(x15r), pch=19, xlab="Average score", ylab="Density", type="o")
abline(v=mean(pred.temp), col=4, lwd=2)
abline(v=1, col=2, lty=2 )

pval.temp = pval.sum[match(mean(pred.temp), pval.sum[,1]),2]

plot( pval.sum, pch=19, xlab="Average score", ylab="Cummulative Density", type="o")
abline(v=mean(pred.temp), col=4, lty=2)
abline(h=pval.temp, col=4, lty=2)
abline(v=1, col=2, lty=2)
```
 

# Coexpression analysis
## Make armadillo aggregate network
```{r, eval=FALSE}
require(matrixStats)
exprs = X.cpm.all
genes =rownames( X.cpm.all)
keep = (rowSums(exprs) > 0) & (rowSds(exprs) > 0)


for( i in 1:n_quads){
  
   for (k in 1:n_times) {
     
      j=(k-1)*n_q;
      f.c = which(pData$Quad==i)[(1:n_q)+j]

      sub = exprs[keep,f.c]
      net = build_coexp_network(sub, genes[keep])
      med = median(net, na.rm=T)
      net[is.na(net)] = med
      
      # Save individual time point and quad ranked networks 
      save(net, file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      
      if(i==1 & k ==1 ) {
        agg = net
        } else {
          agg = net + agg
        }
      
   }
}
net = NULL
gc()
save(agg, file="coexp.agg.Rdata")

agg.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
agg = NULL
gc()
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)
save(agg.rank, file="coexp.agg.rank.Rdata")


```


## Make quad aggregate networks
```{r, eval=FALSE}

for( i in 1:n_quads){
  print(i) 
  sub = sub * 0 
  agg = sub * 0   
  for (k in 1:3) {
     print(k) 
      load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      if (k == 1) agg = rowSums( agg, na.rm=T) 
      else agg = rowSums( cbind(agg, net), na.rm=T) 
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, file=paste("coexp.agg.quad",i, "Rdata", sep=".") )
}   

```



## Make leave-one-quad-out aggregate networks and evaluate
```{r, eval=FALSE}
bottom = row(net) > col(net)

for( j in 1:n_quads){
  print(j) 
  sub = sub * 0 
  agg = sub * 0 
  loocv.quad = (1:n_quads)[-j]
  for( i in loocv.quad){
    print(i) 
     for (k in 1:n_times) {
       print(k) 
        load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
        agg = rowSums( cbind(agg, net), na.rm=T) 
     }
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, bottom,genes, f.zz, file=paste("coexp.quad.loo",j,"Rdata", sep=".") )  
}  
  

for( j in 1:n_quads){
  print(j) 
  sub = sub * 0 
  load(file=paste("coexp.quad.loo",j,"Rdata", sep=".") )  

  sub[bottom] = agg 
  sub = sub + t(sub)
  diag(sub) =  1 
  print("aurocs") 
  aurocs = run_GBA(sub, perfect.pred.clust, min=0)
  aurocs[[2]] = ""
  save(aurocs, file=paste("rocs.perfect.pred.clust.coexp.quad.loo",j,"Rdata", sep=".") )
}
```


## Make quad aggregate networks and evaluate
```{r, eval=FALSE}

for( i in 1:n_quads){
  print(i) 
  sub = sub * 0 
  agg = net * 0 
  for (k in 1:n_times) {
     print(k) 
      load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      agg = rowSums( cbind(agg, net), na.rm=T) 
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, file=paste("coexp.agg.quad",i, "Rdata", sep=".") )
}   

for( i in 1:n_quads){
  print(i) 
  load(file=paste("coexp.agg.quad",i,"Rdata", sep=".") )
  sub = sub * 0   
  sub[bottom] = agg 
  sub = sub + t(sub)
  diag(sub) =  1 
  print("aurocs") 
  aurocs = run_GBA(sub, perfect.pred.clust, min=0)
  aurocs[[2]] = ""
  save(aurocs, file=paste("rocs.perfect.pred.clust.coexp.quad",i,"Rdata", sep=".") )
}

```

```{r}
aurocs.loo = list() 
summary.aurocs.clust.loo = list()

for(j in 1:n_quads){
load(file=paste("coexp/rocs.perfect.pred.clust.coexp.quad.loo",j,"Rdata", sep=".") ) 
  aurocs.loo[[j]] = aurocs[[1]][1:6,]
  summary.aurocs.clust.loo[[j]] = aurocs[[1]][-(1:6),]
}

aurocs.circ = list() 
summary.aurocs.clust.circ = list()

for(j in 1:n_quads){
load(file=paste("coexp/rocs.perfect.pred.clust.coexp.quad",j,"Rdata", sep=".") ) 
  aurocs.circ[[j]] = aurocs[[1]][1:6,]
  summary.aurocs.clust.circ[[j]] = aurocs[[1]][-(1:6),]
}

 

summary.aurocs.clust.circ.mat = summary.aurocs.clust.circ[[j]][,1:3]  * 0 
summary.aurocs.clust.loo.mat = summary.aurocs.clust.loo[[j]][,1:3]   * 0 
quadids = matrix(unlist(strsplit(rownames( summary.aurocs.clust.circ.mat  ), " " )  )  , ncol=2, by=T)[,1]

for(j in 1:n_quads){
  id = quads[j]
  summary.aurocs.clust.circ.mat[quadids==id,] = summary.aurocs.clust.circ[[j]][quadids==id, ]
  summary.aurocs.clust.loo.mat[quadids==id,] = summary.aurocs.clust.loo[[j]][quadids==id, ]
} 

filt  = (colSums(perfect.pred.clust)[-c(1:6)] > 1 )

aurocs.loo.list = lapply(1:n_quads, function(j) summary.aurocs.clust.loo.mat[quadids==quads[j] & filt,1] )
aurocs.dens = lapply(n_quads:1, function(i) density(aurocs.loo.list[[i]],   from=0, to=1)  )   

o = order(sapply(1:n_quads, function(i) mean( aurocs.loo.list[[i]] )), decreasing = F) 
aurocs.dens.sums = sapply(o, function(i) aurocs.dens[[i]]$y )
plot(aurocs.dens[[1]]$x, aurocs.dens[[1]]$y, xlim=c(0,1), ylim=c(0,10), col=0, xlab="AUROCs", ylab="Density") 
for(i in n_quads:2){ 
  polygon( c(0,aurocs.dens[[1]]$x,1), c(0,rowSums(aurocs.dens.sums[,1:i]),0) , col=(candy_colors[o])[i], border=NA)
} 
polygon( c(0,aurocs.dens[[1]]$x,1), c(0,aurocs.dens.sums[,1],0) , col=(candy_colors[o])[1], border=NA)
abline(v= mean( unlist(aurocs.loo.list[[i]] )), col=1, lwd=2)

```



## GO for armadillos
```{r}
load("../data/genome/GO.human.Rdata")
gosums = colSums(GO.human.nonIEA )
annot =  GO.human.nonIEA 

load("coexp.agg.Rdata")
library(EGAD)

# Only use homologs
f.n = attr.arm$name[f.zz] != ""
network.arm =  sub[f.n,f.n]
bottom = row(network.arm) > col(network.arm)
sub = network.arm[bottom]
agg = rank(sub, na.last = "keep", ties.method = "average")
agg = agg/max(agg, na.rm=T)
sub = network.arm * 0
sub[bottom] = agg
sub = sub + t(sub)
diag(sub) = 1
rownames(sub) = attr.arm$name[f.zz]
colnames(sub) = attr.arm$name[f.zz]

aurocs.go1 = run_GBA( sub, annot,min=0 ,max=100 )
aurocs.go1[[2]] <- NULL 
save(aurocs.go1, file="aurocs.go1.arm.clust.Rdata")

```


## Human co-expression evaluation
```{r, eval=FALSE} 

load("gene_annotations_v29.Rdata")
attr.human = attr
load("gene_annotations_v0.95_mod.Rdata")
attr.arm = attr
### this file is too large, need to upload somewhere else
netfile="blood.rerank.Rdata"
label="blood.rerank"
    nettype= label
load(netfile)

    network = diag(length(genes.t))
    bottom = row(network) > col(network)
    colnames(network) = genes.t
    rownames(network) = genes.t
    network[bottom] = temp
    network = network + t(network)
    diag(network) = 1


m = match( rownames(network), attr.human$entrezID )
f.n = !is.na(m)
f.ah = m[f.n]

network = network[f.n,f.n]
rownames(network) = attr.human$name[f.ah]
colnames(network) = attr.human$name[f.ah]


m2 = match( attr.human$name[f.ah], attr.arm$name[f.zz] ) 
f.n2 = !is.na(m2) 
f.aa2 = m2[f.n2]

network.arm =  network[f.n2,f.n2]

perfect.pred.human.clust = perfect.pred.clust[f.zz,][f.aa2,]
aurocs = run_GBA( network.arm, perfect.pred.human.clust,min=0 )




bottom = row(network.arm) > col(network.arm)
sub = network.arm[bottom]
agg = rank(sub, na.last = "keep", ties.method = "average")
agg = agg/max(agg, na.rm=T)
sub = network.arm * 0
sub[bottom] = agg 
sub = sub + t(sub) 
diag(sub) = 1 

aurocs = run_GBA(sub, perfect.pred.human.clust,min=0 )
aurocs[[2]] <- NULL
save(aurocs, perfect.pred.human.clust, file="perfect.pred.clust.human.bloodagg.reranked.Rdata")


annot = GO.human.nonIEA[f.g,]
rownames(annot) = attr.human$name[f.a] 


aurocs.go1 = run_GBA( sub, annot,min=0 ,max=100 )
aurocs.go1[[2]] <- NULL 

```

# Cross-species co-expression 
```{r}
hb = hist( unlist(aurocs.loo.list), breaks=10, freq=F) 
aurocs.hist = sapply(5:1, function(i) hist(aurocs.loo.list[[i]], breaks=hb$breaks, freq=F)$density )   
barplot(t(aurocs.hist/10), col=rev(candy_colors[1:5] ) , ylim=c(0,1.5), border=NA, space=0, names=hb$breaks[-1]) 



load("perfect.pred.clust.human.bloodagg.reranked.Rdata")
filt = (colSums(perfect.pred.human.clust)[colSums(perfect.pred.human.clust) > 0] > 1  )
aurocs.loo.list$human = aurocs[[1]][filt,1]


aurocs.hum = sapply(5:1, function(i) hist(aurocs.loo.list$human , breaks=hb$breaks, freq=F)$density )   


load("/null_modules.blood.Rdata")
hist(aurocs.null, breaks=hb$breaks, col="lightgrey", border=NA, freq=F, main="", xlab="AUROCs")
abline(v=mean(aurocs.loo.list$human), lwd=2, col=2 ) 
abline(v=mean(aurocs.null), lty=2, lwd=2 ) 


bottom = row(sub.arm) > col(sub.arm)
armrank = rank(array(sub.arm[bottom]))
humanrank = rank(array(sub.blood[bottom]))


```



# Identity features: GO
```{r}
load("../data/genome/GO.human.Rdata")

annot = GO.human.nonIEA[f.go,]
rownames(annot) = attr.human$name[f.a]

m = match( attr$name, rownames(annot) ) 
f.a = !is.na(m)
f.go = m[f.a]

genesetsgo = matrix(0, ncol=dim(annot)[2], nrow=length(attr$name))
for(i in 1:dim(annot)[2]){ genesetsgo[f.a,i] = annot[f.go,i]}
colnames(genesetsgo) = colnames(annot)
genesetsgo = genesetsgo[,(colSums(genesetsgo[f.zz,]) >= 5 )]
genesetnamesgo = colnames(genesetsgo)


xlab="GO set size"
nr = 1:dim(genesetsgo)[2]
nr2 = colSums(genesetsgo[f.zz,]) 


pred.scores.list = list() 
pval.scores.list = list() 
pred.indscores.list = list()
setsizes.scores.list  = list() 
genesets.scores.list = list() 

i = 1 
for(k in nr){ 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  genesets.temp = list() 
  
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    kj = k*100
    
    q1 = genesetsgo[f.zz,k] == 1  
    q2 = q1
    q3 = q1
    
    t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
    t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
    t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
    genesets.temp[[j]] = list((q1), (q2), (q3))   
  }
  pred.scores.list[[i]] = pred.temp
  pval.scores.list[[i]] = pval.temp
  pred.indscores.list[[i]] = pred.ind.temp
  setsizes.scores.list[[i]] = setsizes.temp
  genesets.scores.list[[i]] = genesets.temp
  
  i = i + 1
}




cd = sapply(nr, function(i) mean(pred.scores.list[[i]]))
bc = sapply(nr, function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(nr, function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(nr, function(i) (pred.scores.list[[i]]))
de.pval = sapply(nr, function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]


save(
  genesetsgo, pred.scores.list,
  pval.scores.list,
  pred.indscores.list,
  setsizes.scores.list, 
  genesets.scores.list,nr, nr2, 
  cd,ab,de,bc, cd.pval,ab.pval,de.pval,bc.pval, 
  file="task1_results_GOarm.Rdata")

#load("task1_results_GOarm.Rdata")
results.mat = cbind(nr2, cd,cd.pval, t(ab),t(ab.pval),  t(bc), t(bc.pval),  t(de), t(de.pval))
save(results.mat,file="task1_results_GOarm.mat.Rdata")

```

```{r}
load("task1_results_GOarm.mat.Rdata")
filt1 = (results.mat[,1] >= 20  & results.mat[,1] <= 1000)
heatmap.3( -log10(results.mat[filt1,9:13]),Colv=F, col=cols5)
heatmap.3( (results.mat[filt1,4:8]),Colv=F, col=magma(20))

filt1_sig = rowSums((results.mat[filt1,4:8]) >= 2.3 ) > 0 
 
heatmap.3( (results.mat[filt1,4:8][filt1_sig,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
heatmap.3( -log10(results.mat[filt1,9:13][filt1_sig,]),Colv=F, col=cols5, ColSideCol=candy_colors[1:5], labCol=quads)

heatmap.3( (results.mat[filt1,4:8][ ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
heatmap.3( -log10(results.mat[filt1,9:13][ ,]),Colv=F, col=cols5, ColSideCol=candy_colors[1:5], labCol=quads)

filt1_sig = rowSums((results.mat[filt1,4:8]) >= 4 ) > 0 
temp = GO.voc[rownames(results.mat[filt1,][filt1_sig,]),]

heatmap.3( (results.mat[filt1,4:8][filt1_sig,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads, labRow = paste( rownames(results.mat[filt1,][filt1_sig,]),  temp[,2])  )



```

```{r}
load("aurocs.go1.blood.clust.Rdata") 
aurocs.go1.blood = aurocs.go1 
load("aurocs.go1.arm.clust.Rdata") 
aurocs.go1.arm = aurocs.go1 

m = match(rownames(aurocs.go1.blood[[1]]), rownames(aurocs.go1.arm[[1]] ) )
f.h = !is.na(m) 
f.a = m[f.h]

load("annot.agg.Rdata") 
mn = match( names(aurocs.go1.blood[[1]][f.h,1]), names(gosums )) 


main="GO"
xlab="GO AUROCs - human"
ylab="GO AUROCs - armadillo"


fgo = gosums[mn] >= 20 & gosums[mn] <= 1000

filt1_sig = rowSums((results.mat[filt1,4:8]) >= 4 ) > 0 
temp = GO.voc[rownames(results.mat[filt1,][filt1_sig,]),]


x =  aurocs.go1.blood[[1]][f.h,1][fgo]
y =  aurocs.go1.arm[[1]][f.a,1][fgo] 
sub = names(which(filt1_sig)) 
modcols = 0*(1:length(sub) )
for(i in 1:n_quads){ 
  modcols[ results.mat[filt1,i+3][filt1_sig] == 4  ]  = candy_colors[i]
} 

names(modcols) = names(which(filt1_sig))


plot(x,y, pch=19, col= "lightgrey", xlab=xlab, ylab=ylab) 
abline(0,1, lwd=2)
points(x[sub], y[sub], col=modcols, pch=19)
abline(v=mean(x[sub]), col=2, lwd=2 )
abline(h=mean(y[sub]), col=2, lwd=2 )
abline(v=mean(x), col="lightgrey", lty=2 )
abline(h=mean(y), col="lightgrey", lty=2 )

```



# Identity features: MsigDB
```{r}

load(file="../data/genome/msign.Rdata" )

xlab="MsigDB set size"
nr = 1:dim(genesetsgo)[2]
nr2 = colSums(genesetsgo[f.zz,]) 

pred.scores.list = list() 
pval.scores.list = list() 
pred.indscores.list = list()
setsizes.scores.list  = list() 
genesets.scores.list = list() 

i = 1 
for(k in which(nr2>5) ){ 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  genesets.temp = list() 
  
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    kj = k*100
    
    q1 = genesetsgo[f.zz,k] == 1  
    q2 = q1
    q3 = q1
    
    t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
    t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
    t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
    genesets.temp[[j]] = list((q1), (q2), (q3))   
  }
  pred.scores.list[[i]] = pred.temp
  pval.scores.list[[i]] = pval.temp
  pred.indscores.list[[i]] = pred.ind.temp
  setsizes.scores.list[[i]] = setsizes.temp
  genesets.scores.list[[i]] = genesets.temp
  
  i = i + 1
}

index = length(which(nr2>5))

cd = sapply(1:index, function(i) mean(pred.scores.list[[i]]))
bc = sapply(1:index, function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(1:index, function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(1:index, function(i) (pred.scores.list[[i]]))
de.pval = sapply(1:index, function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]



save(
  genesetsgo, pred.scores.list,
  pval.scores.list,
  pred.indscores.list,
  setsizes.scores.list, 
  genesets.scores.list,nr, nr2, 
  file="task1_results_msig.Rdata")



#load("task1_results_msig.Rdata")
nr22 = nr2[nr2>5]
cd = sapply(1:length(nr22), function(i) mean(pred.scores.list[[i]]))
bc = sapply(1:length(nr22), function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(1:length(nr22), function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(1:length(nr22), function(i) (pred.scores.list[[i]]))
de.pval = sapply(1:length(nr22), function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]
results.mat = cbind(nr22, cd,cd.pval, t(ab),t(ab.pval),  t(bc), t(bc.pval),  t(de), t(de.pval))
save(results.mat,file="task1_results_msig.mat.Rdata")


```




```{r} 
load(file="task1_results_msig.mat.Rdata")
filtk = grep("KEGG", rownames(results.mat) )
filtr = grep("REACTOME", rownames(results.mat) )
filtbc = grep("BIOCARTA", rownames(results.mat) )
filthlm = grep("HALLMARK", rownames(results.mat) )

filt_sig = rowSums((results.mat[filtk,4:8]) >= 4 ) > 0 
filt_sig = rowSums((results.mat[filthlm,4:8]) >= 3 ) > 0 
filt_sig = rowSums((results.mat[filtr,4:8]) >= 4 ) > 0 
filt_sig = rowSums((results.mat[filtbc,4:8]) >= 4 ) > 0 


heatmap.3( (results.mat[filtk,4:8][filt_sig ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
heatmap.3( -log10(results.mat[filtk,9:13][ ,]),Colv=F, col=cols5, ColSideCol=candy_colors[1:5], labCol=quads)

heatmap.3( (results.mat[filthlm,4:8][  ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
 
```

