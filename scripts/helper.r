quiet <- suppressPackageStartupMessages 

## Libraries/packages 
if( ! quiet(require("gplots")))
  install.packages("gplots")

if( ! quiet(require("gtools")))
  install.packages("gtools")

if( !quiet(require("plyr")))
  install.packages("plyr")

if( !quiet(require("RColorBrewer")))
  install.packages("RColorBrewer")

if( !quiet(require("viridis")))
  install.packages("viridis")

if( !quiet(require("beanplot")))
  install.packages("beanplot")

if( !quiet(require("beeswarm")))
  install.packages("beeswarm")

if( !quiet(require("vioplot")))
  install.packages("vioplot")


if( !quiet(require("pheatmap")))
  install.packages("pheatmap")

if( !quiet(require("ape")))
  install.packages("ape")

if( !quiet(require("tidyverse")))
  install.packages("tidyverse")

if( !quiet(require("boot")))
  install.packages("boot")

require(EGAD)
require(venn)
require(corrplot)

## Colors
cols = colorpanel(16, "red", "blue")
cols2 = brewer.pal(8, "Spectral")
cols3= rainbow(30)
cols4 = colorpanel(63, "lightgrey", "blue", "darkblue")
cols5 = colorpanel(300, "lightgrey", "red", "darkred")
cols6 = colorpanel(100, "lightgrey", "red", "darkmagenta")
cols7 = c("seagreen", "black", "darkmagenta")
cols8 = viridis(10)
cols9 = colorpanel(100, "white", "red", "darkmagenta")
cols10 = colorpanel(100, "white", "blue", "darkcyan")
cols11 = colorpanel(100, "white", "orange", "deeppink4")
cols12 = magma(100)
cols13 = viridis(100)
col_map= cols13
cols14 = c("grey", colorpanel(100, "black", "darkmagenta", "magenta") )
cols15 = viridis(15)

candy_colors = c("red3","purple4","darkorange1","green3","#26C3F0","brown") 
sex_colors = c("red", "blue")



## Identity prediction functions  
prediction_gene_scores <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
  if( length(geneset ) < 2){
    temp12 = matrix(0, ncol=length(t1), nrow=length(t1))
    
    q.sums.cross =  list(temp12, 0, "")
    return(q.sums.cross)
  }
  temp12 = sapply(t1, function(i) colSums(training[geneset,i] == testing[geneset,t2] )  )
  id_score = count_identity_mod(temp12)
  q.sums.cross =  list(temp12, id_score[[1]], id_score[[2]]  )
  return(q.sums.cross)
}


count_identity_mod <- function( temp1) {
  score = 0
  id.list=list()
  colnames(temp1) = rownames(temp1)
  temp1[is.na(temp1)] = 0
  if(sum(temp1) == 0) {  return (list(0, "")) }
  
  shuff = sample( dim(temp1)[1])  
  temp1 = temp1[shuff,shuff]
  for(i in 1:(dim(temp1)[1]-1)){
    pred1 = t(which(temp1== temp1[which.max(temp1)], arr.ind=T))
    p1 = rownames(temp1)[pred1[1]]
    p2 = colnames(temp1)[pred1[2]]
    
    if( p1 == p2){
      score = score +1
      id.list[[i]] = p1
    }
    p1left = rownames(temp1)[-pred1[1]]
    p2left = colnames(temp1)[-pred1[2]]
    temp1 = temp1[-pred1[1], -pred1[2]]
    
  }
  if(p1left == p2left ) {
    score = score +1
    id.list[[i+1]] = p1left
  }
  id.list = unlist(id.list)
  return (list(score, id.list))
}

count_identities_analytic <- function(n){
  A = permutations(n=n,r=n, v=1:n) 
  scores = matrix(0, ncol=1,nrow=dim(A)[1])
  for( j in 1:dim(A)[1] ){ 
    pred1 = 1:n
    pred2 = A[j,]
    score = 0
    temp1 = diag(n) * 0 
    colnames(temp1) = 1:n
    rownames(temp1) = 1:n
    
    #    pred1 = rbind(pred1,pred2)
    for(i in 1:n){ 
      temp1[pred1[i],pred2[i] ] = 1 
    }
    
    
    for(i in 1:(n-1)) {
      pred1 = t(which(temp1== temp1[which.max(temp1)], arr.ind=T))
      
      # predictor
      p1 = rownames(temp1)[pred1[1]]
      # predictions 
      p2 = colnames(temp1)[pred1[2]]
      
      if( p1 == p2){
        score = score +1
      }
      p1left = rownames(temp1)[-pred1[1]]
      p2left = colnames(temp1)[-pred1[2]]
      temp1 = temp1[-pred1[1], -pred1[2]]
      
      
    }
    if(p1left == p2left ) { 
      score = score +1
    }  
    scores[j,1] = score
  }
  return(scores)
}



prediction_gene_scores_mod <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
  if ( length(geneset) < 1 ){ 
    temp12 = matrix(0, ncol=length(t1), nrow=length(t1))
    q.sums.cross =  list(temp12, 0, "")
    return(q.sums.cross)
  }
  if( length(geneset ) == 1){
   temp12 = sapply(t1, function(i) 1*(training[geneset,i] == testing[geneset,t2] )  )
    
  } else { 
  temp12 = sapply(t1, function(i) colSums(training[geneset,i] == testing[geneset,t2] )  )
 } 
  id_score = count_identity_mod(temp12)
  q.sums.cross =  list(temp12, id_score[[1]], id_score[[2]]  )
  return(q.sums.cross)
  
}

prediction_gene_scores_indiv <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
 temp = list() 
  for(j in geneset) {
   temp12 = sapply(t1, function(i) 1*(training[j,i] == testing[j,t2] )  )
   id_score = count_identity_mod(temp12)
   temp[[j]] = id_score[[1]]
  } 
   return(unlist(temp) ) 
}


convolve_x <- function(x,xc,n){
  xcc =  convolve(xc, rev(x), type="open")
  if (n ==1){
    return(xcc)
  } else {
    return(convolve_x(x, xcc, n-1) ) 
  }
}


## Folded and unfolded normal distributions 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

mle_folded <- function(){ 
  mus = seq(0.5,1, by = 0.001)
  sigmas = seq(0,0.5, by = 0.01)
  
  mles = sapply(mus, function(mu)
    sapply(sigmas, function(sigma)
      -sum( log(dnorm( x,   mean = mu, sd = sigma )
                + dnorm( x,   mean = 1-mu, sd = sigma ) )  )))
  mles[!is.finite(mles)] = NA
  coefs = which(mles==min(mles, na.rm=T), arr.ind=T)
  return (list( mus[coefs[2]] , sigmas[coefs[1]]))
} 



## Variance to N relationship 
foo <- function(i,j) { 
  return( 0.25/(sd(unfold(i[j]))^ 2))
} 


## Gene set enrichment 
gene_set_enrichment <- function(genes, genes.labels, voc){
  
  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p
  
  m = match ( genes, genes.names )
  filt.genes  = !is.na(m)
  filt.labels = m[filt.genes]
  
  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g
  
  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]
  
  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { 
    genes.counts.set = genes.labels[filt.labels,]  
    test =  cbind(   (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set) # check again 
  } else { 
    genes.counts.set = colSums(genes.labels[filt.labels,]) 
    test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  }             ## does weird things with 0 sets
  
  
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  sigs = pvals < ( 0.05/length(pvals) )
  pvals.adj = p.adjust( pvals, method="BH")
  
  results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj", "sig" )
  
  return (results)
  
}

## Visualisations 
heatmap.3 <- function(mat, ...){
  heatmap.2( mat, ..., density="none", trace="none")
}

vioplot.2 <- function(x,colin=cols3,plotpoints=FALSE,...){
  if( class(x) == "list") {
    temp = unlist(x)
    yrange = range(temp, na.rm=T)
    xmax = length(x) + 1
    plot(0,0, xlim=c(0,xmax), ylim=yrange, axes=F, col=0, ..., bty="n")
    axis(1, labels = names(x), at = 1:(xmax-1)  )
    for( i in 1:(xmax-1)) {
      vioplot(x[[i]][ !is.na(x[[i]])], col=colin[i], horizontal=FALSE, at=i, add=TRUE, rectCol="gray", axes=F)
      if(plotpoints==TRUE){  points( jitter( rep(i, length(x[[i]])),amount=0.1 ), x[[i]],   pch=19 , ... ) } 
    }
    axis(2)
  }
}



vioplot.3 <- function(x,colin=cols3,plotpoints=FALSE,xlab="", ylab="",...){
  if( class(x) == "list") {
    temp = unlist(x)
    yrange = range(temp, na.rm=T)
    xmax = length(x) + 1
    plot(0,0, ylim=c(0,xmax), xlim=yrange, axes=F, ylab=ylab, xlab=xlab, col=0, ..., bty="n")
    axis(2, labels = names(x), at = 1:(xmax-1) , las=2 , cex.axis=0.2)
    
    for( i in 1:(xmax-1)) {
      vioplot(x[[i]][ !is.na(x[[i]])], col=colin[i], horizontal=TRUE, at=i, add=TRUE, rectCol="gray", axes=F)
      if(plotpoints==TRUE){  points( x[[i]], jitter( rep(i, length(x[[i]])),amount=0.1 ),   pch=19 , ... ) } 
    }
    axis(1)
  }
}



makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}



## Random functions
rank_std <- function(x)  { r = rank( x, na.last="keep"); r/max(r, na.rm=T)  }
colSD <- function( data){ return(apply( data, 2, sd, na.rm=T))}
rowSD <- function( data){ return(apply( data, 1, sd, na.rm=T))}
colSE <- function( data){ return( apply( data, 2, sd, na.rm=T)/sqrt(dim(data)[2]))}
rowSE <- function( data){ return( apply( data, 1, sd, na.rm=T)/sqrt(dim(data)[1]))}
se    <- function(x){ sd(x,na.rm=T)/sqrt(length(!is.na(x))) }
rmse  <- function(error){sqrt(mean(error^2, na.rm=T) )}
mae   <- function(error){ mean(abs(error), na.rm=T)}
residuals <- function(x,y,A,B,C){ (A*x + B*y + C) }
residuals2 <- function(x,y,A,B,C) { (A*x + B*y + C)/sqrt(A^2+B^2) }
residuals3 <- function(x,y,A,B,C) { abs(A*x + B*y + C)/sqrt(A^2+B^2) }

geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

geo_sd <- function(data) {
  log_data <- log(data)
  gs <- exp(sd(log_data[is.finite(log_data)]))
  return(gs)
}

geo_se <- function(data) {
  gs <- geo_sd(data)
  log_data <- log(data)
  gse <- gs/sqrt(sum(is.finite(log_data)))
  return(gse)
}


lm.studentized  <- function(x,y){
  z = lm(y ~ x )
  z = rstudent(z)
  return( rank(abs(z)) )
}

lm.function  <- function(x,y){
  z = lm(y ~ x )
  return( rank(abs(z$residuals)) )
}


z_scores <- function(x) {
  mean_x = mean(x, na.rm=T)
  sd_x = sd(x, na.rm=T)
  z =  (x - mean_x) / (sd_x)
  return(z)
}

z_scores_mod <- function(x) {
  med_x = median(x, na.rm=T)
  mad_x = median(abs(x-med_x), na.rm=T)
  z =  0.6745 * (x - med_x) / (mad_x)
  return(z)
}

calc_cpm <-function(X){
  K  = colSums(X)
  X.cpm = sapply(1:length(K), function(k) 10^6*X[,k]/K[k] )
  return(X.cpm)
}


# Given x and two points
get_value <- function( x1, x2, y1,y2, x) {
  m = (y2 - y1) / (x2 - x1 )
  y = y1 + m *( x - x1)
  return(y)
}

# Given y and two points
get_value_x <- function( x1, x2, y1,y2, y) {
  m = (y2 - y1) / (x2 - x1 )
  x = x1 + (y - y1)/m
  return(x)
}


## Formats the density distribution from the histogram function
get_density <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$density, hist$density))
  y = c(0,y,0)
  
  return(cbind(x,y))
}


## Formats the counts distribution from the histogram function
get_counts <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$counts, hist$counts))
  y = c(0,y,0)
  
  return(cbind(x,y))
}

convolve_nets <- function(netA,netB,f){
  n <- order(netA)
  temp_netA <- netA[n]
  
  temp_netB = convolve( netB[n], rep(1,f),type="filter")
  convolved = cbind(temp_netA[(f/2):(length(temp_netA)-f/2)],temp_netB/f)
}

calc_fdrs_recur  <- function( data, pp = 0.05) {
  
  temp1 = lapply(1:1000, function(i) shuffle_cols(data) )
  temp2 = sapply(1:1000, function(i) rowSums(temp1[[i]], na.rm=T))
  
  nmax = dim(data)[2] + 1
  
  recur = rowSums(data, na.rm=T)
  ob =  count_recur( as.matrix(recur) , nmax)
  ex = rowSums(sapply(1:1000, function(j) count_recur( as.matrix(temp2[,j] ), nmax ) ) )/1000
  
  test = cbind( ob, ex)
  
  observed = cbind( (rev(cumsum(rev(test[,1])))), 0:(nmax-1))
  expected = cbind( (rev(cumsum(rev(test[,2])))), 0:(nmax-1))
  
  FDR = expected[,1]/observed[,1]
  Pt = expected[ min(which(FDR < pp ) ), 2]
  sig = sum(test[  (which(FDR < pp )) ,1], na.rm=T)
  res = list(FDR, test, Pt, sig, pp)
  names(res) = c("FDR", "test", "Pt", "sig", "pp")
  return( res )
}


shuffle_cols <- function( data ) {
  nc = dim(data)[2]
  nr = dim(data)[1]
  return (sapply(1:nc, function(i) data[sample(nr),i] ))
}



shuffle_rows <- function( data ) {
  nc = dim(data)[2]
  nr = dim(data)[1]
  return (t(sapply(1:nr, function(i) data[i,sample(nc)] )))
}


count_recur <- function(data,nmax){
  freq = plyr::count(data[,1])
  res = matrix(0, nrow=nmax, ncol=1 )
  rownames(res) = 0:(nmax-1)
  m = match(0:(nmax-1), freq[,1])
  f.r = !is.na(m)
  f.f = m[f.r]
  res[f.r,1] = freq[f.f,2]
  return(res)
}

if (!file.exists("armadillo.helper.Rdata")) {

    ## Convolved distributions for armadillos 
    A = count_identities_analytic(4)
    h = hist(A, breaks=c(-1:4))
    x = h$counts 

    x15r = convolve_x(x,x,14)
    x5r = convolve_x(x,x,4)
    x3r = convolve_x(x,x,2)
    x9r = convolve_x(x,x,8)

    # Analytic nulls
    pval.rsum = cbind((0:12)/3, rev(cumsum(rev(x3r/sum(x3r)))) ) 
    pval.csum = cbind((0:20)/5, rev(cumsum(rev(x5r/sum(x5r)))) ) 
    pval.sum = cbind((0:60)/15, rev(cumsum(rev(x15r/sum(x15r)))) ) 
    pval.sum2 = cbind((0:36)/9, rev(cumsum(rev(x9r/sum(x9r)))) ) 
    pval.raw = cbind(0:4, rev(cumsum(rev(x/sum(x)))) )
    pval.sum[pval.sum[,2] < 0.5e-15,2]  = 2e-16
    pval.sum2[pval.sum2[,2] < 0.5e-15,2]  = 2e-16
    save(pval.csum, pval.raw, pval.rsum, pval.sum, pval.sum2, x15r, x9r, x5r, x3r, x, convolve_x, quadlabels, timelabels,sexlabels, labels, quads, tlab, lane, sex, quad, pData, count_identities_analytic, count_identity_mod, prediction_gene_scores, original_colors, tropical_colors, mm_colors,  candy_colors, file="armadillo.helper.Rdata")
}