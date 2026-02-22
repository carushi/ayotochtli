require(viridis)
require(ggplot2)
require(dplyr)
require(parallel)

set.seed(42)

require(parallel)

DIR <<- "../bulk_transcriptome/data/"
ASE = TRUE
if (ASE) {
    load(paste0(DIR, "exprs_all.Rdata"))
    load(paste0(DIR, "metadata.Rdata"))
    load(paste0(DIR, "ase_ratios.test_train.Rdata"))
    load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
} else {
    load(paste0(DIR, "cpm_strand_comb.Rdata"))
    load(paste0(DIR, "counts_strand_comb.Rdata"))
    load(paste0(DIR, "gene_annotations_v0.95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
    load(paste0(DIR, "armadillo_hem.Rdata"))
    load(paste0(DIR, "ref.strand.test_train.Rdata"))
}
source(paste0("../scripts/", "helper.r"))
source("../scripts/sc_clustering/multi_sc_base.R")
source("../scripts/annotation/gene_enrichment.R")
source("../scripts/multi_omics/bulk_pred_utils.R")
source("../scripts/multi_omics/sc_pred_utils.R")

load(paste0(DIR, "Xscaffolds.Rdata"))
scaffoldsX = scaffoldsX.prop3
scaffoldsX = scaffoldsX.sub

# Variables 
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


exprs.all.filt = lapply( 1:5, function(i) exprs.all[[i]][!is.na(exprs.all[[i]][,1]),]) 
ratios = list() 
for(j in 1:5){ 
  X.temp = (exprs.all.filt[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  ratios[[j]] = cbind((X.temp[(1:nj),1:4] / (X.temp[(1:nj)+nj,1:4] +X.temp[(1:nj),1:4] ) ),
                      (X.temp[(1:nj),5:8] / (X.temp[(1:nj)+nj,5:8] +X.temp[(1:nj),5:8] ) ),
                      (X.temp[(1:nj),9:12] / (X.temp[(1:nj)+nj,9:12] +X.temp[(1:nj),9:12] ) ))
  
}

all_ratios = list() 
for(j in 1:5){ 
  X.temp = (exprs.all[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  all_ratios[[j]] = cbind((X.temp[(1:nj),1:4] / (X.temp[(1:nj)+nj,1:4] +X.temp[(1:nj),1:4] ) ),
                      (X.temp[(1:nj),5:8] / (X.temp[(1:nj)+nj,5:8] +X.temp[(1:nj),5:8] ) ),
                      (X.temp[(1:nj),9:12] / (X.temp[(1:nj)+nj,9:12] +X.temp[(1:nj),9:12] ) ))
  
}
target_quad = 5
mat = exprs.all[[target_quad]]
for (gene in c('DNMT3A', 'DNMT3B', 'TET1', 'TET2', 'TET3', 'ASXL1', 'JAK2')) {
    index = which(mat[,'name'] == gene)
    print(index)
    if (length(index) == 0 || any(is.na(index))) next
    index = index[1:(length(index)/2)]
    all_ratios[[target_quad]][index,]
    results <- construct.annot.data(colnames(all_ratios[[target_quad]]), row=TRUE)
    require(reshape2)
    require(ggplot2)
    m <- cbind(position=exprs.all[[target_quad]][index,5], all_ratios[[target_quad]][index,])
    m <- merge(melt(m, id=c('position')), results[[2]], by.x='variable', by.y=0, all.x=T)
    m[,'position'] = factor(m[,'position'])
    temp <- m %>% group_by(position, ident) %>% summarise(value = mean(value))

    g <- ggplot(data.frame(temp), aes(x=position, y=value, color=ident))+geom_point()+scale_color_viridis(discrete=TRUE)+theme_classic()
    pdf(paste0('ASE_', gene, '.pdf'), useDingbats=FALSE)
    plot(g)
    dev.off()
    g <- ggplot(data.frame(m), aes(x=ident, y=value, color=ident))+geom_point()+scale_color_viridis(discrete=TRUE)+theme_classic()+facet_wrap(~ position)
    pdf(paste0('ASE_', gene, '_all.pdf'), useDingbats=FALSE)
    plot(g)
    dev.off()
}

