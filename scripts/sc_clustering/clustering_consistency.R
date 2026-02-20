require(parallel)
require(ggplot2)
require(circlize)

DIR <<- "../bulk_transcriptome/data/"
ASE = FALSE
if (ASE) {
    load(paste0(DIR, "exprs_all.Rdata"))
    load(paste0(DIR, "metadata.Rdata"))
    load(paste0(DIR, "ase_ratios.test_train.Rdata"))
    load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
    source(paste0(DIR, "useful.r"))
} else {
    load(paste0(DIR, "cpm_strand_comb.Rdata"))
    load(paste0(DIR, "counts_strand_comb.Rdata"))
    load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
    load(paste0(DIR, "armadillo_hem.Rdata"))
    load(paste0(DIR, "ref.strand.test_train.Rdata"))
}
source("../scripts/sc_clustering/multi_sc_base.R")
source("../scripts/annotation/gene_enrichment.R")
# source("../scripts/multi_omics/sc_pseudo_profile.R")
source("../scripts/multi_omics/bulk_pred_utils.R")
source("../scripts/multi_omics/sc_pred_utils.R")
load(paste0(DIR, "Xscaffolds.Rdata"))
scaffoldsX = scaffoldsX.prop3
scaffoldsX = scaffoldsX.sub


quad_names = get.quad.names(X.cpm.all)


add.quad.names <- function(odf, sc_flag) {
    ann = t(apply(odf, c(1), function(x) {
        quad1 = substr(x[['group1_name']], 1, 4) 
        ind1 = substr(x[['group1_name']], 1, 6)
        time1 = substr(x[['group1_name']], 8, nchar(x[['group1_name']]))
        quad2 = substr(x[['group2_name']], 1, 4) 
        ind2 = substr(x[['group2_name']], 1, 6)
        time2 = substr(x[['group2_name']], 8, nchar(x[['group2_name']]))
        return(c(quad1, ind1, ind1==ind2, quad1==quad2, quad1 != quad2))
    }))
    colnames(ann) = c('quad', 'ind', 'time', 'within', 'across')
    return(cbind(odf, ann))
}

filtered.mat <- function(mat, pattern) {
    cols = grep(pattern, colnames(mat))
    temp = mat[,cols]
    filtered = temp[!is.na(rowSums(temp, na.rm=TRUE)),]
    filtered[is.na(filtered)] = 0
    return(filtered)
}

get.clust <- function(mat) {
    hc <- fastcluster::hclust(as.dist(1-cor(mat, method="spearman")), method="complete")
    return(hc)
}

extract.heatmap.cluster <- function(mat, plot_output=NULL, each=TRUE) {
    require(mclust)
    require(dplyr)
    require(fastcluster)
    require(ggdendro)
    cluster_size = 15
    results <- construct.annot.data(colnames(mat), row=TRUE)[[2]]
    time_cluster = paste0(results[,1], '_', results[,3])
    ind_cluster = paste0(results[,1], '_', results[,2])
    if (each) {
        df = NULL
        for (quad in unique(results[,1])) {
            cols = rownames(results[results[,1] == quad,])
            print(quad)
            print(cols)
            hc <- get.clust(mat[,cols])
            # HM = Heatmap(mat[,cols], cluster_columns=hc, km=15)
            vec <- cutree(hc, k=4)
            if (!is.null(plot_output)) {
                pdf(paste0(plot_output, '_', quad, '.pdf'), useDingbats=FALSE)
                plot(ggdendrogram(hc))
                dev.off()
            }
            time = adjustedRandIndex(vec, time_cluster[results[,1] == quad])
            print(time)
            ind = adjustedRandIndex(vec, ind_cluster[results[,1] == quad])
            print(ind)
            df = rbind(df, c(time, 'Time', quad))
            df = rbind(df, c(ind, 'ID', quad))
        }
        colnames(df) = c('ARI', 'cluster', 'quad')
        return(df)
    } else {
        hc <- get.clust(mat)
        HM = Heatmap(mat, cluster_columns=hc, km=15)
        vec <- cutree(hc, k=15)
        if (!is.null(plot_output)) {
            pdf(plot_output, useDingbats=FALSE)
            plot(ggdendrogram(hc))
            dev.off()
        }
        print(plot_output)
        print(vec)
        time = adjustedRandIndex(vec, time_cluster)
        print(time)
        ind = adjustedRandIndex(vec, ind_cluster)
        print(ind)
        return(c(time, ind))
    }
}

check.clustering.consistency <- function() {
    require(ComplexHeatmap)
    require(viridis)
    bulk_mat <- data.frame(readRDS('../data/bulk/bulk_exp.rds'))
    bulk_mat <- change.colnames(bulk_mat)
    rank_test <- readRDS('../data/bulk/bulk_exp_rank_test.rds')
    each = TRUE
    for (tail in c('', '_filt')) {
        if (tail != '') {
            bulk_mat = bulk_mat[f.zz,]
            rank_test = rank_test[f.zz,]
        }
        if (each) {
            colnames(rank_test) = colnames(bulk_mat)
            quad_col = c('purple', 'red', 'orange', 'lightgreen', 'skyblue')
            cluster_result_bulk = extract.heatmap.cluster(bulk_mat, paste0('dendro_bulk', tail), each=each)
            cluster_result_bulk = cbind(cluster_result_bulk, exp=rep('raw', dim(cluster_result_bulk)[1]))
            cluster_result_rank = extract.heatmap.cluster(rank_test, paste0('dendro_rank', tail), each=each)
            cluster_result_rank = cbind(cluster_result_rank, exp=rep('rank', dim(cluster_result_rank)[1]))
            df = data.frame(rbind(cluster_result_bulk, cluster_result_rank))
            df[,'ARI'] = as.numeric(as.character(df[,'ARI']))
            print(head(df))
            pdf(paste0('bar_ari_', tail, '_sign_each.pdf'), useDingbats=FALSE)
            g <- ggplot(df, aes(x=cluster, y=ARI, fill=quad))+geom_bar(stat='identity', position=position_dodge())+facet_wrap(~exp, scales="free")+theme_classic()+scale_fill_manual(values=quad_col)
            plot(g)
            dev.off()
            pdf(paste0('bar_ari_', tail, '_each.pdf'), useDingbats=FALSE)
            # g <- ggplot(df, aes(x=cluster, y=ARI, color=quad))+geom_jitter(width=0.05)+facet_wrap(~exp, scales="free")+theme_classic()+scale_color_manual(values=quad_col)+
            #     stat_summary(fun.y=mean, geom="crossbar", shape="\U2014", color="black")+
            #     stat_summary(fun.y=function(d) {mean(d) + sd(d)*c(-1,1)}, geom="crossbar", width=0.5, shape="\U2014", colour="gray")
            g <- ggplot(df, aes(x=cluster, y=ARI, fill=cluster))+geom_boxplot()+geom_jitter(aes(color=quad), width=0.05)+facet_wrap(~exp, scales="free")+theme_classic()+scale_color_manual(values=quad_col)
            plot(g)
            dev.off()
            df = cbind(df, abs_ARI=abs(df[,'ARI']), color=unlist(sapply(df[,'ARI'], function(x) {
                if (x < 0) return('negative');
                return('positive')
            })))

        } else {
            cluster_result_bulk = extract.heatmap.cluster(bulk_mat, paste0('dendro_bulk', tail, '.pdf'))
            cluster_result_rank = extract.heatmap.cluster(rank_test, paste0('dendro_rank', tail, '.pdf'))
            print(cluster_result_bulk)
            print(cluster_result_rank)
            print(length(cluster_result_bulk))
            print(length(cluster_result_rank))
            df = data.frame(ARI=c(cluster_result_bulk, cluster_result_rank), exp=c('raw', 'raw', 'rank', 'rank'), cluster=rep(c('Time', 'ID'), 2))
            print(head(df))
            pdf(paste0('bar_ari_', tail, '.pdf'), useDingbats=FALSE)
            g <- ggplot(df, aes(x=cluster, y=ARI, fill=cluster))+geom_bar(stat='identity')+facet_wrap(~exp, scales="free")+theme_classic()
            plot(g)
            dev.off()
            df = cbind(df, abs_ARI=abs(df[,'ARI']), color=unlist(sapply(df[,'ARI'], function(x) {
                if (x < 0) return('negative');
                return('positive')
            })))
            pdf(paste0('bar_ari_', tail, '_sign.pdf'), useDingbats=FALSE)
            g <- ggplot(df, aes(x=cluster, y=abs_ARI, fill=color))+geom_bar(stat='identity')+facet_wrap(~exp, scales="free")+theme_classic()+scale_fill_manual(values=c('blue', 'red'))
            plot(g)
            dev.off()
        }
        # for (exp in c('raw', 'rank')) {
        #     pdf(paste0('bar_ari_', tail, '_', exp, '.pdf'), useDingbats=FALSE)
        #     g <- ggplot(df[df[,'exp'] == exp,], aes(x=cluster, y=ARI, fill=cluster))+geom_bar(stat='identity')+theme_classic()
        #     plot(g)
        #     dev.off()
        # }
    }
}

check.clustering.consistency()