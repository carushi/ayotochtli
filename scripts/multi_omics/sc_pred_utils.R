
scData <<- read.table('../data/phenotype_data/scPdata.tsv', header=T, sep="\t")


cor.trio.sc <- function(training, test, ident=5, offset=0) {
    N = dim(training)[1]
    j = (ident-1)*n_q
    return(sapply(1:N, function(i) cor( training[i,(1:n_q)+j], test[i,(1:n_q)+j], m="s") ))
}

rank.correlation.sc.trio <- function(df, quad_names, target) {
    pred.t1t2 = cor.trio.sc(df[['training']][['X.tr3']], df[['test']][[target]])
    pred.t1t3 = cor.trio.sc(df[['training']][['X.tr2']], df[['test']][[target]])
    pred.t2t3 = cor.trio.sc(df[['training']][['X.tr1']], df[['test']][[target]])
    predcor.all = do.call(cbind, lapply(1:n_quads, function(j) cbind( pred.t2t3[[j]], pred.t1t3[[j]], pred.t1t2[[j]]))) 
    colnames(predcor.all) <- unlist(lapply(quad_names, function(x){return(c(paste0(x, '_', c('23_1', '13_2', '12_3'))))}))
    rownames(predcor.all) <- rownames(df[['training']][['X.tr1']])
    return(predcor.all)
}

visualize.single.cell.data <- function(df) {
    for (name in names(df[['test']])) {
        if (any(name == c('X.tt1', 'X.tt2', 'X.tt3', 'X.ttNA')))
            next
        print(c('cluster', name))
        plot.heatmap.clustering(df[['test']][[name]], paste0('sc_average_rank_', name), 'rank')
    }
}

get.batch.weight <- function(weight) {
    weight_list = list()
    for (batch in c(-1:3)) {
        if (batch >= 0) {
            bweight = weight[grepl(paste0(batch, '_'), names(weight))]
            names(bweight) = unlist(sapply(names(bweight), function(x) substr(x, 3, nchar(x))))
        } else { # batch == -1
            bweight = NULL
            for (i in 0:3) {
                tweight = weight[grepl(paste0(i, '_'), names(weight))]
                names(tweight) = names(weight[grepl(paste0(i, '_'), names(weight))])
                if (is.null(bweight)) {
                    labels = unique(unlist(sapply(names(weight), function(x) substr(x, 3, nchar(x)))))
                    bweight = rep(0, length(labels))
                    names(bweight) = labels
                }
                total = sum(tweight)
                for (label in names(bweight)) {
                    blabel = paste0(i, '_', label) 
                    if (any(blabel == names(tweight))) {
                        bweight[label] = bweight[label]+tweight[which(blabel == names(tweight))[1]]/total
                    }
                }
            }
            bweight = bweight/4
        }
        weight_list[[as.character(batch)]] = bweight
    }
    return(weight_list)
}

weighted.average <- function(data, bweight, batch, remove.undefined=FALSE) {
    mat = matrix(0, nrow=dim(data)[1], ncol=4)
    print(dim(bweight))
    for (i in 0:3) {
        mat[,i+1] = t(colSums(do.call(rbind, lapply(names(bweight), function(x) {
            if (any(colnames(data) == paste0('g', i, '-', x)))
                return(bweight[x]*data[,paste0('g', i, '-', x)])
            else
                return(0)
        }))))
    }
    rownames(mat) = rownames(data)
    colnames(mat) = 1:4
    stopifnot(dim(mat)[1] < 40000)
    return(mat)
}

average.expression <- function(sc, gene_set, label_list) {
    Idents(sc) <- label_list
    data = AverageExpression(sc,
        assays = 'RNA',
        features = gene_set,
        return.seurat = FALSE,
        group.by = "ident"
    )
    data <- data[[1]]
    saveRDS(data, file='each_sc_cell_type.rds')
    return(data)
}

weighted.expression <- function(data, sc, label, label_list) {
    results = list()
    batch_list = unique(sc@meta.data[,'batch'])
    weight = table(label_list)

    # each cell-type specific profile
    for (name in unique(sc@meta.data[,label])) {
        pseudo_mat = do.call(cbind, lapply(0:3, function(x) {
            return(data[,paste0('g', x, '-', name)])
        }))
        rownames(pseudo_mat) = rownames(sc[["RNA"]])
        results[[as.character(name)]] = pseudo_mat
    }
    # each individual-specific (cell-type specific pseudo-bulk profiles
    # balanced according to the cell type population from one individual
    weight_list <- get.batch.weight(weight)
    for (batch in unique(sc@meta.data[,'batch'])) {
        pseudo_mat = weighted.average(data, weight_list[[as.character(batch)]], batch)
        rownames(pseudo_mat) = rownames(sc[["RNA"]])
        results[[as.character(batch+1)]] = pseudo_mat
    }
    pseudo_mat = weighted.average(data, weight_list[['-1']], -1)
    rownames(pseudo_mat) = rownames(sc[["RNA"]])
    results[[as.character(5)]] = pseudo_mat
    return(results)
}

adjusted.pseudo.profile <- function(sc, gene_set=NULL, label='celltype', plot_flag=TRUE) {
    label_list = unlist(apply(sc@meta.data, c(1), function(x) {
        return(paste0(x[['batch']], '_', x[[label]]))
    }))
    data <- average.expression(sc, gene_set, label_list)
    results <- weighted.expression(as.matrix(data), sc, label, label_list)
    if (plot_flag) {
        plot.heatmap.clustering.annotation(data, paste0('sc_average_', 'all'), 'zscore')
    }
    return(results)
}

get.sc.rank.data <- function(data, rownames) {
    sc.rank = lapply(1:1, function(i)   t(apply(data[,(1:n_q) + n_q*(i-1) ],  1, rank )))[[1]]
    sc.rank = sc.rank[!grepl('gSpikein', rownames(sc.rank)),]
    colnames(sc.rank) <- paste0('16-90', 1:4, '_4')
    if (dim(sc.rank)[1] < 100) exit()
    index = unlist(sapply(rownames, function(x){return(which(x == rownames(sc.rank))[1])}))
    sc.rank <- sc.rank[index,]
    if (dim(sc.rank)[1] < 100) exit()
    rownames(sc.rank) = rownames
    return(sc.rank)
}

add.single.cell.data <- function(df, fname, norm='average', label='celltype', plot_flag=TRUE) {
    require(Seurat)
    sc = readRDS(fname)
    sc = UpdateSeuratObject(sc)
    Idents(sc) <- sc@meta.data$batch
    if (norm == 'average' || norm == 'weighted') {
        data = AverageExpression(sc,
            assays = 'RNA',
            features = NULL,
            return.seurat = FALSE,
            group.by = "batch",
        )
        data = do.call(cbind, data)
        saveRDS(data, file='each_sc_batch.rds')
        print('?????')
        df[['test']][['X.tt4']] = get.sc.rank.data(data, rownames(df[['training']][['X.tr1']]))
        print('??????')
    }
    if (norm == 'adjusted' || norm == 'weighted') {
        results = adjusted.pseudo.profile(sc, plot_flag=plot_flag, label=label)
        saveRDS(results, file='each_weighted_sc_batch.rds')
        print('???')
        if (plot_flag) {
            for (name in 1:5) {
                plot.heatmap.clustering(results[[as.character(name)]], paste0('sc_average_', name), 'zscore')
            }
        }
        for (batch in names(results)) {
            df[['test']][[paste0('X.tt', as.numeric(batch)+4)]] = get.sc.rank.data(results[[batch]], rownames(df[['training']][['X.tr1']]))
        }
        print('???')
        for (name in unique(sc@meta.data[,label])) {
            print(name)
            if (is.na(name)) next
            df[['test']][[paste0('X.tt_', name)]] = get.sc.rank.data(results[[name]], rownames(df[['training']][['X.tr1']]))
            stopifnot(!is.na(df[['test']][[paste0('X.tt_', name)]][1,1]))
        }
    }
    return(df)
}

prediction.sc.comb <- function(df, quad_names, target, ident=5, offset=0) {
    pred1 = prediction.score.each.comb(df[['training']][['X.tr1']], df[['test']][[target]], ident=ident, offset=offset, rowSum_flag=TRUE)
    pred2 = prediction.score.each.comb(df[['training']][['X.tr2']], df[['test']][[target]], ident=ident, offset=offset, rowSum_flag=TRUE)
    pred3 = prediction.score.each.comb(df[['training']][['X.tr3']], df[['test']][[target]], ident=ident, offset=offset, rowSum_flag=TRUE)
    mat = cbind(pred1[[1]], pred2[[1]], pred3[[1]])
    colnames(mat) = unlist(sapply(1:3, function(x) paste0(pData[pData$Quad == ident,]$altID, '_', x)))
    return(mat)
}

prediction.sc.each.comb <- function(df, quad_names, target, ident=5, offset=0) {
    stopifnot(length(ident) == 1)
    pred1 = prediction.score.each.comb(df[['training']][['X.tr1']], df[['test']][[target]], ident=ident, offset=offset)
    pred2 = prediction.score.each.comb(df[['training']][['X.tr2']], df[['test']][[target]], ident=ident, offset=offset)
    pred3 = prediction.score.each.comb(df[['training']][['X.tr3']], df[['test']][[target]], ident=ident, offset=offset)
    mat = cbind(pred1[[1]], pred2[[1]], pred3[[1]])
    print(dim(mat))
    print(as.vector(unlist(sapply(1:3, function(x) { paste0(scData[scData[,'Target'] == target, 'altID'], '_', 4, '_', x) }))))
    colnames(mat) = as.vector(unlist(sapply(1:3, function(x) { paste0(scData[scData[,'Target'] == target, 'altID'], '_', 4, '_', x) })))
    return(mat)
}

add.sc.tr <- function(X.cpm.all, df, q_samp, fname) {
    df = add.single.cell.data(df, fname, 'weighted', plot=FALSE)
    first_for_each = ((q_samp-1)*n_q+n_samp*(0:2)+1)
    all_columns = sapply(first_for_each, function(x){return(x+1:n_q-1)})
    X.rank = lapply(first_for_each, function(i)   t(apply(X.cpm.all[,i+1:n_q-1],  1, rank )))
    X.rank2 = do.call(cbind, X.rank)

    X.t1.q = (X.rank2[,1:n_q] + df[['test']][['X.tt4']])/2
    X.t2.q = (X.rank2[,1:n_q+n_q] + df[['test']][['X.tt4']])/2
    X.t3.q = (X.rank2[,1:n_q+n_q*2] + df[['test']][['X.tt4']])/2

    X.tr1s = t(apply(X.t1.q, c(1), rank))
    X.tr2s = t(apply(X.t2.q, c(1), rank))
    X.tr3s = t(apply(X.t3.q, c(1), rank))

    colnames(X.tr1s) = pData$altID[first_for_each[1]+1:n_q-1]
    colnames(X.tr2s) = pData$altID[first_for_each[2]+1:n_q-1]
    colnames(X.tr3s) = pData$altID[first_for_each[3]+1:n_q-1]

    rownames(X.tr1s) = rownames(X.cpm.all)
    rownames(X.tr2s) = rownames(X.cpm.all)
    rownames(X.tr3s) = rownames(X.cpm.all)

    df[['training']][['X.tr1s']] = X.tr1s
    df[['training']][['X.tr2s']] = X.tr2s
    df[['training']][['X.tr3s']] = X.tr3s

    return(df)
}

plot.heatmap.clustering.annotation <- function(mat, output, zscore='') {
    require(viridis)
    require(ComplexHeatmap)
    require(ggsci)
    require(circlize)
    celltype=unlist(sapply(colnames(mat), function(x) unlist(strsplit(x, '_'))[2] ))
    celltype_vec = pal_ucscgb()(length(unique(celltype)))
    names(celltype_vec) = unique(celltype)
    batch=unlist(sapply(colnames(mat), function(x) unlist(strsplit(x, '_'))[1] ))
    batch_vec = pal_npg("nrc")(4)
    names(batch_vec) = 0:3
    annot_mat <- data.frame(cluster=colnames(mat), celltype=celltype, batch=batch, row.names='cluster')
    col_fun = colorRamp2(seq(min(mat), 3, length = 100), viridis(100))
    print(head(annot_mat))
    ha = HeatmapAnnotation(df=annot_mat, col=list(
        batch=batch_vec, celltype=celltype_vec
    ))
    mat <- mat[rowSums(mat) > 0,]
    if (zscore == 'zscore')
        mat <- t(apply(mat,  1, scale))
    print(head(mat))
    colnames(mat) <- rownames(annot_mat)
    h <- Heatmap(mat, 
        col=col_fun,
        # col=viridis(100),
        name = 'Pseudo-bulk profiles',
        show_row_names = FALSE,
        top_annotation=ha,
        cluster_columns = TRUE, cluster_rows=TRUE)
    pdf(paste0(output, '.pdf'), width=10, height=8)
    draw(h)
    dev.off()
}

plot.heatmap.clustering <- function(mat, output, method='') {
    require(viridis)
    require(ComplexHeatmap)
    mat <- mat[rowSums(mat) > 0,]
    mat <- mat[unlist(sapply(1:dim(mat)[1], function(x) { return(min(mat[x,]) != max(mat[x,])) })),]
    if (method == 'rank') {
        index = unlist(apply(mat, c(1), function(x) x[1]*1000+x[2]*100+x[3]*10+x[4]*1))
        mat <- mat[order(index, decreasing=TRUE),]
        h <- Heatmap(mat, 
            col=viridis(100),
            name = 'Pseudo-bulk profiles',
            show_row_names = FALSE,
            cluster_columns = FALSE, cluster_rows=FALSE)
        pdf(paste0(output, '.pdf'))
        draw(h)
        dev.off()
    } else {
        if (method == "zscore") {
            mat <- t(apply(mat,  1, scale))
            h <- Heatmap(mat, 
                col=viridis(100),
                name = 'Pseudo-bulk profiles',
                show_row_names = FALSE,
                cluster_columns = FALSE, cluster_rows=TRUE)
            pdf(paste0(output, '.pdf'))
            draw(h)
            dev.off()
        }
    }
}
