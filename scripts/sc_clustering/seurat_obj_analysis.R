require(Seurat)
require(scales)
require(viridis)
require(cowplot)
require(ComplexHeatmap)
require(stringr)
require(dplyr)
require(RColorBrewer)
source('../scripts/sc_clustering/multi_sc_base.R')
source('../scripts/multi_omics/stack_plot.R')


plot.dimension.reduction <- function(x.sp, data, tail, cluster, col_codes, alpha=1.0) {
    print(col_codes)
    for (dr in c('umap', 'tsne', 'pca')) {
        pdf(paste0(dr, '_merged_', tail, '_', cluster, '.pdf'), width=10, height=7)
        plot(DimPlot(x.sp, reduction=dr, cols=col_codes, group.by=cluster))
        dev.off()
        pdf(paste0(dr, '_merged_', tail, '_omics.pdf'), width=10, height=7)
        plot(DimPlot(x.sp, reduction=dr, group.by='orig.ident'))
        dev.off()
        pdf(paste0(dr, '_merged_', tail, '_', cluster, '_batch.pdf'), width=25, height=7)
        plot(DimPlot(x.sp, reduction=dr, cols=col_codes, split.by='batch', group.by=cluster))
        dev.off()
    }

    for (dr in c('tsne', 'umap')) {
        if (data == 'both') break
        df <- x.sp@meta.data[,c(paste0(dr, '1'), paste0(dr, '2'), cluster)]
        df[,cluster] = factor(df[,cluster], levels=names(col_codes))
        colnames(df) = c('x', 'y', 'cluster')
        print(head(df))
        g <- ggplot(df, aes(x=x, y=y, color=cluster))+geom_point(alpha=alpha, size=1)+scale_color_manual(values=col_codes)+theme_cowplot()
        pdf(paste0('original_', dr, '_', cluster, '_', tail, '.pdf'), width=10, height=7)
        plot(g)
        dev.off()
    }
}

plot.heatmap.marker.genes <- function(x.sp, data, tail, cluster, codes, id=TRUE) {
    raster_q = 5
    files <- MARKER
    annot = data.frame(x.sp@meta.data[,c(cluster, 'batch')])
    colnames(annot) = c('cluster', 'ident')
    rownames(annot) = as.character(1:dim(annot)[1])
    if (data == 'atac') {
        mat = as.matrix(t(x.sp@assays$ACTIVITY@data))
    } else {
        mat = as.matrix(t(x.sp@assays$RNA@data))
    }
    mat[is.na(mat)] = 0
    rownames(mat) = as.character(1:dim(annot)[1])
    annot = annot[rowSums(mat) > 0,]
    mat = log10(mat[rowSums(mat) > 0,]+1)
    for (i in 1:length(files)) {
        if (i < 10) next
        file = files[i]
        marker.genes <- unlist(read.table(file.path(mdir, file), header=F, stringsAsFactor=F)[,1])
        if (id) {
            column_vec <- unlist(sapply(colnames(mat), function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
            column_vec <- convert.to.symbol(column_vec, rev=FALSE)
        } else {
            column_vec <- colnames(mat)
        }
        ind = unlist(sapply(marker.genes, function(x) {
            pos = which(column_vec == x)
            if (length(pos) >= 1) return(pos[1])
        }))
        index = names(ind)
        tmat = mat[,ind]
        colnames(tmat) = index
        StackVlnPlotGG(tmat, index, annot[,'cluster'], paste0(tail, '_', i, '_', cluster), FALSE)
        tmat = filt.mat(tmat)
        col_fun = viridis(100)
        ident_col = viridis(5)[1:4]
        names(ident_col) = sort(unique(annot[,'ident']))
        row_ha = rowAnnotation(df=annot, col=list(cluster=codes, ident=ident_col))
        pdf(paste0('heatmap_marker_', tail, '_', i, '_', cluster, '.pdf'), width=10, height=15)
        h <- Heatmap(as.matrix(tmat), name="exp", show_row_names = FALSE, left_annotation=row_ha, row_split=annot[,'cluster'], border=TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
        draw(h)
        dev.off()
        for (norm in c('range', 'standardize')) {
            if (norm == 'range')
                normed_mat =  t(BBmisc::normalize(t(tmat), method=norm, range=c(0, 1), margin=2))
            else {
                normed_mat = matrix(0, nrow=dim(tmat)[1], ncol=dim(tmat)[2])
                colnames(normed_mat) = colnames(tmat)
                normed_mat[tmat > 0] = 1
                print(rowSums(normed_mat))
                saveRDS(normed_mat, file='test_mat.rds')
            }
            print(table(normed_mat))
            pdf(paste0('heatmap_marker_', tail, '_', i, '_', cluster, '_', norm, '.pdf'), width=10, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=F, left_annotation=row_ha, row_split=annot[,'cluster'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE, show_column_names=TRUE)
            draw(h)
            dev.off()
        }
    }
}

plot.heatmap.module.strength <- function(x.sp, data, tail, cluster, codes) {
    annot = data.frame(x.sp@meta.data[,c(cluster, 'batch')])
    colnames(annot) = c('cluster', 'ident')
    rownames(annot) = 1:dim(annot)[1]
    ident_col = viridis(5)[1:4]
    names(ident_col) = sort(unique(annot[,'ident']))
    row_ha = rowAnnotation(df=annot, col=list(cluster=codes, ident=ident_col))
    col_fun = NULL
    for (sp in c('human', 'mouse')) {
        mx.sp = read.module.list(sp, x.sp)
        module_mat = mx.sp@meta.data[,grepl('Module', colnames(mx.sp@meta.data))]
        pdf(paste0('heatmap_module_', tail, '_', cluster, '_', sp, '.pdf'), width=12, height=15)
        h <- Heatmap(module_mat, name="exp", show_row_names=FALSE, left_annotation=row_ha, row_split=annot[,'cluster'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 4))
        draw(h)
        dev.off()
        for (norm in c('range', 'standardize')[1]) {
            normed_mat =  t(BBmisc::normalize(t(module_mat), method=norm, range=c(0, 1)))
            pdf(paste0('heatmap_marker_', tail, '_', cluster, '_', norm, '.pdf'), width=12, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=FALSE, left_annotation=row_ha, row_split=annot[,'cluster'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 4))
            draw(h)
            dev.off()
        }
    }
}

plot.feature.module.strength <- function(x.sp, data) {
    files <- MODULE
    for (i in 1:length(files)) {
        file <- files[i]
        marker.genes <- unlist(read.table(file.path(mdir, file), header=F, stringsAsFactor=F)[,1])
        if (length(grep('ENSDNOG', marker.genes)) == length(marker.genes)) {
            features = marker.genes
        } else {
            features = convert.to.symbol(marker.genes, rev=TRUE)
        }
        features = features[unlist(sapply(features, function(x) { return(any(x == rownames(x.sp)))}))]
        key = paste0('module')
        x.sp <- AddModuleScore(object = x.sp, features = list(features), ctrl = 5, name = key)
        pdf(paste0('feature_activity_', data, '_', key, '_', i, '.pdf'))
        g <- FeaturePlot(object = x.sp, features = paste0(key, '1'))+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "OrRd"))
        plot(g)
        dev.off()
    }
}


extract.cells <- function(x.sp, data) {
    if (data == 'rna') {
        index = which(x.sp@meta.data[,'orig.ident'] == 'RNA')
    } else if (data == 'atac') {
        index = which(x.sp@meta.data[,'orig.ident'] == 'ATAC')
    } else {
        index = 1:dim(x.sp)[2]
    }
    if (length(index) == 0) return(NULL)
    else return(x.sp[,index])
}

plot.marker.mean <- function(gmat, cluster, tail, marker_file='') {
    require(ComplexHeatmap)
    require(dplyr)
    require(viridis)
    require(BBmisc)
    require(circlize)
    if (marker_file == '') {
        files <- MARKER
    } else {
        files <- c(marker_file)
    }
    column_vec <- unlist(sapply(colnames(gmat), function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
    column_vec <- convert.to.symbol(column_vec)

    for (i in 1:length(files)) {
        if (i == 5 && grepl('imputed', tail))  break # memory error
        file <- files[i]
        marker.genes <- unlist(read.table(file.path(mdir, file), header=F, stringsAsFactor=F)[,1])
        print(head(marker.genes))
        mat <- do.call(cbind, lapply(marker.genes, function(x) {
            index = which(column_vec == x)
            if (length(index) == 0)
                return(rep(0, dim(gmat)[1]))
            else if (length(index) == 1) {
                return(gmat[,index])
            } else {
                return(Matrix::colSums(gmat[,index]))
            }
        }))
        mat <- data.frame(mat)
        colnames(mat) = marker.genes
        mat <- mat[,Matrix::colSums(mat) > 0]
        mat <- cbind(mat, cluster=as.numeric(as.character(cluster)))
        data <- mat %>% group_by(cluster) %>% summarise_all(mean, na.rm = TRUE)
        mat <- data.matrix(as.data.frame(data))
        rownames(mat) <- as.numeric(as.character(mat[,1]))
        mat <- mat[,-c(1)]
        pdf(paste0('heamap_mean_', i, tail, '.pdf'))
        draw(Heatmap(mat, col=col_fun, show_row_names = TRUE))
        dev.off()
        for (norm in c('range', 'standardize')) {
            normed_mat <- t(BBmisc::normalize(t(mat), method=norm, range=c(0, 1)))
            colnames(normed_mat) = colnames(mat)
            col_fun = NULL
            pdf(paste0('heatmp_mean_norm_', i, '_', norm, tail, '.pdf'))
            draw(Heatmap(normed_mat, col=col_fun))
            dev.off()
            pdf(paste0('heatmp_mean_norm_', i, '_', norm, tail, '_ordered.pdf'))
            draw(Heatmap(normed_mat, col=col_fun, cluster_rows = FALSE, cluster_column=FALSE))
            dev.off()
        }
    }
}

extract.cluster.markers <- function(x.sp, cluster, tail) {
    require(viridis)
    if (tail == 'rna') {
        Idents(object = x.sp) = x.sp@meta.data[[cluster]]
    } else {
        Idents(object = x.sp) = x.sp[[cluster]]
    }
    x.sp.markers <- FindAllMarkers(object = x.sp, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, return.thresh = 1.01)
    top10 <- x.sp.markers %>% group_by(cluster) %>% top_n(40000, avg_logFC)
    write.table(top10, file=paste0('rank_gene_', cluster, '_', tail, '.tsv'), sep="\t")
    top10 <- x.sp.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    pdf(paste0('heatmap_top10_', tail, '.pdf'))
    if (tail == '_rna')
        d <- DoHeatmap(object = x.sp, assay='RNA', features = top10$gene, group.bar=TRUE, label = FALSE, slot='scale.data')+scale_fill_viridis()
    else
        d <- DoHeatmap(object = x.sp, assay='RNA', features = top10$gene, group.bar=TRUE, label = FALSE, slot='counts')+scale_fill_viridis()
    plot(d)
    dev.off()
}

filter.chr <- function(x.sp, max_idx) {
    chr.exclude = seqlevels(x.sp@feature)[grep("chrM", seqlevels(x.sp@feature))];
    idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
    if (length(idy) > 0) {
        x.sp@bmat = x.sp@bmat[,-idy]
    }
    x.sp = makeBinary(x.sp, mat="bmat");
    bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
    bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
    idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
    x.sp@bmat = x.sp@bmat[, idy]
    x.sp@bmat = x.sp@bmat[,order(Matrix::colSums(x.sp@bmat), decreasing=TRUE)]
    print(dim(x.sp@bmat))
    return(x.sp)
}

plot.violinplot <- function(x.sp, symbol=FALSE) {
    require(Seurat)
    se_obj <- CreateSeuratObject(counts=t(x.sp@gmat), assay="RNA")
    se_obj <- AddMetaData(object = se_obj, metadata = factor(x.sp@metaData[,'cluster']), col.name='cluster')
    se_obj <- se_obj[,se_obj@meta.data[,'nCount_RNA'] > 2500] # to extract high coverage cells
    se_obj[['RNA']]@counts <- as(log2(se_obj[['RNA']]@counts+1), "dgCMatrix")
    for (i in 1:length(files)) {
        file <- files[i]
        marker.genes <- unlist(read.table(file.path(mdir, file), header=F, stringsAsFactor=F)[,1])
        marker.genes <- convert.gene.id(marker.genes)
        label = 'gene_id'
        if (symbol == label) label = 'symbol'
        pdf(paste0('plot_violin_', i, '.pdf'), width=20, height=20)
        g <- VlnPlot(se_obj, features = marker.genes[,'gene_id'], group.by='cluster', assay='RNA', pt.size=0)
        plot(g)
        dev.off()
        pdf(paste0('plot_ridge_', i, '.pdf'), width=20, height=20)
        g <- RidgePlot(se_obj, features = marker.genes[,'gene_id'], group.by='cluster', assay='RNA', log=FALSE)
        plot(g)
        dev.off()
    }
}


single.omics.analysis <- function(fname, id=TRUE, cluster='cluster', data_label='') {
    fhead = unlist(strsplit(fname, '.', fixed=TRUE))[1]
    x.sp = readRDS(fname)
    colnames(x.sp@meta.data)[colnames(x.sp@meta.data) == cluster] = 'cluster'
    if (!id) {
        x.sp@meta.data[,'cluster'] = str_pad(x.sp@meta.data[,'cluster'], 2, pad='0')
    }
    hex_code_list = set.hex.code(x.sp)
    for (data in c('rna', 'atac', 'both')) {
        if (data_label != '') {
            if (data != data_label) next
            temp = x.sp
            print(head(temp@meta.data))
        } else {
            temp = extract.cells(x.sp, data)
            if (data == 'rna') {
                temp2 = CreateSeuratObject(counts=temp@assays$RNA@counts)
                temp2@meta.data = temp@meta.data
                print(head(temp2))
                temp = temp2
            }
        }
        if (is.null(temp)) next
        header = paste0(fhead, '_', data)
        plot.feature.module.strength(x.sp, header)
        for (cluster in c('cluster', 'rna_celltype', 'batch')) {
            if (!any(colnames(x.sp@meta.data) == cluster)) next
            hex_codes = hex_code_list[[cluster]]
            plot.heatmap.module.strength(temp, data, header, cluster, hex_codes)
            plot.heatmap.marker.genes(temp, data, header, cluster, hex_codes, id)
            extract.cluster.markers(temp, cluster, data)
        }
    }
}


multi.omics.analysis <- function(fname, cluster='cluster')
{
    require(viridis)
    x.sp = readRDS(fname)
    hex_code_list = set.hex.code(x.sp)
    fhead = unlist(strsplit(fname, '.', fixed=TRUE))[1]    
    for (data in c('atac', 'rna', 'both')) {
        temp = extract.cells(x.sp, data)
        if (data == 'atac') {
            temp = RenameAssays(object=temp, RNA='ACTIVITY')
        } else if (data == 'rna') {
            temp = manual.curation.rna.celltype(temp, input_label='cluster', output_label='rna_cluster')
        }
        for (cluster in c('cluster', 'rna_celltype', 'batch')) {
            hex_codes = hex_code_list[[cluster]]
            header = paste0(fhead, '_', data)
            plot.dimension.reduction(temp, data, header, cluster, hex_code_list[[cluster]], alpha=0.8)
        }
   }
}


args = commandArgs(trailingOnly=TRUE)
fname = args[1]
print(fname)
if (length(args) > 1) {
    if (args[2] == 'refseq') {
        id = FALSE
    } else { # ensembl
        id = TRUE
    }
}
if (length(args) > 2) {
    cluster = args[3]
} else {
    cluster = 'cluster'
}
if (length(args) > 3) {
    data = args[4]
} else {
    data = ''
}

single.omics.analysis(fname, TRUE, cluster, data)
# multi.omics.analysis(fname, cluster)



