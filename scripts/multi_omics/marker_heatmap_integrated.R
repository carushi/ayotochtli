require(parallel)
require(ggplot2)
require(circlize)
require(BBmisc)

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
load(paste0(DIR, "Xscaffolds.Rdata"))

SCRIPT_DIR <<- "../scripts/"
source(paste0(SCRIPT_DIR, "sc_clustering/multi_sc_base.R"))
source(paste0(SCRIPT_DIR, "annotation/gene_enrichment.R"))
source(paste0(SCRIPT_DIR, "multi_omics/bulk_pred_utils.R"))
source(paste0(SCRIPT_DIR, "multi_omics/sc_pred_utils.R"))

GENE_DIR <<- "../data/gene/"
ENCODE_MAT <<- file.path(GENE_DIR, "ENCODE_ChIP-seq_mat.tsv")

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

summarize.predictable.genes.recursive <- function(pred.each, gene_ann, tail) {
    bulk.each = pred.each[,!grepl('16-9_._4', colnames(pred.each))]
    results <- construct.annot.data(colnames(bulk.each), row=TRUE)
    column_ha = results[[1]]
    annot = results[[2]]
    col_list = results[[3]]
    quads = sort(unique(substr(colnames(bulk.each), 1, 4)))
    inds = sort(unique(substr(colnames(bulk.each), 1, 6)))

    df = NULL
    for (quad in quads) {
        for (ind in inds[grep(quad, inds)]) {
            print(c(quad, ind))
            filtered = filtered.mat(bulk.each, ind)
            print(colnames(filtered))
            for (i in 1:(dim(filtered)[2])) {
                df = rbind(df, c(quad, paste0(ind, '_', 1), paste0(quad, '_', 1), length(which(rowSums(filtered[,c(i),drop=F]) == 1))))
            }
            for (i in 1:(dim(filtered)[2]-1)) {
                for (j in (i+1):dim(filtered)[2]) {
                    df = rbind(df, c(quad, paste0(ind, '_', 2), paste0(quad, '_', 2), length(which(rowSums(filtered[,c(i,j)]) == 2))))
                }
            }
            print(dim(filtered)[2])
            df = rbind(df, c(quad, paste0(ind, '_', 3), paste0(quad, '_', 3), length(which(rowSums(filtered) == dim(filtered)[2]))))
        }
        filtered = filtered.mat(bulk.each, quad)
        df = rbind(df, c(quad, paste0(quad, '_5'), paste0(quad, '_', 4), length(which(rowSums(filtered) == dim(filtered)[2]))))
    }
    filtered = filtered.mat(bulk.each, '1')
    df = data.frame(df)
    print(df)
    colnames(df) = c('color', 'x', 'sum_x', 'y')
    df[,'y'] = as.numeric(as.character(df[,'y']))
    color = col_list$quad
    color[['all']] = 'black'
    pdf(paste0('barplot_recursive_overlap', tail, '.pdf'), width=20, height=5, useDingbats=F)
    g <- ggplot(df, aes(x=x, y=y, fill=color))+geom_point()+geom_boxplot()+scale_fill_manual(values=color)+theme_classic()+ facet_wrap(vars(color), scales="free", nrow=1)
    plot(g)
    dev.off()
    pdf(paste0('barplot_recursive_overlap_sum', tail, '.pdf'), width=20, height=5, useDingbats=F)
    g <- ggplot(df, aes(x=sum_x, y=y, fill=color))+geom_point()+geom_boxplot()+scale_fill_manual(values=color)+theme_classic()+ facet_wrap(vars(color), scales="free", nrow=1)
    plot(g)
    dev.off()
    pdf(paste0('barplot_recursive_overlap_sum_fix', tail, '.pdf'), width=20, height=5, useDingbats=F)
    g <- ggplot(df, aes(x=sum_x, y=y, fill=color))+geom_point()+geom_boxplot()+scale_fill_manual(values=color)+theme_classic()+ facet_wrap(vars(color), scale="free_x", nrow=1)
    plot(g)
    dev.off()
    pdf(paste0('barplot_recursive_overlap_sum_fix', tail, '_log.pdf'), width=20, height=5, useDingbats=F)
    g <- ggplot(df, aes(x=sum_x, y=y, fill=color))+geom_point()+geom_boxplot()+scale_fill_manual(values=color)+theme_classic()+ facet_wrap(vars(color), scale="free_x", nrow=1)+scale_y_log10()+coord_cartesian(ylim=(c(10, 10000)))
    plot(g)
    dev.off()
}

summarize.predictable.genes <- function(pred.each, gene_ann) {
    require(dplyr)
    for (data_label in c('sc', 'bulk')) {
        if (data_label == 'sc') {
            bulk.each = pred.each[,grepl('16-9_', colnames(pred.each))]
        } else {
            bulk.each = pred.each[,!grepl('16-9_._4', colnames(pred.each))]
        }
        temp = bulk.each
        results <- construct.annot.data(colnames(temp), row=TRUE)
        column_ha = results[[1]]
        annot = results[[2]]
        col_list = results[[3]]
        odf = NULL
        gdf = NULL
        for (i in 1:(dim(temp)[2]-1)) {
            for (j in (i+1):dim(temp)[2]) {
                filtered = temp[!is.na(temp[,i]) & !is.na(temp[,j]),c(i,j)]
                overlap = length(which(rowSums(filtered) == 2))
                group1 = length(which(filtered[,1] > 0))
                group2 = length(which(filtered[,2] > 0))
                total = dim(filtered)[1]
                pvalue = phyper(overlap-1, group2, total-group2, group1, lower.tail= FALSE)
                pfisher = fisher.test(matrix(c(overlap, group2-overlap, group1-overlap, total-group2-group1+overlap), 2, 2), alternative='greater')
                odf = rbind(odf, c(overlap/total, overlap, group1, group2, total, pvalue))
                odf = rbind(odf, c(overlap/total, overlap, group2, group1, total, pvalue))
                gdf = rbind(gdf, c(colnames(temp)[i], colnames(temp)[j]))
                gdf = rbind(gdf, c(colnames(temp)[j], colnames(temp)[i]))
            }
        }
        odf = data.frame(odf, gdf)
        colnames(odf) = c('jaccard', 'overlap', 'group1', 'group2', 'total', 'pvalue', 'group1_name', 'group2_name')
        odf = add.quad.names(odf)
        print(head(odf))
        for (comparison in c('time', 'within', 'across')) {
            if (data_label == 'sc') {
                if (comparison == 'across') next
                xlabel =  'ind'
                color = col_list$ident
                names(color) = paste0('16-9_', names(color))
                tdf = data.frame(odf[odf[,comparison] == TRUE, c('ind', 'quad', 'jaccard', 'overlap')])
            } else {
                xlabel = 'quad'
                color = col_list$quad
                tdf = data.frame(odf[odf[,comparison] == TRUE, c('ind', 'quad', 'jaccard', 'overlap')]) %>% group_by(ind, quad) %>% summarise_all(mean)
            }
            pdf(paste0('boxplot_', data_label, '_', comparison, '_jaccard.pdf'), useDingbats=FALSE)
            g <- ggplot(tdf, aes_string(x=xlabel, y='jaccard', fill=xlabel))+geom_boxplot()+geom_point()+scale_fill_manual(values=color)+theme_classic()
            plot(g)
            dev.off()
            pdf(paste0('boxplot_', data_label, '_', comparison, '_overlap.pdf'), useDingbats=F)
            g <- ggplot(tdf, aes_string(x=xlabel, y='overlap', fill=xlabel))+geom_boxplot()+geom_point()+scale_fill_manual(values=color)+theme_classic()
            plot(g)
            dev.off()
        }
    }
}

plot.heatmap.gene.set <- function(mat, tail, rank_flag, gene_ann, index, ind, top_annotated=FALSE, gene_flag=TRUE) {
    raster_q = 5
    results <- construct.annot.data(rownames(mat), row=TRUE)
    column_ha = results[[1]]
    annot = results[[2]]
    width=8
    height=12
    tmat = mat[,ind]
    colnames(tmat) = index
    print(dim(tmat))
    if (any(dim(tmat) == 0)) return()
    tmat = filt.mat(tmat)
    col_fun = viridis(100)
    ident_col = viridis(5)[1:4]
    names(ident_col) = sort(unique(annot[,'ident']))
    if (top_annotated) {
        top_ann = gene.color.annotation(gene_ann, colnames(tmat))
    }
    if (rank_flag) {
        if (top_annotated) {
            pdf(paste0('heatmap_marker_', tail, '.pdf'), width=10, height=15)
            h <- Heatmap(as.matrix(tmat), name="exp", show_row_names = FALSE, show_column_names=FALSE, top_annotation=top_ann, left_annotation=column_ha, row_split=annot[,'quad'], border=TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
            pdf(paste0('heatmap_marker_', tail, '_ordered.pdf'), width=10, height=15)
            h <- Heatmap(as.matrix(tmat), name="exp", show_row_names = FALSE, show_column_names=FALSE, top_annotation=top_ann, cluster_rows=FALSE, left_annotation=column_ha, row_split=annot[,'quad'], border=TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
        } else {
            if (gene_flag) {
                column_split = NULL
                top_ann = NULL
            } else {
                top_ann = results[[4]]
                column_split=annot[,'quad']
            }
            pdf(paste0('heatmap_marker_', tail, '_wo.pdf'), width=10, height=15)
            # col_fun = colorRamp2(seq(min(tmat), 1.0, length.out=8), viridis(8))
            h <- Heatmap(as.matrix(tmat), name="exp", show_row_names = FALSE, show_column_names=(dim(tmat)[2] < 60), top_annotation=top_ann, left_annotation=column_ha, row_split=annot[,'quad'], column_split=column_split, border=TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
            pdf(paste0('heatmap_marker_', tail, '_wo_ordered.pdf'), width=10, height=15)
            h <- Heatmap(as.matrix(tmat), name="exp", show_row_names = FALSE, show_column_names=(dim(tmat)[2] < 60), top_annotation=top_ann, cluster_rows=FALSE, cluster_columns=TRUE, left_annotation=column_ha, row_split=annot[,'quad'], column_split=column_split, border=TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
        }
    }
    if (rank_flag) return()
    for (norm in c('log', 'range', 'standardize')) {
        if (norm == 'log') {
            normed_mat = log10(tmat+1)
        } else if (norm == 'range')
            normed_mat =  BBmisc::normalize(tmat, method=norm, range=c(0, 1), margin=2)
        else {
            normed_mat = matrix(0, nrow=dim(tmat)[1], ncol=dim(tmat)[2])
            colnames(normed_mat) = colnames(tmat)
            normed_mat[tmat > 0] = 1
        }
        if (top_annotated) {
            pdf(paste0('heatmap_marker_', tail, '_', norm, '.pdf'), width=10, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=F, left_annotation=column_ha, show_column_names=(dim(normed_mat)[2] < 60), top_annotation=top_ann, row_split=annot[,'quad'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
            pdf(paste0('heatmap_marker_', tail, '_', norm, '_ordered.pdf'), width=10, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=F, left_annotation=column_ha, show_column_names=(dim(normed_mat)[2] < 60), top_annotation=top_ann, cluster_rows=FALSE, row_split=annot[,'quad'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
        } else {
            pdf(paste0('heatmap_marker_', tail, '_', norm, '_wo.pdf'), width=10, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=F, show_column_names=(dim(normed_mat)[2] < 60), left_annotation=column_ha, row_split=annot[,'quad'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE)
            draw(h)
            dev.off()
            pdf(paste0('heatmap_marker_', tail, '_', norm, '_ordered.pdf'), width=10, height=15)
            h <- Heatmap(normed_mat, name="exp", show_row_names=F, show_column_names=(dim(normed_mat)[2] < 60), left_annotation=column_ha, row_split=annot[,'quad'], border = TRUE, col=col_fun, column_names_gp = gpar(fontsize = 8), raster_quality = raster_q, use_raster=TRUE, cluster_rows=FALSE, cluster_columns=FALSE)
            draw(h)
            dev.off()
        }
    }
}

plot.heatmap.marker.genes <- function(mat, header, id=TRUE, rank_flag=FALSE, gene_ann=NULL, marker_list=c('all', 'marker', 'annotation'), gene_flag=TRUE) {
    files <- MARKER
    mdir <- "../data/gene/"
    mat = t(mat)
    print(head(mat))
    for (marker_set in marker_list) {
        if (marker_set == 'all') {
            i = 0
            ind = 1:dim(mat)[2]
            index = colnames(mat)
            plot.heatmap.gene.set(mat, paste0(header, '_', i), rank_flag, gene_ann, index, ind, ((marker_set != 'marker') & !is.null(gene_ann)), gene_flag=gene_flag)
        } else if (marker_set == 'marker') {
            for (i in 1:length(files)) {
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
                plot.heatmap.gene.set(mat, paste0(header, '_', i), rank_flag, gene_ann, index, ind, ((marker_set != 'marker') & !is.null(gene_ann)))
            }
        } else {
            for (col in colnames(gene_ann)) {
                index = rownames(gene_ann)[(gene_ann[,col] == '1') & (rownames(gene_ann) %in% colnames(mat))]
                ind = which(colnames(mat) %in% index)
                i = col
                plot.heatmap.gene.set(mat, paste0(header, '_', i), rank_flag, gene_ann, index, ind, ((marker_set != 'marker') & !is.null(gene_ann)))
            }
        }
    }
}

celltype.color.code <- function() {
    gene_table = read.table(file.path(GENE_DIR, 'rank_gene_rna_celltype_rna.tsv'), header=T, sep="\t", stringsAsFactor=FALSE)
    ct_list = sort(unique(as.character(gene_table[,6])))
    hex_codes2 <- rev(hue_pal()(length(ct_list)))
    names(hex_codes2) = ct_list
    return(hex_codes2)
}

obtain.gene.annotation.tf <- function(pred.each, gt=1, all_types=FALSE, monocyte_only=FALSE, convert2factor=TRUE) {
    predictability = pred.each[,grep('16-9_2_[1-3]', colnames(pred.each))]
    colnames(predictability) = paste0('predictability_', 1:3)
    tf_list = c('EGR1_MCF7_HG19', 'TCF12_SKNSH_HG19', 'RELA_GM18526_HG19', 'NFIC_ECC1_HG19', 'CEBPD_K562_HG19', 'E2F4_MELCELLLINE_MM9', 'CEBPB_MYOCYTE_MM9', 
            'NFIC_SKNSH_HG19', 'FOXM1_ECC1_HG19', 'CEBPB_GM12878_HG19')
    tf_mat = read.table(ENCODE_MAT, header=T, sep="\t")
    gene_ann = predictability
    for (tf in tf_list) {
        tf_target = tf_mat[,tf, drop=F]
        tf_filter = rownames(tf_target)[tf_target[,1] == 1]
        gene_ann = cbind(gene_ann, tf=unlist(sapply(rownames(predictability), function(x) {
            return(any(x == tf_filter))
        })))
        colnames(gene_ann)[dim(gene_ann)[2]] = paste0('tf_', tf)
    }
    if (all_types) {
        gene_ann = cbind(gene_ann, tf_alltypes=unlist(sapply(rowSums(gene_ann[,grep('tf_', colnames(gene_ann))], na.rm=T), function(x) min(1, x))))
    }
    for (omics in c('atac', 'rna')) {
        if (omics == 'atac') {
            gene_table = read.table(file.path(GENE_DIR, 'atac_rna_celltype_peak_mat.tsv'), header=T, sep="\t", stringsAsFactor=FALSE)
            gene_table = gene_table[rownames(gene_ann),]
            if (monocyte_only) gene_table = gene_table[,c('CD14+ Monocytes'), drop=F]
            colnames(gene_table) = paste0('dar_', colnames(gene_table))
            gene_ann = cbind(gene_ann, gene_table)
            if (!monocyte_only && all_types) {
                gene_ann = cbind(gene_ann, dar_alltypes=unlist(sapply(rowSums(gene_table, na.rm=T), function(x) {if(x > 0) return(1); return(0)})))
            }
        } else {
            gene_table = read.table(file.path(GENE_DIR, 'rank_gene_rna_celltype_rna.tsv'), header=T, sep="\t", stringsAsFactor=FALSE)
            if (monocyte_only) {
                gene_table = gene_table[gene_table[,6] == 'CD14+ Monocytes' | gene_table[,6] == 'Monocytes',]
                gene_filter = unique(gene_table[,7])
                gene_ann = cbind(gene_ann, ct=unlist(sapply(rownames(predictability), function(x) {
                    return(any(x == gene_filter))
                })))
            } else {
                for (ct in sort(unique(as.character(gene_table[,6])))) {
                    if (ct == 'Undefined') next
                    ttable = gene_table[gene_table[,6] == ct,]
                    gene_filter = unique(ttable[,7])
                    gene_ann = cbind(gene_ann, unlist(sapply(rownames(predictability), function(x) {
                        return(any(x == gene_filter))
                    })))
                    colnames(gene_ann)[dim(gene_ann)[2]] = paste0('ct_', gsub("_+", "_", gsub("\\+", "_", gsub("\\s+", "_", as.character(ct)))))
                }
                if (!monocyte_only && all_types) {
                    gene_ann = cbind(gene_ann, ct_alltypes=unlist(sapply(rowSums(gene_ann[,grep(paste0('ct_'), colnames(gene_ann))], na.rm=T), function(x) {if(x > 0) return(1); return(0)})))
                }
            }
        }
    }
    gene_ann = data.frame(gene_ann)
    gene_ann[is.na(gene_ann)] = 0
    if (convert2factor) {
        for (i in 1:dim(gene_ann)[2]) {
            gene_ann[,i] = factor(as.character(gene_ann[,i]), levels=c('0', '1'))
        }
    }
    return(gene_ann)
}

obtain.gene.annotation <- function(pred.each, target=1, gt=1, all_ct=FALSE) {
    predictability = pred.each[,grep('16-9_2_[1-3]', colnames(pred.each))]
    colnames(predictability) = paste0('predictability_', 1:3)
    print(target)
    tf = c('EGR1_MCF7_HG19', 'TCF12_SKNSH_HG19', 'RELA_GM18526_HG19', 'NFIC_ECC1_HU19', 'CEBPD_K562_HG19', 'E2F4_MELCELLLINE_MM9', 'CEBPB_MYOCYTE_MM9', 
            'NFIC_SKNSH_HG19', 'FOXM1_ECC1_HG19', 'CEBPB_GM12878_HG19')[target]
    tf_mat = read.table(ENCODE_MAT, header=T, sep="\t")
    print(tf)
    print(any(colnames(tf_mat) == tf))
    tf_target = tf_mat[,tf, drop=F]
    print(head(tf_target))
    print(sum(tf_target))
    tf_filter = rownames(tf_target)[tf_target[,1] == 1]
    gene_ann = cbind(predictability, tf=unlist(sapply(rownames(predictability), function(x) {
        return(any(x == tf_filter))
    })))
    if (all_ct) { # ATAC DAR
        gene_table = read.table(file.path(GENE_DIR, 'rank_gene_rna_celltype_rna.tsv'), header=T, sep="\t", stringsAsFactor=FALSE)
        for (ct in sort(unique(as.character(gene_table[,6])))) {
            if (ct == 'Undefined') next
            ttable = gene_table[gene_table[,6] == ct,]
            gene_filter = unique(ttable[,7])
            gene_ann = cbind(gene_ann, unlist(sapply(rownames(predictability), function(x) {
                return(any(x == gene_filter))
            })))
            colnames(gene_ann)[dim(gene_ann)[2]] = paste0('ct_', gsub("_+", "_", gsub("\\+", "_", gsub("\\s+", "_", as.character(ct)))))
        }
    } else {
        gene_table = read.table(file.path(GENE_DIR, 'rank_gene_rna_celltype_rna.tsv'), header=T, sep="\t", stringsAsFactor=FALSE)
        gene_table = gene_table[gene_table[,6] == 'CD14+ Monocytes' | gene_table[,6] == 'Monocytes',]
        gene_filter = unique(gene_table[,7])
        gene_ann = cbind(gene_ann, ct=unlist(sapply(rownames(predictability), function(x) {
            return(any(x == gene_filter))
        })))
    }
    gene_ann = data.frame(gene_ann)
    gene_ann[is.na(gene_ann)] = 0
    for (i in 1:dim(gene_ann)[2]) {
        gene_ann[,i] = factor(as.character(gene_ann[,i]), levels=c('0', '1'))
    }
    return(gene_ann)
}

gene.color.annotation <- function(gene_ann, gene_list) {
    col_list = list()
    index = 1
    call = 5
    for (i in 1:3) {
        type_col = c('white', rainbow(call)[index])
        index = index+1
        names(type_col) = c('0', '1')
        col_list[[paste0('predictability_', i)]] = type_col
    }
    for (label in c('tf', 'ct')) {
        type_col = c('white', rainbow(call)[index])
        index = index+1
        if (any(colnames(gene_ann) == label)) {
            # names(type_col) = unique(gene_ann[,label])
            names(type_col) = c('0', '1')
            col_list[[label]] = type_col
        }
    }
    ct_list = colnames(gene_ann)[grep('ct_', colnames(gene_ann))]
    ct_list = ct_list[order(ct_list)]    
    for (i in 1:length(ct_list)) {
        if (length(ct_list) == 0) break
        #  hex_codes2 <- rev(hue_pal()(length(unique(coembed@meta.data[,'rna_celltype']))))
        hex_codes2 <- rev(hue_pal()(length(ct_list)))
        ct = ct_list[i]
        type_col = c('white', hex_codes2[i])
        names(type_col) = c('0', '1')
        col_list[[ct]] = type_col
    }
    tf_list = colnames(gene_ann)[grep('tf_', colnames(gene_ann))]
    for (i in 1:length(tf_list)) {
        if (length(tf_list) == 0) break
        tf = tf_list[i]
        type_col = c('white', rainbow(length(tf_list))[i])
        names(type_col) = c('0', '1')
        col_list[[tf]] = type_col
    }
    tf_list = colnames(gene_ann)[grep('dar_', colnames(gene_ann))]
    for (i in 1:length(tf_list)) {
        if (length(tf_list) == 0) break
        tf = tf_list[i]
        type_col = c('white', rainbow(length(tf_list))[i])
        names(type_col) = c('0', '1')
        col_list[[tf]] = type_col
    }
    for (col in colnames(gene_ann)) {
        print(col)
        print(col_list[[col]])
        print(unique(gene_ann[gene_list,col]))
    }
    column_ann = HeatmapAnnotation(df=gene_ann[gene_list,], col=col_list, which='column', show_annotation_name=FALSE)
    return(column_ann)
}

plot.hist.similarity.within <- function(cor_bulk) {
    makeTransparent <- function(col_name) {
        rgb_value = t(col2rgb(col_name, alpha = TRUE))
        rgb_value[4] = 80
        rgb_value = rgb_value/255
        return(rgb(rgb_value[1], rgb_value[2], rgb_value[3], rgb_value[4]))
    }
    sim.same.all = c()
    sim.within.all = c()
    sim.across.all = c()
    for (quad in 1:5) {
        start = (quad-1)*12+1
        end = quad*12
        for (ind in 1:4) {
            istart = start+(ind-1)*3
            iend = istart+2
            for (i in istart:iend) {
                for (j in (i+1):dim(cor_bulk)[2]) {
                    if (i == dim(cor_bulk)[2]) next
                    # print(c(i, j, (j<=iend), (j<=end)))
                    if (j <= iend) sim.same.all <- c(sim.same.all, cor_bulk[i,j])
                    else if (j <= end) sim.within.all <- c(sim.within.all, cor_bulk[i,j])
                    else sim.across.all <- c(sim.across.all, cor_bulk[i,j])
                }
            }
        }
    }
    print(length(sim.same.all))
    print(length(sim.within.all))
    print(length(sim.across.all))

    h1 = hist(sim.across.all, breaks=30, plot=F) 
    h2 = hist(sim.within.all, breaks=20, plot=F) 
    h3 = hist(sim.same.all, breaks=20, plot=F)
    print('breaks')
    print(length(h1$breaks))
    print(length(h2$breaks))
    print(length(h3$breaks))
    print(h3)
    h3$density = h3$counts/sum(h2$counts)
    h2$density = h2$counts/sum(h2$counts)
    h1$density = h1$counts/sum(h1$counts)
    ylim= range(c( h1$density, h2$density, h3$density) )
    xlim= range(c( h1$breaks, h2$breaks, h3$breaks) )
    pdf('cor_within_dist.pdf')
    plot(h1, freq=F, col=makeTransparent("blue"), border=NA, xlim=xlim, ylim=ylim, xlab="Sample-sample correlations", main="")
    plot(h2, freq=F,add=T, border=NA, col=makeTransparent("purple") )
    plot(h3, freq=F,add=T, border=NA, col=makeTransparent("darkgreen") )
    abline(v= mean(sim.across.all), lwd=2, col="blue" )
    abline(v= mean(sim.within.all), lwd=2, col="purple" )
    abline(v= mean(sim.same.all), lwd=2, col="darkgreen" )
    dev.off()
    print(c('across.all', mean(sim.across.all)))
    print(c('within.all', mean(sim.within.all)))
    print(c('same.all',   mean(sim.same.all)))
    print(wilcox.test( sim.across.all, sim.within.all ))
    print(wilcox.test( sim.within.all, sim.same.all ))
    print(wilcox.test( sim.across.all, sim.same.all ))
}

debug.rank <- function(bulk_mat, all_bulk_mat, rank_test, rank_train) {
    temp = all_bulk_mat[,grep('16.9_', colnames(all_bulk_mat))]
    print(temp[1,grep('train', colnames(temp))])
    print(temp[2,grep('train$', colnames(temp))])
    print(temp[3,grep('train$', colnames(temp))])
    head(temp[,grep('train', colnames(temp))])
    for (row in 1:dim(temp)[1]) {
        print(temp[row,1:4])
        for (time in 1:3) {
            others = c(1:3)
            others = others[others != time]
            print(bulk_mat[row,grep(paste0('16.90[1-4]_', time), colnames(bulk_mat))])
            print(rank_test[row,grep(paste0('16.9_[1-4]_', time), colnames(rank_test))])
            exp_order = rank(bulk_mat[row,grep(paste0('16.90[1-4]_', time), colnames(bulk_mat))])
            print(exp_order)
            for (i in 1:4) {
                stopifnot(exp_order[i] == rank_test[row,grep(paste0('16.9_[1-4]_', time), colnames(rank_test))][i])
            }
            print(rank(bulk_mat[row,grep(paste0('16.90[1-4]_', others[1]), colnames(bulk_mat))]))
            print(rank(bulk_mat[row,grep(paste0('16.90[1-4]_', others[2]), colnames(bulk_mat))]))
            ave_rank = (rank(bulk_mat[row,grep(paste0('16.90[1-4]_', others[1]), colnames(bulk_mat))])+rank(bulk_mat[row,grep(paste0('16.90[1-4]_', others[2]), colnames(bulk_mat))]))/2
            ave_rank = rank(ave_rank)
            print(rank(ave_rank))
            for (i in 1:4) {
                stopifnot(ave_rank[i] == rank_train[row,grep(paste0('16.9_[1-4]_', time), colnames(rank_train))][i])
            }
        }
    }
}

combine.rank.data <- function(rank_test, rank_train)
{
    rank_combined <- do.call(cbind, lapply(1:dim(rank_test)[2], function(x) {
        mat = cbind(rank_test[,x], rank_train[,x])
        colnames(mat) = paste0(rep(colnames(rank_test)[x], 2), '_', c('test', 'train'))
        return(mat)
    }))
    rank_combined = data.frame(rank_combined)
    colnames(rank_combined) = unlist(sapply(colnames(rank_combined), function(x) {
        return(substr(x, 2, nchar(x)))
    }))
    rownames(rank_combined) = rownames(rank_test)

    rank_combined = rank_combined[,order(colnames(rank_combined))]
    return(rank_combined)
}

plot.pie.chart <- function(seta, setb, group1, overlap, header, pvalue) {
    df = data.frame(value=c(group1-overlap, overlap), group=c(paste0('1 Unique'), paste0('2 Intesrect')))
    pdf(paste0('pie_chart_', header, '.pdf'), useDingbats=FALSE)
    g <- ggplot(df, aes(x="", y=value, fill=group))+geom_col()+coord_polar(theta="y")+scale_fill_manual(values=c('gray', 'red'))+theme_classic()+ggtitle(paste0('pvalue=', pvalue))
    plot(g)
    dev.off()
}

obtain.tf.target <- function(gene_list, tf_list, all_types=TRUE) {
    tf_mat = read.table(ENCODE_MAT, header=T, sep="\t")
    mat = NULL
    for (tf in tf_list) {
        tf_target = tf_mat[,tf, drop=F]
        tf_filter = rownames(tf_target)[tf_target[,1] == 1]
        mat = cbind(mat, tf=unlist(sapply(gene_list, function(x) {
            return(any(x == tf_filter))
        })))
        colnames(mat)[dim(mat)[2]] = paste0('tf_', tf)
    }
    if (all_types) {
        mat = cbind(mat, tf_alltypes=unlist(sapply(rowSums(mat[,grep('tf_', colnames(mat))], na.rm=T), function(x) min(1, x))))
    }
    return(mat)

}

plot.overlap.bar <- function(gene_ann, exp_mat, col1, color_pattern, header, prefix)
{
    results = NULL
    total = NULL
    for (col2 in colnames(gene_ann)) {
        if (col1 == col2) next
        if (!grepl(header, col2)) next
        print(c(col1, col2))
        index = (!is.na(gene_ann[,col1]) & !is.na(gene_ann[,col2]))
        stopifnot(length(which(index)) == dim(pred.each)[1])
        filtered = data.matrix(gene_ann[index, c(col1, col2)])
        mat = cbind(overlap=c('1 not', '2 common')[as.numeric(rowSums(filtered) == 2)+1], rank=exp_mat[index], gene_set=rep(col2, dim(filtered)[1]))
        mat = mat[mat[,'overlap'] != '1 not',,drop=F]
        # print(mat)
        if (color_pattern == 'overlap') {
            freq = data.frame(mat) %>% group_by(overlap, gene_set) %>% summarise(count=n())
        } else {
            freq = data.frame(mat) %>% group_by(overlap, rank, gene_set) %>% summarise(count=n())
        }
        # print('combine')
        # print(header)
        # print(head(results))
        # print(head(freq))
        results = rbind(results, freq)
        total = rbind(total, c(col2, length(which(rowSums(filtered) == 2))))
        print(total)
    }
    results = data.frame(results)
    gene_set_list = total[order(as.numeric(total[,2]), decreasing=TRUE),1]
    results = cbind(results, x=unlist(sapply(results[,'gene_set'], function(x) {
        return(paste0(str_pad(which(gene_set_list == x)[1], 2, pad='0'), '_', x))
    })))
    results[,'gene_set'] = factor(results[,'gene_set'], levels=gene_set_list)
    pdf(paste0('barplot_overlap_', prefix, '.pdf'), width=10, height=5)
    if (color_pattern == 'overlap') {
        g <- ggplot(results, aes(x=x, y=count, fill=overlap))+scale_fill_manual(values=c('violetred'), drop=F)
    } else {
        results[,'rank'] = factor(results[,'rank'], levels=sort(unique(exp_mat)))
        print(unique(results[,'rank']))
        # results[,'rank'] = factor(results[,'rank'], levels=sort(as.numeric(unique(results[,'rank']))))
        g <- ggplot(results, aes(x=x, y=count, fill=rank))+scale_fill_viridis(discrete=TRUE)
    }
    g <- g+geom_bar(stat='identity', position='stack', color='black')+theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot(g)
    dev.off()
    return(results)
}

plot.overlap.bar.all <- function(gene_ann, exp_mat, col1, color_pattern, header, results) {
    mat = cbind(overlap=rep('2 common', sum(gene_ann[,col1])), rank=exp_mat[gene_ann[,col1] > 0], gene_set=rep('all', sum(gene_ann[,col1])))
    freq = data.frame(mat) %>% group_by(overlap, rank, gene_set) %>% summarise(count=n())
    freq = cbind(data.frame(freq), x=paste0('00_', data.frame(freq)[,'gene_set']))
    print(head(freq))
    print(head(results))
    # results[,'rank'] = as.numeric(as.character(results[,'rank']))
    results = rbind(freq, results)
    # results[,'rank'] = factor(results[,'rank'], levels=sort(unique(results[,'rank'])))
    pdf(paste0('barplot_overlap_', header, '_', color_pattern, '_all.pdf'), width=10, height=5)
    g <- ggplot(results, aes(x=x, y=count, fill=rank))+scale_fill_viridis(discrete=TRUE)
    g <- g+geom_bar(stat='identity', position='stack', color='black')+theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot(g)
    dev.off()
}

visualize.annotation.overlap.whole.data <- function(pred.each, exp_all)
{
    require(dplyr)
    require(viridis)
    require(stringr)
    gene_ann = obtain.gene.annotation.tf(pred.each, monocyte_only=FALSE, all_types=TRUE, convert2factor=FALSE)
    gene_ann[is.na(gene_ann)] = 0
    gene_ann = gene_ann[,!grepl('tf', colnames(gene_ann))]
    gene_ann = gene_ann[,!grepl('predictability', colnames(gene_ann))]
    print(head(pred.each))
    print(head(exp_all))
    print(pred.each[!(rownames(pred.each) %in% rownames(exp_all)),])
    stopifnot(all(unlist(sapply(1:dim(pred.each)[1], function(x) {
        return(rownames(pred.each)[x] == rownames(exp_all)[x])
    }))))
    freq_matrix = NULL
    for (quad in c('15-5', '12-1', '16-2', '16-3', '16-9')) {
        for (id in 1:4) {
            for (time in 1:3) {
                exp_mat = exp_all[,paste0(quad, '0', id, '_', time)]
                exp_mat = paste0('rank_', as.character(exp_mat*10))
                col1 = paste0(quad, '_', id, '_', time)
                tf_list = read.table(paste0(GENE_DIR, "ENCODE_ChIP-seq_marker_", quad, "_", id, "_", time, "_", ".tsv"), header=T, sep="\t")
                top_tfs = tf_list[,2][1:10]
                tf_mat = obtain.tf.target(rownames(pred.each), top_tfs)
                temp_ann = cbind(gene_ann, tf_mat, pred.each[,col1])
                temp_ann[is.na(temp_ann)] = 0
                colnames(temp_ann)[dim(temp_ann)[2]]  = paste0('predictability_', time)
                for (color_pattern in c('overlap', 'exp')) {
                    # for (header in c('tf', 'ct', 'dar')) {
                    #     results = plot.overlap.bar(temp_ann, exp_mat, colnames(temp_ann)[dim(temp_ann)[2]], color_pattern, header, paste0(header, '_', quad, '_', id, '_', time, '_', color_pattern))
                    #     if (color_pattern == 'exp')
                    #         plot.overlap.bar.all(temp_ann, exp_mat, colnames(temp_ann)[dim(temp_ann)[2]], color_pattern, paste0(header, '_', quad, '_', id, '_', time), results)
                    # }
                    if (color_pattern == 'overlap') {
                        selected = temp_ann[,c('ct_CD14_Monocytes', 'ct_alltypes', 'tf_alltypes', 'dar_alltypes')]
                        selected = selected[temp_ann[,colnames(temp_ann)[dim(temp_ann)[2]]] == 1,]
                        for (i in 1:(dim(selected)[2]-1)) {
                            selected[selected[,i] > 0, (i+1):(dim(selected)[2])] = 0
                        }
                        selected = cbind(selected, others=unlist(sapply(rowSums(selected), function(x){return(c(1, 0)[x+1])})))
                        colnames(selected) = paste0(rev(1:dim(selected)[2]), '_', colnames(selected))
                        data = data.frame(group=colnames(selected), value=colSums(selected))
                        if (quad == '16-9' && id == 2) {
                            print(head(data))
                            pdf(paste0('pie_chart_overlap_', quad, '_', id, '_', time, '.pdf'), useDingbats=FALSE)
                            g <- ggplot(data, aes(x="", y=value, fill=group))+geom_col()+coord_polar(theta="y")+theme_classic()
                            plot(g)
                            dev.off()
                        }
                        data[,'value'] = data[,'value']/sum(data[,'value'])
                        # data[,'group'] = str_replace(data[,'group'], '[0-9]')
                        freq_matrix = rbind(freq_matrix, cbind(data, quad=rep(quad, dim(data)[1]), id=rep(id, dim(data)[1]), time=rep(time, dim(data)[1])))
                    }
                }
                print(head(freq_matrix))
            }
        }
        quad_freq = freq_matrix[freq_matrix[,'quad'] == quad,] %>% select(c(group, value)) %>% group_by(group) %>% summarize_all(mean)
        quad_freq = data.frame(quad_freq)
        pdf(paste0('pie_chart_overlap_', quad, '.pdf'), useDingbats=FALSE)
        g <- ggplot(quad_freq, aes(x="", y=value, fill=group))+geom_col()+coord_polar(theta="y")+theme_classic()
        plot(g)
        dev.off()
    }
    quad_freq = freq_matrix %>% select(c(group, value)) %>% group_by(group) %>% summarize_all(mean)
    quad_freq = data.frame(quad_freq)
    pdf(paste0('pie_chart_overlap_all.pdf'), useDingbats=FALSE)
    g <- ggplot(quad_freq, aes(x="", y=value, fill=group))+geom_col()+coord_polar(theta="y")+theme_classic()
    plot(g)
    dev.off()
    write.table(freq_matrix, file='pred_gene_ratio.tsv', sep="\t")
}

visualize.annotation.overlap <- function(pred.each, exp_mat) {
    require(dplyr)
    require(viridis)
    require(stringr)
    # exp_mat[is.na(exp_mat)] = 0
    exp_mat = paste0('rank_', as.character(exp_mat*10))
    print(unique(exp_mat))
    gene_ann = obtain.gene.annotation.tf(pred.each, monocyte_only=FALSE, all_types=TRUE, convert2factor=FALSE)
    gene_ann[is.na(gene_ann)] = 0
    stopifnot(all(unlist(sapply(1:dim(pred.each)[1], function(x) {
        return(rownames(pred.each)[x] == rownames(exp_mat)[x])
    }))))
    for (color_pattern in c('overlap', 'exp')[2]) {
        for (time in 1:3) {
            for (header in c('tf', 'ct', 'dar')) {
                results = NULL
                total = NULL
                col1 = paste0('predictability_', time)
                temp_ann = gene_ann[,grep(paste0(header, '_'), colnames(gene_ann))]
                print(c(col1, sum(gene_ann[,col1])))
                for (col2 in colnames(temp_ann)) {
                    if (col1 == col2) next
                    print(c(col1, col2))
                    index = (!is.na(gene_ann[,col1]) & !is.na(gene_ann[,col2]))
                    stopifnot(length(which(index)) == dim(pred.each)[1])
                    filtered = data.matrix(gene_ann[index, c(col1, col2)])
                    mat = cbind(overlap=c('1 not', '2 common')[as.numeric(rowSums(filtered) == 2)+1], rank=exp_mat[index], gene_set=rep(col2, dim(filtered)[1]))
                    mat = mat[mat[,'overlap'] != '1 not',]
                    if (color_pattern == 'overlap') {
                        freq = data.frame(mat) %>% group_by(overlap, gene_set) %>% summarise(count=n())
                    } else {
                        freq = data.frame(mat) %>% group_by(overlap, rank, gene_set) %>% summarise(count=n())
                    }
                    results = rbind(results, freq)
                    total = rbind(total, c(col2, length(which(rowSums(filtered) == 2))))
                    print(total)
                }
                results = data.frame(results)
                gene_set_list = total[order(as.numeric(total[,2]), decreasing=TRUE),1]
                results = cbind(results, x=unlist(sapply(results[,'gene_set'], function(x) {
                    return(paste0(str_pad(which(gene_set_list == x)[1], 2, pad='0'), '_', x))
                })))
                all_header = paste0(header, '_', 'alltypes')
                results[,'gene_set'] = factor(results[,'gene_set'], levels=gene_set_list)
                results[,'rank'] = factor(results[,'rank'], levels=sort(unique(exp_mat)))
                print(unique(results[,'rank']))
                # exit()
                # results[,'rank'] = factor(as.numeric(as.character(results[,'rank'])))
                pdf(paste0('barplot_overlap_', header, '_', time, '_', color_pattern, '.pdf'), width=10, height=5)
                if (color_pattern == 'overlap') {
                    g <- ggplot(results, aes(x=x, y=count, fill=overlap))+scale_fill_manual(values=c('violetred'), drop=F)
                } else {
                    # results[,'rank'] = factor(results[,'rank'], levels=sort(as.numeric(unique(results[,'rank']))))
                    g <- ggplot(results, aes(x=x, y=count, fill=rank))+scale_fill_viridis(discrete=TRUE)
                }
                g <- g+geom_bar(stat='identity', position='stack', color='black')+theme_classic()+
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                plot(g)
                dev.off()
                if (color_pattern == 'exp') {
                    mat = cbind(overlap=rep('2 common', sum(gene_ann[,col1])), rank=exp_mat[gene_ann[,col1] > 0], gene_set=rep('all', sum(gene_ann[,col1])))
                    freq = data.frame(mat) %>% group_by(overlap, rank, gene_set) %>% summarise(count=n())
                    freq = cbind(data.frame(freq), x=paste0('00_', data.frame(freq)[,'gene_set']))
                    print(head(freq))
                    print(head(results))
                    # results[,'rank'] = as.numeric(as.character(results[,'rank']))
                    results = rbind(freq, results)
                    # results[,'rank'] = factor(results[,'rank'], levels=sort(unique(results[,'rank'])))
                    pdf(paste0('barplot_overlap_', header, '_', time, '_', color_pattern, '_all.pdf'), width=10, height=5)
                    g <- ggplot(results, aes(x=x, y=count, fill=rank))+scale_fill_viridis(discrete=TRUE)
                    g <- g+geom_bar(stat='identity', position='stack', color='black')+theme_classic()+
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                    plot(g)
                    dev.off()
                }
            }
            if (color_pattern == 'overlap') {
                selected = gene_ann[,c('ct_CD14_Monocytes', 'ct_alltypes', 'tf_alltypes', 'dar_alltypes')]
                selected = selected[gene_ann[,paste0('predictability_', time)] == 1,]
                for (i in 1:(dim(selected)[2]-1)) {
                    selected[selected[,i] > 0, (i+1):(dim(selected)[2])] = 0
                }
                selected = cbind(selected, others=unlist(sapply(rowSums(selected), function(x){return(c(1, 0)[x+1])})))
                colnames(selected) = paste0(rev(1:dim(selected)[2]), '_', colnames(selected))
                data = data.frame(group=colnames(selected), value=colSums(selected))
                # print(head(data))
                pdf(paste0('pie_chart_overlap_', time, '.pdf'), useDingbats=FALSE)
                g <- ggplot(data, aes(x="", y=value, fill=group))+geom_col()+coord_polar(theta="y")+theme_classic()
                plot(g)
                dev.off()
            }
        }
    }
}

test.annotation.overlap <- function(pred.each) {
    gene_ann = obtain.gene.annotation.tf(pred.each, monocyte_only=FALSE, all_types=TRUE, convert2factor=FALSE)
    gdf = NULL
    odf = NULL
    for (col1 in colnames(gene_ann)) {
        for (col2 in colnames(gene_ann)) {
            if (col1 == col2) next
            header = paste0(col1, '_', col2)
            temp = gene_ann
            print(c(col1, col2))
            filtered = data.matrix(temp[!is.na(temp[,col1]) & !is.na(temp[,col2]),c(col1,col2)])
            overlap = length(which(rowSums(filtered) == 2))
            group1 = length(which(filtered[,1] > 0))
            group2 = length(which(filtered[,2] > 0))
            total = dim(filtered)[1]
            print(c(overlap, group1, group2, total))
            pvalue = phyper(max(0, overlap-1), group2, total-group2, group1, lower.tail= FALSE)
            pfisher = fisher.test(matrix(c(overlap, group2-overlap, group1-overlap, total-group2-group1+overlap), 2, 2), alternative='greater')
            # if (pfisher$p.value < 0.05)
                plot.pie.chart(col1, col2, group1, overlap, header, pvalue)
            odf = rbind(odf, c(overlap/total, overlap, group1, group2, total, pvalue, pfisher$p.value))
            gdf = rbind(gdf, c(header, col1, col2))
        }
    }
    output = data.frame(gdf, odf)
    print(head(output))
    print(dim(output))
    colnames(output) = c('pair', 'pair1', 'pair2', 'ratio', 'intersect', 'group1', 'group2', 'total', 'pvalue', 'fpvalue')
    write.table(output, file='gene_annotation_enrichment.tsv', sep="\t")
}


all_cell_type_marker = FALSE
trim = TRUE

pred.each = readRDS('prediction_each_comb.rds')
pred.each = cbind(pred.each, readRDS('prediction_each_sc.rds'))
pred.each <- pred.each[!is.na(rowSums(pred.each)),]

for (fig_iter in 1:6) {
    break
    if (fig_iter == 1) {
        plot.heatmap.marker.genes(as.matrix(pred.each), 'pred', id=TRUE)
    } else if (fig_iter == 2) {
        sc_fname = '../data/bulk/data_exp_sc_rank.rds'
        df = readRDS(sc_fname)
        rank_mat = df[['test']][['X.tt3']]
        rank_mat = rank_mat[rownames(pred.each),]
        visualize.annotation.overlap(pred.each, rank_mat[,'16-902_3'])
        rank_mat = do.call(cbind, list(df[['test']][['X.tt1']], df[['test']][['X.tt2']], df[['test']][['X.tt3']]))
        rank_mat = rank_mat[rownames(pred.each),]
        visualize.annotation.overlap.whole.data(pred.each, rank_mat)
    } else {
        tf = 1
        gt = 1
        if (FALSE) {
            gene_ann = obtain.gene.annotation(pred.each, tf, gt, all_cell_type_marker)
            tail = ''
        } else {
            gene_ann = obtain.gene.annotation.tf(pred.each, df[['test']][['X.tt3']])
            tail = '_tf10'
        }
        print(head(gene_ann))
        bulk_mat <- data.frame(readRDS('../data/bulk/bulk_exp.rds'))
        bulk_mat <- change.colnames(bulk_mat)

        if (trim) {
            tail = paste0(tail, '_filt')
        } else {
            tail = tail
        }

        if (trim) {
            print(f.zz)
            bulk_mat = bulk_mat[f.zz,]
            pred.each = readRDS('prediction_each_comb.rds')
            pred.each = cbind(pred.each, readRDS('prediction_each_sc.rds'))
            pred.each = pred.each[f.zz,]
        }

        rank_test <- readRDS('../data/bulk/bulk_exp_rank_test.rds')
        rank_test = data.frame(rank_test)
        colnames(rank_test) = unlist(sapply(colnames(rank_test), function(x) {
            return(substr(x, 2, nchar(x)))
        }))
        rank_train <- readRDS('../data/bulk/bulk_exp_rank_train.rds')
        rank_train = data.frame(rank_train)
        colnames(rank_train) = unlist(sapply(colnames(rank_train), function(x) {
            return(substr(x, 2, nchar(x)))
        }))

        if (trim) {
            print('filter')
            print(dim(rank_test))
            print(dim(rank_train))
            rank_test = rank_test[f.zz,]
            rank_train = rank_train[f.zz,]
            print(dim(rank_test))
            print(dim(rank_train))
            # debug.rank(temp, rank_train) 
        }

        rank_combined <- combine.rank.data(rank_test, rank_train)

        if (fig_iter == 3) {
            bulk_mat = filt.mat(bulk_mat)
            # plot.heatmap.marker.genes(bulk_mat, paste0('exp', tail), id=TRUE, rank_flag=FALSE, gene_ann=NULL, marker_list=c('all'))
            summarize.predictable.genes(pred.each, gene_ann)
            summarize.predictable.genes.recursive(pred.each, gene_ann, tail)
        } else if (fig_iter == 4) {
                cor_bulk = cor(bulk_mat, method='spearman')
                cor_bulk = cor_bulk[order(colnames(cor_bulk)),order(colnames(cor_bulk))]
                print('before')
                plot.hist.similarity.within(cor_bulk)
                print(cor_bulk[1:10, 1:10])
                time_order <- as.numeric(unlist(sapply(1:5, function(x) {
                    vec <- c()
                    for (i in 1:3) {
                        print(c((0:3)*3+i+(x-1)*12))
                        vec <- c(vec, c((0:3)*3+i+(x-1)*12))
                    }
                    return(vec)
                })))
                print(time_order)
                cor_bulk <- cor_bulk[time_order, time_order]
                plot.heatmap.marker.genes(cor_bulk, paste0('exp_cor', tail), id=TRUE, rank_flag=TRUE, gene_ann=NULL, marker_list=c('all'), gene_flag=FALSE)
        } else if (fig_iter == 5) {
            cindex = grep('16.90._', colnames(bulk_mat))
            # plot.heatmap.marker.genes(bulk_mat[,cindex], paste0('exp_169', tail), id=TRUE, rank_flag=FALSE, gene_ann=gene_ann, marker_list=c('annotation'))
            # plot.heatmap.marker.genes(bulk_mat[,cindex], paste0('exp_169', tail), id=TRUE, rank_flag=FALSE, gene_ann=gene_ann, marker_list=c('all'))
            # plot.heatmap.marker.genes(bulk_mat[,cindex], paste0('exp_169', tail), id=TRUE, gene_ann=gene_ann, marker_list=c('marker'))
            cindex = grep(paste0('16.9_._[1-3]'), colnames(rank_test))
            plot.heatmap.marker.genes(rank_test[,cindex], paste0('exp_rank_169_raw', tail), id=TRUE, rank_flag=TRUE, gene_ann=gene_ann, marker_list=c('annotation'))
            plot.heatmap.marker.genes(rank_test[,cindex], paste0('exp_rank_169_raw', tail), id=TRUE, rank_flag=TRUE, gene_ann=gene_ann, marker_list=c('all'))
            # cindex = grep(paste0('16.9_._[1-3]'), colnames(rank_combined))
            # plot.heatmap.marker.genes(rank_combined[,cindex], paste0('exp_rank_169', tail), id=TRUE, rank=TRUE, gene_ann=gene_ann, marker_list=c('annotation'))
            # plot.heatmap.marker.genes(rank_combined[,cindex], paste0('exp_rank_169', tail), id=TRUE, rank=TRUE, gene_ann=NULL, marker_list=c('all'))
            # plot.heatmap.marker.genes(rank_combined, paste0('exp_rank_comb', tail), id=TRUE, rank=TRUE, gene_ann=NULL, marker_list=c('all'))
            # for (i in 1:3) {
            #     cindex = grep(paste0('16.9_._', i), colnames(rank_combined))
            #     plot.heatmap.marker.genes(rank_combined[,cindex], paste0('exp_rank_169_', i), id=TRUE, rank=TRUE, skip=TRUE, gene_ann=gene_ann)
            # }
            # cindex = grep(paste0('16.9_._[1-3]'), colnames(rank_combined))
            # plot.heatmap.marker.genes(rank_combined[,cindex], paste0('exp_rank_169', tail), id=TRUE, rank=TRUE, gene_ann=gene_ann, marker_list=c('annotation'))

        } else if (fig_iter == 6) {
            sc_fname = 'data_exp_sc_rank.rds'
            df = readRDS(sc_fname)
            sc_mat = df[['test']][['X.tt4']]
            colnames(sc_mat) = paste0('16.9_', 1:4, '_4')
            rank_test = cbind(rank_test, sc_mat)    
            cindex = grep(paste0('16.9_._[1-4]$'), colnames(rank_test))
            # plot.heatmap.marker.genes(rank_test[,cindex], paste0('exp_rank_169_sc', tail), id=TRUE, rank_flag=TRUE, gene_ann=gene_ann, marker_list=c('annotation'))
            # plot.heatmap.marker.genes(rank_test[,cindex], paste0('exp_rank_169_sc', tail), id=TRUE, rank_flag=TRUE, gene_ann=NULL, marker_list=c('all'))
            cor_bulk = cor(rank_test[,cindex], method='spearman')
            print(rownames(cor_bulk))
            print(colnames(cor_bulk))
            cor_bulk = cor_bulk[order(colnames(cor_bulk)),order(colnames(cor_bulk))]
            cor_bulk = cor_bulk[grep('16.9_._4$', colnames(cor_bulk)),!grepl('16.9_._4$', colnames(cor_bulk))]
            # print(cor_bulk[1:10, 1:10])
            plot.heatmap.marker.genes(cor_bulk, paste0('exp_rank_169_sc_cor', tail), id=TRUE, rank_flag=FALSE, gene_ann=NULL, marker_list=c('all'), gene_flag=FALSE)
        }
    }
}



