require(GenomicRanges)
require(scales)
require(stringr)


GTF_DIR <<- "../data/genome/"
# Select gene definition
GTF_FILE <<- paste0(GTF_DIR, "16-90_transcript_gb.gtf")
GTF_FILE <<- paste0(GTF_DIR, "genome_cd4_transcript_gb.gtf")
GTF_FILE <<- paste0(GTF_DIR, "16-90_transcript_tss.gtf")

MARKER_DIR <<- "../data/gene/"

MARKER <<- c("marker_list1.txt", "marker_list2.txt", "marker_list3.txt", "marker_integrated.txt", "marker_list5.txt", "marker_integrated2.txt", "marker_integrated_sm.txt",
             "cd_gene.txt", "il_gene.txt", "ccr_gene.txt", "marker_egr_bob1.txt", "marker_egr_bob2.txt", "marker_egr_bob3.txt")
MODULE <<- c('module_1.txt', 'module_2.txt', 'module_pred_16-9_2_1.txt', 'module_pred_16-9_2_2.txt', 'module_pred_16-9_2_3.txt', 'module_pred_16-9_2_4_1.txt', 'module_pred_16-9_2_4_2.txt', 'module_pred_16-9_2_4_3.txt')

MITO_GENE <<- "../scripts/annotation/chrM_gene.txt"

BIN_DIR <<- "~/miniforge/env/genome/bin/"


change.colnames <- function(mat) {
    colnames(mat) = unlist(sapply(colnames(mat), function(x) {
        return(substr(x, 2, nchar(x)))
    }))
    return(mat)
}

colMax <- function(X) apply(X, 2, max)

gene.region.file <- function(region='gb', genome='genome_cd4') {
    # region = c('tss', 'gb')
    # file_list = paste0(c('16-90', 'genome_cd4', 'genome', 'genome_refseq'), '_transcript_', region, '.gtf')
    file = paste0(genome, '_transcript_', region, '.gtf')
    return(file.path(GTF_DIR, file))
}

filt.mat <- function(mat) {
    mat <- mat[,Matrix::colSums(mat) > 0]
    min_vec = unlist(apply(mat, c(2), function(x) min(x, na.rm=TRUE)))
    max_vec = unlist(apply(mat, c(2), function(x) max(x, na.rm=TRUE)))
    print(dim(mat))
    print(length(min_vec))
    print(length(max_vec))
    mat <- mat[, unlist(sapply(1:length(min_vec), function(x) {return(min_vec[x] != max_vec[x])}))]
    return(mat)
}

add.cocluster <- function(cluster, rna=TRUE) {
    if (rna) {
        original_cluster <- (unlist(sapply(cluster, function(x) { # 0-origin in python
            return(c(1, 2, 2, 3, 4, 5, 6, 7, 8, 7, 4, 7, 9, 4, 10, 1)[as.integer(as.character(x))+1])
        })))
    } else {
        original_cluster <- unlist(sapply(cluster, function(x) { # 1-origin 
            return(c(8, 3, 1, 1, 5, 2, 3, 5, 7, 6, 9, 4, 10)[as.integer(as.character(x))])
        }))
    }
    ordered_cluster <- unlist(sapply(original_cluster, function(x){return(c(4, 10, 1, 2, 5, 9, 6, 3, 7, 8)[x])}))
    print(c(original_cluster[1:5], ordered_cluster[1:5]))
    return(ordered_cluster)
}

add.mcluster <- function(cluster, umap1, rna=TRUE) {
    stopifnot(rna)
    mcluster <- unlist(sapply(1:length(cluster), function(x) { # 0-origin in python
        if (cluster[x] != 14) return(cluster[x])
        if (umap1[x] < 0) return(16)
        else return(14)
    }))
    return(mcluster)
}

add.mcocluster <- function(cluster, umap1, rna=TRUE) {
    print(unique(cluster))
    if (rna) {
        original_cluster <- (unlist(sapply(cluster, function(x) { # 0-origin in python
            return(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 10, 11, 12, 13, 14, 15, 16)[as.integer(as.character(x))+1])
            })))
    } else {
        original_cluster <- unlist(sapply(cluster, function(x) { # 1-origin 
            return(c(6, 4, 4, 15, 1, 3, 4, 4, 8, 1, 12, 5, 3)[as.integer(as.character(x))])
        }))
        print(unique(original_cluster))
    }
    ordered_cluster <- original_cluster
    ord <- c(4, 6, 9, 1, 5, 10, 8, 13, 14, 7, 16, 12, 11, 2, 3, 15)
    ord_cluster <- unlist(sapply(seq(1, 16), function(x) {which(ord == x)}))
    ordered_cluster <- unlist(sapply(original_cluster, function(x){return(ord_cluster[x])}))
    return(ordered_cluster)
}

add.celltype <- function(cluster, umap1, rna=TRUE) {
    stopifnot(rna)
    print(unique(cluster))
    celltype_dict <- c('CD4+ T and Tregs', 'CD4+ T', 'CD8+ T', 'Naive T', 'Naive T', 'Naive T', 'B', 'Plasma B', 
    'DC', 'NK', 'HSC', 'CD14+ Monocytes', 'Monocytes', 'Monocytes', 'Monocytes', 'Undefined')
    majortype_dict <- c('T', 'T', 'T', 'T', 'T', 'T', 'B', 'B', 'DC', 'NK', 'HSC', 'M', 'M', 'M', 'M', 'No signal')
    celltype <- unlist(sapply(cluster, function(x){return(celltype_dict[x])}))
    majortype <- unlist(sapply(cluster, function(x){return(majortype_dict[x])}))
    return(cbind(celltype=celltype, majortype=majortype))
}

celltype.dict <- function(cluster) {
    if (cluster == 'celltype')
        celltype_dict <- c('CD4+ T and Tregs', 'CD4+ T', 'CD8+ T', 'Naive T', 'B', 'Plasma B', 'DC', 'NK', 'HSC', 'CD14+ Monocytes', 'Monocytes', 'Undefined')
    else
        celltype_dict <- c('T', 'B', 'DC', 'NK', 'HSC', 'M', 'No signal')
    return(celltype_dict)
}

manual.curation.rna.celltype <- function(rna, input_label='mcluster', output_label='celltype', debug=FALSE) {
    require(stringr)
    print(table(rna@meta.data[,input_label]))
    rna@meta.data[,output_label] = unlist(sapply(rna@meta.data[,input_label], function(x) {
        # cluster = c('10', '04', '08', '13', '06', '05', '00', '03', '01', '02', '11', '07', '09', '12', '15', '14') 
        celltype = list('00'='CD4+ T', '01'='CD14+ Monocytes', '02'='CD14+ Monocytes', '03'='CD4+ Treg', '04'='CD8+ Treg', '05'='CD4+ naive T', '06'='NK',  '07'='B', '08'='CD8+ T',  '09'='B', '10'='DN T',
        '11'='Monocytes',  '12'='Plasma B',  '13'='Monocytes + T',  '14'='DC', '15'='Undefined', '16'='HSC')
        if (is.na(x)) return(NA)
        return(celltype[[str_pad(x, 2, pad='0')]])
    }))
    print(table(rna@meta.data[,output_label]))
    if (debug) {
        gmat = t(rna@assays[['RNA']]@data)
        cluster = rna@meta.data[,output_label]
        plot.marker.mean(gmat, cluster, '_nrna')
    }
    return(rna)
}

add.additional.metadata <- function(x.sp, rna) {
    if (rna) {
        x.sp@meta.data[['cocluster']] <- add.cocluster(x.sp@meta.data[['cluster_leiden']], rna=rna)
        x.sp@meta.data[['mcluster']] <- add.mcluster(x.sp@meta.data[['cluster_leiden']], x.sp@meta.data[['umap1']], rna=rna)
        x.sp@meta.data[['mcocluster']] <- add.mcocluster(x.sp@meta.data[['mcluster']], rna=rna)
        celltype <- add.celltype(x.sp@meta.data[['mcocluster']])
        x.sp@meta.data[['celltype']] <- celltype[,1]
        x.sp@meta.data[['majortype']] <- celltype[,2]
    } else {
        x.sp@metaData <- cbind(x.sp@metaData, cocluster=add.cocluster(x.sp@cluster, rna=rna))
        x.sp@metaData <- cbind(x.sp@metaData, mcocluster=add.mcocluster(x.sp@cluster, rna=rna))
        celltype <- add.celltype(x.sp@metaData[,'mcocluster'])
        x.sp@metaData <- cbind(x.sp@metaData, celltype)
    }
    return(x.sp)
}

colPanel <<- c(
    "grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
    "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
    "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
) # SnapATAC color scheme

read.gtf.file <- function(matrix=TRUE) {
    gtf_name <- GTF_FILE
    gtf <- read.table(gtf_name, header=F, sep="\t")
    colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'name')
    gtf[,'name'] <- apply(gtf, c(1), function(x){if(is.na(x['name']))return(x['gene_id']); return(x['name'])})
    gtf <- gtf[!duplicated(gtf),]
    if (matrix) return(gtf)
    genes.gr <- GRanges(gtf[,1], IRanges(gtf[,2], gtf[,3]), strand=gtf[,4], gene_id=gtf[,5], name=gtf[,6]);
    return(genes.gr)
}

id2symbol <- function(genes, symbol=NULL, id=NULL) {
    if (!is.null(id) || !is.null(symbol)) {
        if (is.null(id)) {
            from = 6
            to = 5
            vec = symbol
        } else {
            from = 5
            to = 6
            vec = id
        }
        new_vec = unlist(sapply(vec, function(x) {
            index = which(genes[,from] == x)
            return(genes[which(genes[,from] == x)[1], to])
        }))
        return(new_vec)
    }
    return(c())
}

convert.to.symbol <- function(vec, rev=FALSE) {
    gtf_name <- GTF_FILE
    gtf <- read.table(gtf_name, header=F, sep="\t", stringsAsFactor=FALSE)
    colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'name')
    from = 5
    to = 6
    if (rev) {
        from = 6
        to = 5
    }
    symbols <- sapply(as.character(vec), function(x) {
        ind <- which(gtf[,from] == x)
        if (length(ind) == 0) return(x)
        else {
            if (rev) {
                ind <- ind[unlist(sapply(gtf[ind,to], function(x) {
                    return(grepl('ENSD', x, fixed=TRUE))
                }))]
            }
            return(gtf[ind[1],to])
        }
    })
    symbols <- as.character(symbols)
    return(symbols)      
}

convert.gene.id <- function(marker_genes) {
    marker_genes <- as.character(marker_genes)
    for (i in 1:length(marker_genes)) {
        if (marker_genes[i] == "IL32") {
            marker_genes[i] = "RF00003"
        }
        if (marker_genes[i] == "CD8B") {
            marker_genes[i] = "RF00017"
        }
        if (marker_genes[i] == "CST3") {
            marker_genes[i] = "STOML2"
        }
    }
    return(marker_genes)
}

gene.region.data <- function() {
    gtf_name <- GTF_FILE
    gtf <- read.table(gtf_name, header=F, sep="\t", stringsAsFactor=FALSE)
    colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'symbol')
    gtf <- gtf[!duplicated(gtf),]
    gtf[,'symbol'] <- apply(gtf, c(1), function(x){if(is.na(x['symbol']))return(x['gene_id']); return(x['symbol'])})
    return(gtf)
}

extract.all.marker.genes <- function(files) {
	all_marker_genes <- NULL
    for (fname in files) {
        marker.genes <- unlist(read.table(file.path(mdir, fname), header=F)[,1])
        marker.genes <- convert.gene.id(marker.genes)
		all_marker_genes <- rbind(all_marker_genes, marker.genes)
	}
    return(all_marker_genes)
}


plot.marker.genes <- function(x.sp, marker_file='', symbol=FALSE) {
    if (marker_file == '') {
	    files <- MARKER
    } else {
        files <- c(marker_file)
    }
	marker_genes <- extract.all.marker.genes(files)
	column_vec <- unlist(sapply(colnames(x.sp@gmat), function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
	for (dim_reduct in c('umap', 'tsne')) {
		for (m in unique(marker.genes[,'symbol'])) {
			if (symbol)  {
				j = which(m == column_vec)
			} else {
				j = unlist(sapply(marker.genes[marker.genes[,'symbol'] == m,'gene_id'], function(x){which(column_vec == x)}))
			}
			if (length(j) == 0) next
			pdf(paste0('marker_', dim_reduct, '_', m, '.pdf'))
				plotGene(
					obj=x.sp,
					gene=colnames(x.sp@gmat)[j],
					viz.method=dim_reduct, 
					point.size=0.1, 
					point.shape=19, 
					down.sample=10000,
					quantiles=c(0.01, 0.99),
				)
				dev.off()
		}
	}
}


plot.marker.mean <- function(gmat, cluster, tail, marker_file='') {
    require(ComplexHeatmap)
    require(dplyr)
    require(viridis)
    require(BBmisc)
    require(circlize)
    require(Matrix)
    if (marker_file == '') {
        files <- MARKER
    } else {
        files <- c(marker_file)
    }
    column_vec <- unlist(sapply(colnames(gmat), function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
    column_vec <- convert.to.symbol(column_vec)

    for (i in 1:length(files)) {
        # if (i == 5 && grepl('imputed', tail))  break # memory error
        file <- files[i]
        marker.genes <- extract.all.marker.genes(c(file))
        mat <- do.call(cbind, lapply(marker.genes, function(x) {
            index = which(column_vec == x)
            if (length(index) == 0)
                return(rep(0, dim(gmat)[1]))
            else if (length(index) == 1) {
                return(gmat[,index])
            } else {
                return(colMax(gmat[,index]))
            }
        }))
        mat <- data.frame(mat)
        colnames(mat) = marker.genes
        mat <- mat[,Matrix::colSums(mat) > 0]
        mat <- cbind(mat, cluster=as.character(cluster))
        data <- mat %>% group_by(cluster) %>% summarise_all(mean, na.rm = TRUE)
        print(head(data))
        cluster_names = unlist(data[,1])
        mat <- data.matrix(as.data.frame(data[,-c(1)]))
        print(dim(mat))
        print(length(cluster_names))
        print(cluster_names)
        rownames(mat) <- as.vector(cluster_names)
        for (norm in c('range', 'standardize')) {
            normed_mat <- as.matrix(t(BBmisc::normalize(t(mat), method=norm, range=c(0, 1))))
            colnames(normed_mat) = colnames(mat)
            rownames(normed_mat) = rownames(mat)
            col_fun = colorRamp2(seq(0, 1, 0.01), viridis(101))
            pdf(paste0('heatmp_mean_norm_', i, '_', norm, tail, '.pdf'))
            draw(Heatmap(normed_mat, col=col_fun))
            dev.off()
            pdf(paste0('heatmp_mean_norm_', i, '_', norm, tail, '_ordered.pdf'))
            draw(Heatmap(normed_mat, col=col_fun, cluster_rows = FALSE, cluster_columns=FALSE))
            dev.off()
        }
        pdf(paste0('heamap_mean_', i, tail, '.pdf'))
        draw(Heatmap(mat, col=viridis(100), show_row_names = TRUE))
        dev.off()
    }
}

scale.seurat <- function(x.sp) {
    x.sp <- ScaleData(object = x.sp)
    x.sp <- RunPCA(object = x.sp)
    x.sp <- FindNeighbors(object = x.sp)
    x.sp <- RunTSNE(object = x.sp, dims=1:15, check_duplicates = FALSE)
    x.sp <- RunUMAP(object = x.sp, dims=1:15)
    return(x.sp)
}

normalize.seurat <- function(x.sp, tail, reclustering=FALSE) {    
    x.sp <- subset(x.sp, features=rownames(x.sp@assays[["RNA"]])[!grepl("^gSpikein-", rownames(x.sp@assays[["RNA"]]))])
    chrM_gene <- as.vector(read.table(MITO_GENE, header=F, stringsAsFactor=F))
    x.sp <- subset(x.sp, features=rownames(x.sp@assays[["RNA"]])[!(rownames(x.sp@assays[["RNA"]]) %in% chrM_gene)])
    x.sp <- NormalizeData(object = x.sp)
    x.sp <- FindVariableFeatures(object = x.sp, selection.method = "vst", nfeatures=5000)
    x.sp <- ScaleData(object = x.sp)
    x.sp <- RunPCA(object = x.sp)
    x.sp <- FindNeighbors(object = x.sp)
    if (reclustering) {
        x.sp <- FindClusters(object = x.sp)
    }
    x.sp <- RunTSNE(object = x.sp, dims=1:15, check_duplicates = FALSE)
    x.sp <- RunUMAP(object = x.sp, dims=1:15)
    x.sp <- JackStraw(x.sp, num.replicate = 100)
    x.sp <- ScoreJackStraw(x.sp, dims = 1:20)    
    pdf(paste0('pca_elbow_', tail, '.pdf'))
    plot(JackStrawPlot(x.sp, dims = 1:20))
    dev.off()
    return(x.sp)
}



strsplit.bracket <- function(seq) {
    seq = gsub('[', '"', seq, fixed=TRUE)
    seq = gsub(']', '"', seq, fixed=TRUE)
    return(gsub(', ', ' ', seq))
}

read.module.list <- function(sp, x.sp) {
    marker_file = file.path(MARKER_DIR, paste0('cell_marker_', sp, '_blood.tsv'))
    df = read.table(marker_file, header=T, sep="\t", quote="\"", stringsAsFactor=F)
    ng_word = c('progen', 'eryth', 'Baso', 'eosino', 'Megaka')
    for (word in ng_word) {
        df = df[!grepl(word, df[,'celltype'], ignore.case=TRUE), ]
    }
    marker_list = list()
    for (i in 1:dim(df)[1]) {
        symbol = strsplit.bracket(df[i,'symbol'])
        id = strsplit.bracket(df[i, 'gene_id'])
        symbols = unlist(strsplit(gsub('"', '', symbol), ' '))
        ids = c(unlist(strsplit(id, ',')))
        marker_list[[as.character(paste0(df[i,'celltype'], '_', length(ids), '_', i))]] = unlist(sapply(ids, function(x) { if(any(x == rownames(x.sp))) return(x) }))
    }
    for (n in names(marker_list)) {
        x.sp <- AddModuleScore(
            x.sp,
            list(marker_list[[n]]),
            pool = NULL,
            nbin = 24,
            ctrl = 100,
            k = FALSE,
            assay = NULL,
            name = gsub(' ', '_', paste0("Module_", n), fixed=TRUE),
            seed = 1,
            search = FALSE,
        )
    }
    return(x.sp)
}

set.hex.code <- function(x.sp) {
    require(viridis)
    columns = colnames(x.sp@meta.data)
    hex_code_list = list()
    for (col in c('cluster', 'rna_celltype', 'batch')) {
        if (!any(col == columns)) next
        if (col == 'batch')
            hex_codes1 <- viridis_pal()(length(unique(x.sp@meta.data[,col])))
        else
            hex_codes1 <- rev(hue_pal()(length(unique(x.sp@meta.data[,col]))))
        names(hex_codes1) = sort(unique(x.sp@meta.data[,col]))
        hex_code_list[[col]] = hex_codes1
    }
    return(hex_code_list)
}

genes <<- read.gtf.file()