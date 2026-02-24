require(Seurat)
require(SnapATAC)
require(ggplot2)
require(stringr)
set.seed(42)
source("../scripts/sc_clustering/multi_sc_base.R")
dir = "./" # seurat scRNA-seq data without refseq CD4
dir_ref = "./" # seurat scRNA-seq data with refseq CD4

add.cd4.from.refseq <- function(rna, remove=FALSE) {
    new_rna = readRDS(paste0(dir_ref, "seurat_scrna_integrated.rds"))
    cd4_vec = new_rna@assays[['RNA']][c('CD4'),]
    en_id = 'ENSDNOG00000002205'
    if (remove) {
        index = ((rownames(rna) != en_id) & (rownames(rna) != "CD4"))
        count_mat = rbind(rna@assays[['RNA']]@counts[index, ], cd4_vec)
    } else {
        index = (rownames(rna) != "CD4")
        count_mat = as.matrix(rna@assays[['RNA']]@counts[index,])
        print(dim(count_mat))
        rownames(count_mat)[which(rownames(count_mat) == en_id)[1]] = paste0(en_id, '.en')
        count_mat = rbind(count_mat, cd4_vec)
    }
    rownames(count_mat)[dim(count_mat)[1]] = en_id
    se_obj <- CreateSeuratObject(counts=count_mat, assay="RNA")
    for (col in colnames(rna@meta.data)) {
        se_obj <- AddMetaData(object = se_obj, metadata = rna@meta.data[,col], col.name=col)
    }
    return(se_obj)    
}


convert.to.seurat <- function(x.sp, remove=FALSE) {
    mat = x.sp@gmat
    colnames(mat) = genes[,5]
    if (remove) {
        en_id = 'ENSDNOG00000002205'
        mat = mat[,colnames(mat) != en_id]
        colnames(mat)[colnames(mat) == 'CD4'] = en_id
    } else {
        en_id = 'ENSDNOG00000002205'
        colnames(mat)[colnames(mat) == en_id] = paste0(en_id, '.en')
        colnames(mat)[colnames(mat) == 'CD4'] = en_id
    }
    se_obj <- CreateSeuratObject(counts=t(mat), assay="RNA")
    for (cluster in c('cluster')) {
        print(cluster)
        se_obj <- AddMetaData(object = se_obj, metadata = factor(x.sp@metaData[,cluster]), col.name=cluster)
    }
    se_obj <- AddMetaData(object = se_obj, metadata = factor(unlist(sapply(x.sp@metaData[,'sample'], function(x){return(x-1)}))), col.name='batch')
    for (column in c('tsne', 'umap')) {
        for (i in 1:2) {
            se_obj <- AddMetaData(object = se_obj, metadata = slot(x.sp, column)[,i], col.name=paste0(column, i))
        }
    }
    return(se_obj)    
}

CD4 = TRUE

# rna = readRDS(file.path(dir, paste0('final_pseudo_bulk_', 'rna', '.rds')))
# if (CD4) {
#     rna = add.cd4.from.refseq(rna)
#     atac_snap = readRDS('../snap_211207/snap_obj_clust.rds')
#     atac = convert.to.seurat(atac_snap)
# } else {
#     atac_snap = readRDS('snap_obj_clust.rds')
#     atac = convert.to.seurat(atac_snap)
# }

rna = readRDS(file.path(dir, paste0('seurat_scrna_integrated.rds')))
rna = add.cd4.from.refseq(rna, remove=FALSE) # Original ENSDNOG00000002205 count is maintained as ENSDNOG00000002205.en if remove == FALSE
atac_snap = readRDS('snap_obj_clust.rds')
atac = convert.to.seurat(atac_snap, remove=FALSE)

rna <- add.additional.metadata(rna, TRUE)
rna <- manual.curation.rna.celltype(rna)
rna <- normalize.seurat(rna, 'rna', reclustering=FALSE)

saveRDS(rna, 'final_pseudo_bulk_rna.rds')


atac = RenameAssays(object=atac, RNA='ACTIVITY')
VariableFeatures(atac) <- row.names(atac)
atac <- RunLSI(atac, n=50, scale.max=NULL)

transfer.anchors <- FindTransferAnchors(reference=rna, query=atac, features=VariableFeatures(object=rna), 
    reference.assay="RNA", query.assay="ACTIVITY", reduction="cca")

genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype,
    weight.reduction = atac[["lsi"]], dims = 2:30)
atac <- AddMetaData(atac, metadata= celltype.predictions)
atac <- AddMetaData(atac, metadata=celltype.predictions[,'predicted.id'], col.name='rna_celltype')
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]], dims = 2:30)
atac[["RNA"]] <- imputation

saveRDS(atac, file='final_coembed_atac.rds')

colnames(rna@meta.data)[which(colnames(rna@meta.data) == 'mcluster')] = 'cluster'
rna <- AddMetaData(rna, metadata=rna@meta.data[,'celltype'], col.name='rna_celltype')

coembed <- merge(x = rna, y = atac)
coembed <- normalize.seurat(coembed, 'both', reclustering=FALSE)
coembed@meta.data[,'cluster'] = str_pad(coembed@meta.data[,'cluster'], 2, pad='0')
coembed@meta.data[,'orig.ident'] = unlist(sapply(coembed@meta.data[,'orig.ident'], function(x) {
    if (x == 'SeuratProject') return('ATAC')
    else return('RNA')
}))
saveRDS(coembed, 'final_coembed.rds')


anchors <- FindIntegrationAnchors(
    object.list = list(rna, atac),
    anchor.features = genes.use,
    assay = c('RNA', 'ACTIVITY'),
    k.filter = NA
)
integrated <- IntegrateData(
    anchorset = anchors,
    weight.reduction = atac[['lsi']],
    dims = 2:30,
    preserve.order = TRUE
)

saveRDS(integrated, 'first_integrated.rds')

integrated <- scale.seurat(integrated)

integrated@meta.data[,'cluster'] = str_pad(integrated@meta.data[,'cluster'], 2, pad='0')
integrated@meta.data[,'orig.ident'] = unlist(sapply(integrated@meta.data[,'orig.ident'], function(x) {
    if (x == 'SeuratProject') return('ATAC')
    else return('RNA')
}))

saveRDS(integrated, 'final_integrated.rds')

hex_code_list = set.hex.code(integrated)
tail = 'integrated'
for (cluster in c('cluster', 'batch', 'rna_celltype')) {
    tail = 'both'
    for (dr in c('umap', 'tsne', 'pca')) {
        pdf(paste0(dr, '_merged_', tail, '_', cluster, '.pdf'), width=10, height=7)
        plot(DimPlot(integrated, reduction=dr, cols=hex_code_list[[cluster]], group.by=cluster))
        dev.off()
        pdf(paste0(dr, '_merged_', tail, '_', cluster, '_batch.pdf'), width=24, height=7)
        plot(DimPlot(integrated, reduction=dr, cols=hex_code_list[[cluster]], split.by='batch', group.by=cluster))
        dev.off()
    }
}
