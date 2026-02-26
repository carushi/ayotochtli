library(Seurat)
library(Matrix)

raw_data = c('original', 'cd4', 'refseq')[2]
if (raw_data == 'original') {
    bam_dir = 'bam'
} else if (raw_data == 'cd4') {
    bam_dir = 'bam_p_cd4'
} else {
    bam_dir = 'bam_refseq'
}
obs_mat <- read.table(file.path('scarmadillo_filtered_obs.csv'), header=T, sep=",")

print(unique(obs_mat[,'cluster_leiden']))
if (raw_data != 'refseq') {
    obs_mat <- obs_mat[,-c(1)]
}

for (i in 1:4) {
    id = i+3
    sample_id = paste0("30526", id)
    solo_dir <- file.path('../data/scRNA/bam', bam_dir, paste0(sample_id, 'Solo.out'), 'Gene', 'filtered')
    print(solo_dir)
    mat <- readMM(file.path(solo_dir, 'matrix.mtx'))
    print(dim(mat))
    barcodes <- read.table(file.path(solo_dir, 'barcodes.tsv'))
    features <- read.table(file.path(solo_dir, 'features.tsv'))
    colnames(features) <- c('gene_id', 'gene_name')
    print(dim(features))
    print(dim(subset(obs_mat, obs_mat[,'batch'] == i-1)))
    print(dim(barcodes))
    temp_mat <- obs_mat[obs_mat[,'batch'] == i-1,]
    barcodes <- cbind(barcodes, temp_mat)
    print(head(barcodes))
    barcodes <- barcodes[,-c(1)]
    print(head(obs_mat))
    print(head(features))
    rownames(features) <- features[,1]
    rownames(mat) <- rownames(features)
    colnames(mat) <- rownames(barcodes)
    obj <- CreateSeuratObject(counts = mat, project = "scrna")
    obj <- AddMetaData(object = obj, metadata = barcodes)
    obj[["RNA"]] <- AddMetaData(obj[['RNA']], metadata=features)
    saveRDS(obj, paste0('seurat_scrna_', i, '.rds'))
}

obj <- NULL
for (i in 1:4) {
    tobj <- readRDS(paste0('seurat_scrna_', i, '.rds'))
    if (is.null(obj)) obj <- tobj
    else {
        obj <- merge(obj, y = tobj, project = "scrna")
    }

}
saveRDS(obj, paste0('seurat_scrna_integrated.rds'))

    