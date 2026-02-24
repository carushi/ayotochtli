require(SnapATAC)

x.sp <- readRDS("snap_obj_clust.rds")
for (batch in unique(x.sp@metaData[,'batch'])) {
    barcode = x.sp@metaData[x.sp@metaData[,'batch'] == batch,'barcode']
    barcode = unique(barcode)
    write.table(barcode, file=paste0('barcode_atac_', batch, '.tsv'), quote=F, row.names=F, col.names=F)
    for (sample in unique(x.sp@metaData[,'sample'])) {
        barcode = x.sp@metaData[x.sp@metaData[,'batch'] == batch & x.sp@metaData[,'sample'] == sample,'barcode']
        barcode = unique(barcode)
        write.table(barcode, file=paste0('barcode_atac_', batch, '_', sample, '.tsv'), quote=F, row.names=F, col.names=F)
    }
}