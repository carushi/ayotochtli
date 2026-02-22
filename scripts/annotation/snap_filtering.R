cov_threshold <<- 1000 # for snap_mismatch_align
cov_threshold <<- 500  # for snap

filter_label <<- "UQ"

save.filtered.data <- function(x.sp, output) {
    cat('original dimension ')
    print(dim(x.sp@bmat))
    print(c('-> coverage cutoff', cov_threshold))
    idx = which(Matrix::rowSums(x.sp@bmat) >= cov_threshold);
    x.sp_filt <- x.sp[idx,]
    if (!is.na(threshold)) {
        print(paste('-> external cutoff', filter_label, '>=', 10**(threshold)))
        index = which(x.sp_filt@metaData[,filter_label] >= 10**(threshold))
        print(c('before', dim(x.sp_filt@metaData)[1]))
        print(c('after', length(index)))
        x.sp_filt <- x.sp_filt[index,]
    }
    saveRDS(x.sp_filt, file=paste0(output, '_filtered.rds'))
    print(dim(x.sp_filt@bmat))
    print(dim(x.sp_filt@metaData))
}


genes_gtf <<- NULL

read.gtf.file <- function(gene='gb') {
    if (gb)
        gtf_name <- "../data/genome/16-90_transcript_gb.gtf"
    else
        gtf_name <- "../data/genome/16-90_transcript_tss.gtf"
    gtf <- read.table(gtf_name, header=F, sep="\t")
    colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'name')
    gtf[,'name'] <- apply(gtf, c(1), function(x){if(is.na(x['name']))return(x['gene_id']); return(x['name'])})
    genes_gtf <<- makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE)
}

read.snap.files <- function(fname, sample, list_flag=TRUE) {
    x1.sp <- createSnap(file=fname, sample=sample)
    if (list_flag) {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
        x1.sp = makeBinary(x1.sp, mat="bmat");
        if (!is.null(genes)) {
            index = which(Matrix::rowSums(x1.sp@bmat) > 0)
            x1.sp <- filter.rows(x1.sp, index)
            x1.sp <- createGmatFromMat(obj=x1.sp, input.mat="bmat", genes=genes, do.par=TRUE, num.cores=5)
        }
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=5000)
    } else {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    }
    return(x1.sp)
}