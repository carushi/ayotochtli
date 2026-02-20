library(SnapATAC)
library(Matrix)
require(GenomicRanges)


dir <- ""
relaxed <- FALSE
args <- commandArgs(trailingOnly=TRUE)
if (length(args) >= 2) {
    input_dir <- args[1]
    output_dir <- args[2]
} else {
    if (relaxed) {
    input_dir <- "../data/snap_relaxed/"
    output_dir <- "../data/matrix_relaxed/"
    } else {
    input_dir <- "../data/snap/"
    output_dir <- "../data/matrix/"
    }
}
print(input_dir)
print(output_dir)

# gtf_name <- "../data/16-90_transcript_tss.gtf"
gtf_name <- "../data/genome_cd4_transcript_gb.gtf"

gtf <- read.table(gtf_name, header=F, sep="\t")
colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'name')
gtf[,'name'] <- apply(gtf, c(1), function(x){if(is.na(x['name']))return(x['gene_id']); return(x['name'])})
genes <- makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE)

files <- list.files(path = input_dir, pattern="*.snap$")

count <- 1
global_samples <- 0
list_flag <- TRUE # make snap object
x.sp.ls = list()
for (file in files) {
    base_name <- unlist(strsplit(gsub('.snap', '', file), '_'))[4:5]
    header <- paste(base_name, collapse='_')
    batch <- header
    if (header != 'S1_L003') next
    print(c(base_name[1], batch))
    x1.sp <- createSnap(file=file.path(input_dir, file), sample=base_name[1])

    print(x1.sp)
    print(genes)
    if (list_flag) {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=5000)
        x1.sp = makeBinary(x1.sp, mat="bmat");
        x1.sp <- createGmatFromMat(obj=x1.sp, input.mat="bmat", genes=genes, do.par=TRUE, num.cores=5)
        x.sp.ls[batch] <- x1.sp
        count <- count+1
        print(x.sp.ls)
        next
    } else {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    }
    write.table(data.frame(barcode=paste0(batch, ".", x1.sp@barcode)), paste0(output_dir, "barcodes_", header, "_1000.tsv"), quote=F, col.names=F, row.names=F, sep=" ")
    df <- data.frame(x1.sp@feature)
    print(head(df))
    write.table(data.frame(x1.sp@feature), paste0(output_dir, "bin_", header, "_1000.tsv"), quote=F, col.names=F, row.names=F, sep=" ")
    data = x1.sp@metaData
    data[,"local_index"] = 1:dim(data)[1]
    data[,"global_index"] = 1:dim(data)[1]+global_samples
    data[,"batch"] = batch
    write.table(data, paste0(output_dir, "qc_", header, "_1000.tsv"), quote=F, sep=" ")
    writeMM(x1.sp@bmat, file=paste0(output_dir, "sparse_mat_", header, "_1000.mtx"))
    count <- count+1
    global_samples <- global_samples+dim(data)[1]
}

if (!list_flag) exit()
x.sp = Reduce(snapRbind, x.sp.ls);
saveRDS(x.sp, file='snap_obj.rds')
