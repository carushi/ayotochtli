library(SnapATAC)
library(Matrix)
require(GenomicRanges)
require(ggplot2)
require(ggsci)
library(GGally)
source("../scripts/sc_clustering/multi_sc_base.R")
cov_threshold <<- 500
threshold <<- NA
filter_label <<- 'UQ'
# genes <<- NULL

genes.gr <<- read.gtf.file(FALSE)

get.promoter.bins <- function(query) {
    ov <- countOverlaps(query, genes.gr, maxgap=-1L, minoverlap=0L, type="any", ignore.strand=TRUE)
    return(unlist(sapply(ov, function(x){return(x > 0)})))
}

filter.rows <- function(x1.sp, index) {
    x.sp = x1.sp[index,];
    return(x.sp)
}

plot.qc.stats <- function(data, header='') {
    for (s in unique(data[,"sample"])) {
        temp <- subset(data, data[,'sample'] == s)
        temp <- temp[,2:6]
        g <- ggpairs(temp, columns = 1:ncol(temp), title = "",  
            axisLabels = "show")+theme_bw()
        pdf(paste0('pairplot_', s, header, '.pdf'), width=10, height=10)
        print(g)
        dev.off()
    }
    data[,'sample'] <- factor(data[,'sample'])
    temp <- NULL
    if (!is.na(threshold)) {
        temp <- data[data[,filter_label] >= threshold,]
    }
    lwd = 2.5
    for (col in colnames(data)) {
        if (any(col == c('barcode', 'sample', 'batch'))) next
        if (!any(col == c('cov', 'chrM_ratio', 'duplicate'))) {
            data[,col] <- log10(data[,col]+1)
            if (!is.null(temp)) temp[,col] <- log10(temp[,col]+1)
        }
        g <- ggplot(data, aes_string(x=col, color='sample'))+geom_freqpoly(aes(group=sample), size=lwd)+scale_color_npg()+theme_bw()
        xlim_vec <- quantile(data[,col],  probs = c(0.1, 99.9)/100, na.rm=TRUE)
        g <- g+xlim(xlim_vec)
        if (col == filter_label && !is.na(threshold))
            g <- g+geom_vline(aes_string(xintercept=threshold), col='grey', linewidth=lwd)
        pdf(paste0('histplot_', col, header, '.pdf'))
        plot(g)
        dev.off()
        if (!is.na(threshold) && header == '') {
            g <- ggplot()+geom_freqpoly(data=data, aes_string(color='sample', group='sample', x=col, y='after_stat(density)'), size=lwd-1)+geom_freqpoly(data=temp, aes_string(x=col, color='sample', group='sample', y='after_stat(density)'), linetype='dashed', size=lwd)
            g <- g+scale_color_npg()+theme_bw()
            g <- g+xlim(xlim_vec)
            pdf(paste0('histplot_', col, '_filt_uq=', threshold, '.pdf'))
            plot(g)
            dev.off()
        }
    }
}

plot.after.integration <- function(x.sp, header='') {
    bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
    pdf(paste0('bin_coverage', header, '.pdf'))
    g <- hist(
        bin.cov[bin.cov > 0], 
        xlab="log10(Bin cov)", 
        main="log10(Bin Cov)", 
        col="lightblue", 
        xlim=c(0, 5)
    );
    dev.off()
    data <- x.sp@metaData
    data[,'UM'] <- log10(data[,'UM']+1)
    pdf(paste0('umi_promoter', header, '.pdf'))
    p1 = ggplot(x.sp@metaData, aes(x=UM, y=promoter_ratio, color=sample)) + 
            geom_point(size=0.3, col="grey", alpha=0.2) +
            theme_classic()	+
            labs(x = "log UMI", y="promoter ratio")
    plot(p1)
    dev.off()
    df <- data.frame(table(x.sp@metaData[,'sample']))
    print(df)
    pdf(paste0('cell_number', header, '.pdf'))
    g <- ggplot(df, aes(x=Var1, y=Freq, fill=Var1))+geom_bar(stat='identity')+theme_classic()+scale_color_npg()
    plot(g)
    dev.off()
}

plot.before.integration <- function(snap.files) {
    df <- NULL
    for (s in names(snap.files)) {
        vec <- c(s, substr(snap.files[[s]]@metaData[1,'batch'],5,5), dim(snap.files[[s]]@bmat))
        stopifnot(length(vec) == 4)
        if (is.null(df)) df <- vec
        else df <- rbind(df, vec)
    }
    colnames(df) <- c("name", "sample", "cell", "bin")
    df <- data.frame(df)
    for (col in c('cell', 'bin')) {
        g <- ggplot(df, aes_string(x='sample', y=col, color='sample'))+geom_point()+geom_boxplot()+theme_classic()+scale_color_npg()
        pdf(paste0('plot_hist_', col, '.pdf'))
        plot(g)
        dev.off()
    }
}

plot.metadata <- function() {
    x.sp <- readRDS('snap_obj.rds')
    plot.qc.stats(x.sp@metaData)
    plot.after.integration(x.sp)
    if (is.na(threshold)) return()
    save.filtered.data(x.sp, 'snap_obj')
    x.sp <- readRDS('snap_obj_filtered.rds')
    plot.qc.stats(x.sp@metaData, '_filtered')
    plot.after.integration(x.sp, '_filtered')
}



read.snap.files <- function(fname, sample, list_flag=FALSE) {
    x1.sp <- createSnap(file=fname, sample=sample)
    if (list_flag) {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
        x1.sp = makeBinary(x1.sp, mat="bmat")
        print(dim(genes.gr))
        print(length(genes.gr))
        if (!is.null(genes.gr)) {
            print(head(genes.gr))
            index = which(Matrix::rowSums(x1.sp@bmat) > 0)
            x1.sp <- filter.rows(x1.sp, index)
            x1.sp <- createGmatFromMat(obj=x1.sp, input.mat="bmat", genes=genes.gr, do.par=TRUE, num.cores=5)
        }
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=5000)
    } else {
        x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    }
    return(x1.sp)
}

add.meta.info <- function(x.sp) {
    print('Adding meta data')
    x.sp@metaData$sample = unlist(sapply(x.sp@sample, function(x){return(as.integer(substr(x, 2,2)))}))
    x.sp@metaData$batch = x.sp@sample
    x.sp@metaData$cov = log10(Matrix::rowSums(x.sp@bmat)+1)
    x.sp@metaData$duplicate <- (1 - (x.sp@metaData[,'UQ']+1)/(x.sp@metaData[,'PP']+1))
    x.sp@metaData$chrM_ratio = x.sp@metaData[,'CM']/x.sp@metaData[,'UQ']
    col_bin <- get.promoter.bins(x.sp@feature)
    x.sp@metaData$promoter_ratio = (Matrix::rowSums(x.sp@bmat[,col_bin])+1)/(Matrix::rowSums(x.sp@bmat)+1)*100
    x.sp@metaData$UMI_ratio = (x.sp@metaData[,'UM']+1)/(x.sp@metaData[,'TN']+1)*100
    return(x.sp)
}

save.filtered.data <- function(x.sp, output) {
    cat('original dimension ')
    print(dim(x.sp@bmat))
    print(c('-> coverage cutoff', cov_threshold))
    idx = which(Matrix::rowSums(x.sp@bmat) >= cov_threshold)
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

make.snap.list <- function(ann, output='snap_obj') {
    x1.sp.ls <- list()
    for (i in 1:dim(ann)[1]) {
        print(ann[i,])
        x1.sp <- read.snap.files(ann[i,1], ann[i,2], list_flag=TRUE)
        x1.sp@metaData[,'batch'] <- ann[i,3]
        x1.sp.ls[ann[i,3]] <- x1.sp
    }
    saveRDS(x1.sp.ls, file='snap_list_obj.rds')
    plot.before.integration(x1.sp.ls)
    x.sp = Reduce(snapRbind, x1.sp.ls);
    x.sp = add.meta.info(x.sp)
    saveRDS(x.sp, file=paste0(output, '.rds'))
    if (!is.na(threshold)) save.filtered.data(x.sp, output)
}

convert.snap.files <- function(ann) {
    global_samples <- 0
    for (i in 1:dim(ann)[1]) {
        x1.sp <- read.snap.files(ann[i,1], ann[i,2], list_flag=TRUE)
        batch <- ann[i,3]
        header <- batch
        write.table(data.frame(barcode=paste0(batch, ".", x1.sp@barcode)), paste0("barcodes_", header, "_1000.tsv"), quote=F, col.names=F, row.names=F, sep=" ")
        df <- data.frame(x1.sp@feature)
        write.table(data.frame(x1.sp@feature), paste0(output_dir, "bin_", header, "_1000.tsv"), quote=F, col.names=F, row.names=F, sep=" ")
        data = x1.sp@metaData
        data[,"local_index"] = 1:dim(data)[1]
        data[,"global_index"] = 1:dim(data)[1]+global_samples
        data[,"batch"] = batch
        write.table(data, paste0("qc_", header, "_1000.tsv"), quote=F, sep=" ")
        writeMM(x1.sp@bmat, file=paste0("sparse_mat_", header, "_1000.mtx"))
        global_samples <- global_samples+dim(data)[1]
    }
}

args <- commandArgs(trailingOnly=TRUE)
ann <- read.table(args[2], header=T, sep=" ", stringsAsFactor=F)
cov_threshold <<- as.numeric(args[3])
threshold <<- as.numeric(args[4])
if (length(args) >= 5) filter_label <<- args[5]
if (args[1] == 'make_list') {
    make.snap.list(ann)
} else if (args[1] == 'convert_text') {
    convert.snap.files(ann)
} else {
    plot.metadata()
}


