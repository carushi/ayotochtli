library(SnapATAC)
library(Matrix)
library(ggplot2)
library(ggsci)
library(GGally)
library(GenomicRanges)
source("../scripts/sc_clustering/multi_sc_base.R")


genes <- gene.region.data()

add.dimension.reduction <- function(x.sp) {
    set.seed(42);
    x.sp@metaData[,"logUMI"] = log10(x.sp@metaData[,"UM"]+1)
    row.covs.dens <- density(
        x = x.sp@metaData[,"logUMI"],
        bw = 'nrd', adjust = 1
    );
    sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps);
    idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob));    
    x.landmark.sp = x.sp[idx.landmark.ds,];
    x.query.sp = x.sp[-idx.landmark.ds,];

    x.landmark.sp = runDiffusionMaps(
        obj= x.landmark.sp,
        input.mat="bmat", 
        num.eigs=50
    );
    x.landmark.sp@metaData$landmark = 1;
    x.query.sp = runDiffusionMapsExtension(
        obj1=x.landmark.sp, 
        obj2=x.query.sp,
        input.mat="bmat"
    );
    x.query.sp@metaData$landmark = 0;
    x.sp = snapRbind(x.landmark.sp, x.query.sp);
    x.sp = x.sp[order(x.sp@metaData[,"sample"])]; #IMPORTANT
    pdf(paste0('dim_reduct.pdf'))
    g <- plotDimReductPW(
        obj=x.sp, 
        eigs.dims=1:50,
        point.size=0.3,
        point.color="grey",
        point.shape=19,
        point.alpha=0.6,
        down.sample=5000,
        pdf.file.name=NULL, 
        pdf.height=7, 
        pdf.width=7
    );
    print(g)
    dev.off()
    return(x.sp)
}

add.clustering <- function(x.sp) {
    x.sp = runKNN(
        obj=x.sp,
        eigs.dims=1:20,
        k=15
    );
    require(leiden);
    x.sp=runCluster(
        obj=x.sp,
        tmp.folder=tempdir(),
        louvain.lib="leiden",
        seed.use=10,
        resolution=1
    );
    x.sp@metaData$cluster = x.sp@cluster;
    return(x.sp)
}

plot.clustering <- function(x.sp, header='') {
    for (dr in c('tsne', 'umap')) {
        if (dr == 'tsne') {
            x.sp = runViz(
                obj=x.sp, 
                tmp.folder=tempdir(),
                dims=2,
                eigs.dims=1:20, 
                method="Rtsne",
                seed.use=10
            );
        } else {
            library(umap);
            x.sp = runViz(
                obj=x.sp, 
                tmp.folder=tempdir(),
                dims=2,
                eigs.dims=1:20, 
                method="umap",
                seed.use=42
            );
        }
        print('run vis')
        print('dr')
        print(length(x.sp@metaData$sample))

        pdf(paste0(header, dr, '_cluster.pdf'))
        plotViz(
            obj=x.sp,
            method=dr, 
            main=dr,
            point.color=x.sp@cluster, 
            point.size=1, 
            point.shape=19, 
            point.alpha=0.8, 
            text.add=TRUE,
            text.size=1.5,
            text.color="black",
            text.halo.add=TRUE,
            text.halo.color="white",
            text.halo.width=0.2,
            down.sample=10000,
            legend.add=TRUE
        );
        dev.off()

        pdf(paste0(header, dr, '_cov.pdf'))
        plotFeatureSingle(
            obj=x.sp,
            feature.value=x.sp@metaData$cov,
            method=dr, 
            main="Armadillo quadruplet cell coverage",
            point.size=0.2, 
            point.shape=19, 
            down.sample=10000
            # quantiles=c(0.01, 0.99)    
        ); 
        dev.off()
        pdf(paste0(header, dr, '_sample.pdf'))
        plotViz(
            obj= x.sp,
            method=dr,
            main="Sample",
            point.size=0.2, 
            point.shape=19,
            point.color=x.sp@metaData$sample, 
            text.add=FALSE,
            text.size=1.5,
            text.color="black",
            down.sample=10000,
            legend.add=TRUE
            );
        dev.off()
        pdf(paste0(header, dr, '_landmark.pdf'))
        plotViz(
            obj= x.sp,
            method=dr, 
            main="Landmark",
            point.size=0.2, 
            point.shape=19, 
            point.color=x.sp@metaData[,"landmark"], 
            text.add=FALSE,
            text.size=1.5,
            text.color="black",
            down.sample=10000,
            legend.add=TRUE
        );
        dev.off()

    }
    print('end')
    return(x.sp)
}

plot.violinplot <- function(x.sp, symbol=FALSE) {
    require(Seurat)
    files <- c("marker_list1.txt", "marker_list2.txt", "marker_list3.txt", "marker_integrated.txt", "marker_list5.txt", "marker_integrated2.txt", "marker_integrated_sm.txt")
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

plot.marker.genes <- function(x.sp, marker_file='', symbol=FALSE) {
    if (marker_file == '') {
	    files <- c("marker_list1.txt", "marker_list2.txt", "marker_list3.txt", "marker_integrated.txt", "marker_list5.txt", "marker_integrated2.txt", "marker_integrated_sm.txt")
    } else {
        files <- c(marker_file)
    }
	all_marker_genes <- NULL
    for (i in 1:length(files)) {
        file <- files[i]
        marker.genes <- unlist(read.table(file.path(mdir, file), header=F)[,1])
        marker.genes <- convert.gene.id(marker.genes)
		all_marker_genes <- rbind(all_marker_genes, marker.genes)
        print(marker.genes)
        print(dim(marker.genes))
	}
	marker_genes <- all_marker_genes
	column_vec <- unlist(sapply(colnames(x.sp@gmat), function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
	for (dim_reduct in c('umap', 'tsne')) {
		for (m in unique(marker.genes[,'symbol'])) {
			if (symbol)  {
				j = which(m == column_vec)
			} else {
				j = unlist(sapply(marker.genes[marker.genes[,'symbol'] == m,'gene_id'], function(x){which(column_vec == x)}))
			}
			print(c('marker genes', m))
			print(j)
			print(column_vec[1:10])
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

filter.chr <- function(x.sp, max_idx) {
    chr.exclude = seqlevels(x.sp@feature)[grep("chrM", seqlevels(x.sp@feature))];
    idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
    if (length(idy) > 0) {
        x.sp@bmat = x.sp@bmat[,-idy]
    };
    x.sp = makeBinary(x.sp, mat="bmat");
    bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
    bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
    idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
    x.sp@bmat = x.sp@bmat[, idy]
    x.sp@bmat = x.sp@bmat[,order(Matrix::colSums(x.sp@bmat), decreasing=TRUE)]
    return(x.sp)
}

args <- commandArgs(trailingOnly=TRUE)
header <- args[2]
if (is.na(header)) header <- ''
cluster_file <- "snap_obj_clust.rds"

if (!file.exists(cluster_file)) {
    x.sp <- readRDS(args[1])
    x.sp <- filter.chr(x.sp)
    print('DR')
    x.sp <- add.dimension.reduction(x.sp)
    print('clustering')
    x.sp <- add.clustering(x.sp)
    saveRDS(x.sp, 'snap_obj_clust.temp')
    x.sp <- plot.clustering(x.sp, header)
    saveRDS(x.sp, 'snap_obj_clust.rds')
} else {
    x.sp <- readRDS(cluster_file)
}

