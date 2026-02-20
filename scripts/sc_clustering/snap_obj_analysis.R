require(SnapATAC)
require(parallel)
library(GGally)
library(GenomicRanges)
source("../scripts/sc_clustering/multi_sc_base.R")

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
    return(x.sp)
}

runMitLineage <- function(
	obj, 
	output.prefix,
	path.to.snaptools,
	path.to.macs,
	gsize,
	tmp.folder,
	buffer.size=500,
	num.cores=1,
	keep.minimal=TRUE
) {
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
		if((x=length(obj@barcode))==0L) {
			stop("obj@barcode is empty");			
		}
		if((x=length(obj@file))==0L){
			stop("obj@file is empty");			
		}
	}
	
	fileList = as.list(unique(obj@file));
	
	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if FM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following nsap files do not contain FM session")
		print(fileList[idx])
		stop()
	}
	
	if(missing(output.prefix)){
		stop("output.prefix is missing");
	}
	
	if(missing(path.to.snaptools)){
		stop("path.to.snaptools is missing");
	}else{
		if(!file.exists(path.to.snaptools)){
			stop("path.to.snaptools does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.snaptools);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.snaptools is not an excutable file");
		}
	}

	if(missing(path.to.macs)){
		stop("path.to.macs is missing");
	}else{
		if(!file.exists(path.to.macs)){
			stop("path.to.macs does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.macs);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.macs is not an excutable file");
		}
	}
	
	if(missing(gsize)){
		stop("gsize is missing");
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
		
	# write the following barcodes down
	barcode.files = lapply(fileList, function(file){
		tempfile(tmpdir = tmp.folder, fileext = ".barcode.txt");
	})

	bed.files = lapply(fileList, function(file){
		tempfile(tmpdir = tmp.folder, fileext = ".bed.gz");
	})
	
	# write down the barcodes
	cat("Epoch: extracting fragments from each snap files ...\n", file = stderr())
	flag.list = lapply(seq(fileList), function(i){
		file.name = fileList[[i]];
		idx = which(obj@file == file.name);
		barcode.use = obj@barcode[idx]
		write.table(barcode.use, file = barcode.files[[i]], append = FALSE, quote = FALSE, sep = "\t",
		                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
		                 col.names = FALSE, qmethod = c("escape", "double"),
		                 fileEncoding = "")
		
	})
	
	# extract the fragments belong to the barcodes	
	flag.list = mclapply(seq(fileList), function(i){
		flag = system2(command=path.to.snaptools, 
			args=c("dump-fragment", 
				   "--snap-file", fileList[[i]], 
				   "--output-file", bed.files[[i]], 
				   "--barcode-file", barcode.files[[i]],
				   "--buffer-size", buffer.size
				   )		
			)				
	}, mc.cores=num.cores);
	
	# combine these two bed files
	combined.bed = tempfile(tmpdir = tmp.folder, fileext = ".bed.gz");
	flag = system2(command="cat", 
		args=c(paste(bed.files, collapse = ' '),
			   ">", combined.bed
			   )		
		)				
	
	# call peaks using MACS2	
	flag = system2(command=path.to.macs, 
		args=c("callpeak", 
			   "-t", combined.bed, 
			   "-f", "BED",
			   "-g", gsize,
			   macs.options,
			   "-n", output.prefix
			   )		
		)				
	if (flag != 0) {
	   	stop("'MACS' call failed");
	}	

	if(keep.minimal){
		system(paste("rm ", output.prefix, "_control_lambda.bdg", sep=""));
		system(paste("rm ", output.prefix, "_peaks.xls", sep=""));
		system(paste("rm ", output.prefix, "_summits.bed", sep=""));
	}

	return(read.table(paste(output.prefix, "_peaks.narrowPeak", sep="")));
}


#' Call Peaks Using MACS2 For All Clusters
#'
#' Identify peaks for all clusters. Fragments belonging to each subset or cluster of cells 
#' are extracted and used to identify peaks using MACS2. This function requires
#' "MACS2" and "snaptools" preinstalled and excutable. 
#' 
#' @param obj A snap object.
#' @param output.prefix Prefix of output file which will be used to generate output file names.
#' @param path.to.snaptools Path to snaptools excutable file.
#' @param path.to.macs Path to macs2 excutable file.
#' @param min.cells min number of cells to perform peak calling [100]. Clusters with cells less 
#' than num.cells will be excluded.
#' @param num.cores number of cpus to use [1].
#' @param gsize effective genome size. 'hs' for human, 'mm' for mouse, 'ce' for C. elegans, 'dm' for fruitfly (default: None)
#' @param buffer.size Buffer size for incrementally increasing internal array size to store reads alignment information. In
#' most cases, you don't have to change this parameter. However, if there are very high coverage dataset that each barcode has
#' more than 10000 fragments, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read snap files) [1000].
#' @param macs.options String indicate options you would like passed to macs2. strongly do not recommand to change unless you know what you are doing. the default is '--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits'.
#' @param tmp.folder Directory to store temporary files generated by runMACSForAll.
#' @param keep.minimal Keep minimal version of output [TRUE].
#' @return Return a GRanges object that contains the non-overlapping combined peaks
#' @importFrom parallel mclapply 
#' @importFrom GenomicRanges GRanges reduce
#' @export
#' 
runMACS <- function(
	obj, 
	output.prefix,
	path.to.snaptools,
	path.to.macs,
	gsize,
	tmp.folder,
	buffer.size=500,
	macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
	num.cores=1,
	keep.minimal=TRUE
){
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
		if((x=length(obj@barcode))==0L){
			stop("obj@barcode is empty");			
		}
		if((x=length(obj@file))==0L){
			stop("obj@file is empty");			
		}
	}
	
	fileList = as.list(unique(obj@file));
	
	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if FM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following nsap files do not contain FM session")
		print(fileList[idx])
		stop()
	}
	
	if(missing(output.prefix)){
		stop("output.prefix is missing");
	}
	
	if(missing(path.to.snaptools)){
		stop("path.to.snaptools is missing");
	}else{
		if(!file.exists(path.to.snaptools)){
			stop("path.to.snaptools does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.snaptools);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.snaptools is not an excutable file");
		}
	}

	if(missing(path.to.macs)){
		stop("path.to.macs is missing");
	}else{
		if(!file.exists(path.to.macs)){
			stop("path.to.macs does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.macs);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.macs is not an excutable file");
		}
	}
	
	if(missing(gsize)){
		stop("gsize is missing");
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
		
	# write the following barcodes down
	barcode.files = lapply(fileList, function(file){
		tempfile(tmpdir = tmp.folder, fileext = ".barcode.txt");
	})

	bed.files = lapply(fileList, function(file){
		tempfile(tmpdir = tmp.folder, fileext = ".bed.gz");
	})
	
	# write down the barcodes
	cat("Epoch: extracting fragments from each snap files ...\n", file = stderr())
	flag.list = lapply(seq(fileList), function(i){
		file.name = fileList[[i]];
		idx = which(obj@file == file.name);
		barcode.use = obj@barcode[idx]
		write.table(barcode.use, file = barcode.files[[i]], append = FALSE, quote = FALSE, sep = "\t",
		                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
		                 col.names = FALSE, qmethod = c("escape", "double"),
		                 fileEncoding = "")
		
	})
	
	# extract the fragments belong to the barcodes	
	flag.list = mclapply(seq(fileList), function(i){
		flag = system2(command=path.to.snaptools, 
			args=c("dump-fragment", 
				   "--snap-file", fileList[[i]], 
				   "--output-file", bed.files[[i]], 
				   "--barcode-file", barcode.files[[i]],
				   "--buffer-size", buffer.size
				   )		
			)				
	}, mc.cores=num.cores);
	
	# combine these two bed files
	combined.bed = tempfile(tmpdir = tmp.folder, fileext = ".bed.gz");
	flag = system2(command="cat", 
		args=c(paste(bed.files, collapse = ' '),
			   ">", combined.bed
			   )		
		)				
	
	# call peaks using MACS2	
	flag = system2(command=path.to.macs, 
		args=c("callpeak", 
			   "-t", combined.bed, 
			   "-f", "BED",
			   "-g", gsize,
			   macs.options,
			   "-n", output.prefix
			   )		
		)				
	if (flag != 0) {
	   	stop("'MACS' call failed");
	}	

	if(keep.minimal){
		system(paste("rm ", output.prefix, "_control_lambda.bdg", sep=""));
		system(paste("rm ", output.prefix, "_peaks.xls", sep=""));
		system(paste("rm ", output.prefix, "_summits.bed", sep=""));
	}

	return(read.table(paste(output.prefix, "_peaks.narrowPeak", sep="")));
}



runMACSForAll <- function(
	obj, 
	output.prefix,
	path.to.snaptools,
	path.to.macs,
	gsize,
	tmp.folder,
	num.cores=1,
	min.cells=100,
	buffer.size=500,
	macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
	keep.minimal=TRUE
) {
	cat("Epoch: checking input parameters ... \n", file = stderr())
		
    peak.ls <-  parallel::mclapply(as.list(levels(obj@cluster)), function(x){
       num.cells = length(which(obj@cluster == x));
       if(num.cells < min.cells){
           return(GenomicRanges::GRanges())
       }
	   peaks.df <- runMACS(
		   obj=obj[which(obj@cluster==x),], 
		   output.prefix=paste(output.prefix, x, sep="."),
		   path.to.snaptools=path.to.snaptools,
		   path.to.macs=path.to.macs,
		   gsize=gsize,
		   buffer.size=buffer.size, 
		   macs.options=macs.options,
		   tmp.folder=tmp.folder,
		   keep.minimal=keep.minimal,
		   num.cores=1
		   );
	   if((x=nrow(peaks.df)) > 0L){
	       peaks.gr = GenomicRanges::GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]));
	   }else{
	   	   peaks.gr = GenomicRanges::GRanges();
	   }
     }, mc.cores=num.cores)
	
	peaks.gr = GenomicRanges::reduce(do.call(c, peak.ls));
	return(peaks.gr);
}

filter.computed.cluster <- function(clusters, header) {
	filtered <- c()
	for (i in clusters) {
		if (!file.exists(paste0(header, ".", gsub(" ", "_", i), "_peaks.narrowPeak"))) {
			filtered <- c(filtered, i)
		}
	}
	return(filtered)
}

combine.peaks <- function(header, clusters) {
    peaks.names = system(paste0("ls | grep narrowPeak | grep ", header), intern=TRUE);
    peak.gr.ls = lapply(peaks.names, function(x){
        peak.df = read.table(x)
        GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
    })
    peak.gr = reduce(Reduce(c, peak.gr.ls));
    peaks.df = as.data.frame(peak.gr)[,1:3];
    write.table(peaks.df, file = paste0("peaks_", header, ".combined.bed"),append=FALSE,
		quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
		row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
		fileEncoding = "")
}

cluster.specific.peaks <- function(x.sp, cluster='cluster') {
    system("which snaptools");
    system("which macs2");
	header = cluster
    clusters.all = names(table(x.sp@metaData[,cluster]))[which(table(x.sp@metaData[,cluster]) > 10)]
	clusters.sel = filter.computed.cluster(clusters.all, header)
    peaks.ls = mclapply(seq(clusters.sel), function(i){
        print(clusters.sel[i])
        peaks = runMACS(
            obj=x.sp[which(x.sp@metaData[,cluster]==clusters.sel[i]),], 
            output.prefix=paste0(header, ".", gsub(" ", "_", clusters.sel)[i]),
            path.to.snaptools=paste0(BIN_DIR, "snaptools"),
            path.to.macs=paste0(BIN_DIR, "macs2"),
            gsize='3.17e9', # mm, hs, etc
            buffer.size=500, 
            num.cores=1,
            macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
            tmp.folder="./"
       );
 	}, mc.cores=5)
	combine.peaks(header, clusters.all)
}

batch.specific.peaks <- function(x.sp) {
    clusters.all = unique(x.sp@metaData[,'batch'])
	header = 'batch'
	clusters.sel = filter.computed.cluster(clusters.all, header)
    peaks.ls = mclapply(clusters.sel, function(i){
        peaks = runMACS(
            obj=x.sp[which(x.sp@metaData[,'batch']==i),], 
            output.prefix=paste0(header, ".", gsub(" ", "_", i)),
            path.to.snaptools=paste0(BIN_DIR, "snaptools"),
            path.to.macs=paste0(BIN_DIR, "macs2"),
            gsize='3.17e9', # mm, hs, etc
            buffer.size=500, 
            num.cores=1,
            # macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
            macs.options="--nomodel --shift 100 --ext 200 --qval 0.1 -B --SPMR",
            tmp.folder="./"
       );
 	}, mc.cores=5)
	combine.peaks(header, clusters.all)
}

analyze.pmat <- function(x.sp, cluster_label, bcv=0.4) {
	print(levels(x.sp@cluster))
	x.sp@cluster = factor(x.sp@metaData[,cluster_label])
	print(levels(x.sp@cluster))

	for (cluster.pos in levels(x.sp@cluster)) {
		print(cluster.pos)
		print(unique(x.sp@cluster))
		if (is.na(cluster.pos)) next
		print((!any(cluster.pos %in% levels(x.sp@cluster))))
		print(levels(x.sp@cluster))
		if (!any(cluster.pos == x.sp@cluster)) next
		print(table(x.sp@cluster))
		print('all cells')
		print(length(which(x.sp@cluster == cluster.pos)))
		DARs = findDAR(
			obj=x.sp,
			input.mat="pmat",
			cluster.pos=cluster.pos,
			cluster.neg.method="random",
			test.method="exactTest",
			bcv=bcv, #0.4 for human, 0.1 for mouse
			seed.use=10
		);
		DARs$FDR = p.adjust(DARs$PValue, method="BH");
		idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
		print(c('detected peaks', length))
		if (length(idy) <= 1) next
		covs = Matrix::rowSums(x.sp@pmat);
		vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
		vals.zscore = (vals - mean(vals)) / sd(vals);
		pdf(paste0('plot_feature_', cluster, '_', cluster.pos, '.pdf'))
		plotFeatureSingle(
			obj=x.sp,
			feature.value=vals.zscore,
			method="umap", 
			main=cluster.pos,
			point.size=0.1, 
			point.shape=19, 
			down.sample=5000,
			quantiles=c(0.01, 0.99)
		);
		dev.off()
		pdf(paste0('DAR_', cluster, '_', cluster.pos, '.pdf'))
		par(mfrow = c(1, 2));
		plot(DARs$logCPM, DARs$logFC, 
			pch=19, cex=0.1, col="grey", 
			ylab="logFC", xlab="logCPM",
			main=cluster.pos
		);
		points(DARs$logCPM[idy], 
			DARs$logFC[idy], 
			pch=19, 
			cex=0.5, 
			col="red"
		);
		abline(h = 0, lwd=1, lty=2);
		dev.off()
 		write.table(data.frame(x.sp[,idy,"pmat"]@peak), file=paste0('DAR_', cluster, '_', gsub(' ', '_', cluster.pos), '.tsv'), sep="\t", quote=F, row.names=F)
		next
		for (type in c('known', 'denovo')) {
			if (type == 'known') {
				only.known = TRUE
				only.denovo = FALSE
			} else {
				only.known = FALSE
				only.denovo = TRUE
			}
			out <- tryCatch(
				{			
					motifs = runHomer(
						x.sp[,idy,"pmat"], 
						mat = "pmat",
						path.to.homer = paste0(BIN_DIR, "findMotifsGenome.pl"),
						result.dir = paste0("./homer/", cluster, '_', gsub(' ', '_', cluster.pos), '_', type),
						num.cores=5,
						genome = 'dasNov3',
						motif.length = 10,
						scan.size = 300,
						optimize.count = 2,
						background = 'automatic',
						local.background = FALSE,
						only.known = only.known,
						only.denovo = only.denovo,
						fdr.num = 5,
						cache = 100,
						overwrite = TRUE,
						keep.minimal = FALSE
					)
				}, error = function(cond) {
					message(cond)
				}
    		)    
		}
	}
}

cluster_file <- "snap_obj_clust.rds"
x.sp <- readRDS(cluster_file)
x.sp <- add.additional.metadata(x.sp, FALSE)

print(x.sp)

annotated <- readRDS("final_coembed_atac.rds")
x.sp@metaData <- cbind(x.sp@metaData, rna_celltype=annotated@meta.data[,'rna_celltype'])

print(head(x.sp@metaData))


batch.specific.peaks(x.sp)
cluster.specific.peaks(x.sp, 'cluster')
cluster.specific.peaks(x.sp, 'rna_celltype')

x.sp = addPmatToSnap(x.sp)
tmp_gr = as.data.frame(x.sp@peak)
tmp_seqnames = as.character(tmp_gr[,1])
tmp_seqnames_new = str_replace( str_replace(tmp_seqnames, "b'", ""), "'", "")
tmp_gr[,1] = tmp_seqnames_new
tmp_names = as.character(tmp_gr[,6])
tmp_names_new = str_replace( str_replace(tmp_names, "b'", ""), "'", "")
tmp_gr[,6] = tmp_names_new
tmp_gr_new = makeGRangesFromDataFrame(tmp_gr, keep.extra.columns=TRUE)
x.sp@peak = tmp_gr_new

for (cluster in c('cluster', 'rna_celltype', 'batch')) {
    analyze.pmat(x.sp, cluster, 0.01)
}

