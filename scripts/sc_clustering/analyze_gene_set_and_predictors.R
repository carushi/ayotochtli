require(parallel)

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
source("../scripts/sc_clustering/multi_sc_base.R")
source("../scripts/annotation/gene_enrichment.R")
source("../scripts/multi_omics/bulk_pred_utils.R")
source("../scripts/multi_omics/sc_pred_utils.R")
load(paste0(DIR, "Xscaffolds.Rdata"))
scaffoldsX = scaffoldsX.prop3
scaffoldsX = scaffoldsX.sub


quad_names = get.quad.names(X.cpm.all)

pred.each = readRDS('prediction_each_comb.rds')
pred.each = cbind(pred.each, readRDS('prediction_each_sc.rds'))


dname = './'
files = list.files(path=dname, pattern='^rank_*')
for (file in files) {
    break
    df = read.table(file.path(dname, file), header=T, sep="\t", stringsAsFactor=FALSE)
    mat = matrix(0, nrow=dim(pred.each)[1], ncol=length(unique(df[,'cluster'])))
    rownames(mat) = rownames(pred.each)
    colnames(mat) = sort(unique(df[,'cluster']))
    for (col in colnames(mat)) {
        gene_list = df[df[,'cluster'] == col,'gene']
        mat[as.vector(sapply(gene_list, function(x) which(rownames(mat) == x)[1])),col] = 1
    }
    prefix = gsub('.tsv', '', gsub('rank_gene_', '', file), fixed=TRUE)
    plot.heatmap.binary.two.mat(pred.each, mat, paste0('predictability_', prefix), col_fun=c('darkgray', 'lightgray', 'red'))
    write.table(mat, file=paste0(prefix, '_mat.tsv'), sep="\t")
    # plot.corrplot(mat, paste0('predictability_', prefix))
}


dname = '../data/gene/'
files = list.files(path=dname, pattern='*.gmt.txt')
for (file in files) {
    break
    print(file.path(dname, file))
    df = read.table(file.path(dname, file), header=F, sep="\t", stringsAsFactor=FALSE, fill=NA, quote="\"", col.names=1:7000)
    prefix = gsub('.gmt.txt', '', file)
    fname = paste0(prefix, '_mat.tsv')
    if (file.exists(fname) && FALSE) {
        mat = read.table(fname, sep="\t", header=T)
        print(head(mat))
    } else {
        mat = matrix(0, nrow=dim(pred.each)[1], ncol=dim(df)[1])
        rownames(mat) = rownames(pred.each)
        colnames(mat) = sort(df[,1])
        rownames(df) = df[,1]
        df = df[,-c(1)]
        symbol = id2symbol(genes, id=rownames(mat))
        for (col in colnames(mat)) {
            print(col)
            gene_list = df[rownames(df) == col,]
            gene_list = gene_list[!is.na(gene_list)]
            mat_index = as.vector(unlist(sapply(gene_list, function(x) {
                index = which(symbol == x)
                if (length(index) > 0) return(index[1])
            })))
            mat[mat_index,col] = 1
        }
        print(dim(mat))
        mat = mat[,colSums(mat) > 0]
        print(dim(mat))
        write.table(mat, file=fname, sep="\t")
    }
}
files = list.files(path=dname, pattern='*._mat.tsv')
for (file in files) {
    print(file.path(dname, file))
    mat = read.table(file.path(dname, file), header=F, sep="\t", stringsAsFactor=FALSE, fill=NA, quote="\"", col.names=1:7000)
    print(head(mat))
    # marker.set.enrichment(pred.each, data.frame(mat), paste0(prefix), category=25, height=12)
    plot.heatmap.binary.two.mat(pred.each, mat, paste0('predictability_', prefix), col_fun=c('darkgray', 'lightgray', 'red'))
}
exit()
dname = './'
window = 10000
for (cluster in c('cluster', 'rna_celltype', 'batch')) {
    # break
    file_list = paste0(dname, 'peaks_', window, '/', cluster, '_list.txt')
    fname_list = read.table(file_list, header=F, stringsAsFactor=F)[,1]
    mat = matrix(0, nrow=dim(pred.each)[1], ncol=length(fname_list))
    rownames(mat) = rownames(pred.each)
    colnames(mat) = unlist(gsub('.bed', '', fname_list))
    fname = paste0('atac_', cluster, '_mat.tsv')
    if (file.exists(fname) && FALSE) {
        mat = read.table(fname, sep="\t", header=T)
        print(head(mat))
    } else {
        for (fname in fname_list) {
            prefix = gsub('.bed', '', fname)
            df = read.table(file.path(dname, 'peaks_10000', fname), header=F, stringsAsFactor=F)
            colnames(df) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'gchr', 'gstart', 'gend', 'gstrand', 'gene_id', 'symbol')
            gene_list = df[,'gene_id']
            mat_index = as.vector(unlist(sapply(gene_list, function(x) {
                index = which(rownames(mat) == x)
                if (length(index) > 0) return(index[1])
            })))
            mat[mat_index,prefix] = 1
            gene.set.enrichment(data.frame(gene_id=gene_list), data.frame(gene_id=rownames(mat)), paste0('atac_', cluster, '_', prefix))
        }
    }
    # plot.heatmap.binary.two.mat(pred.each, data.frame(mat), paste0('predictability_atac_', cluster), col_fun=c('darkgray', 'lightgray', 'red'))
    marker.set.enrichment(pred.each, data.frame(mat), paste0('atac_', cluster), category=10, height=5)
    write.table(mat, file=paste0('atac_', cluster, '_peak_mat.tsv'), sep="\t")
}


dname = './'
thres = 0.01
for (cluster in c('cluster', 'rna_celltype', 'batch')) {
    file_list = paste0(dname, 'DAR_', thres, '/DAR_', cluster, '_list.txt')
    fname_list = read.table(file_list, header=F, stringsAsFactor=F)[,1]
    mat = matrix(0, nrow=dim(pred.each)[1], ncol=length(fname_list))
    rownames(mat) = rownames(pred.each)
    colnames(mat) = unlist(gsub('.bed', '', fname_list))
    fname = paste0('dar_', cluster, '_mat.tsv')
    if (file.exists(fname) && FALSE) {
        mat = read.table(fname, sep="\t", header=T)
        print(head(mat))
    } else {
        for (fname in fname_list) {
            prefix = gsub('.bed', '', fname)
            print(prefix)
            df = read.table(file.path(dname, paste0('DAR_', thres), fname), header=F, stringsAsFactor=F)
            colnames(df) = c('chr', 'start', 'end', 'width', 'strand', 'peak', 'gchr', 'gstart', 'gend', 'gstrand', 'gene_id', 'symbol')
            gene_list = df[,'gene_id']
            mat_index = as.vector(unlist(sapply(gene_list, function(x) {
                index = which(rownames(mat) == x)
                if (length(index) > 0) return(index[1])
            })))
            mat[mat_index,prefix] = 1
            gene.set.enrichment(data.frame(gene_id=unique(gene_list)), data.frame(gene_id=rownames(mat)), paste0('dar_', cluster, '_', prefix))
        }
    }
    marker.set.enrichment(pred.each, data.frame(mat), paste0('dar_', cluster), category=10, height=5)
    write.table(mat, file=paste0('dar_', cluster, '_peak_mat.tsv'), sep="\t")
}