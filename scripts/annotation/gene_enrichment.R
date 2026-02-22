require(GO.db)
require(GOfuncR)
require(clusterProfiler)
require(enrichR)
require(parallel)
require(ggplot2)
load("../bulk_transcriptome/data/GO.arm.Rdata")
SCRIPT_DIR <<- "../scripts/"
source(paste0(SCRIPT_DIR, "sc_clustering/multi_sc_base.R"))


make.go.mat <- function() {
    fname = "go.all.mat.tsv"
    if (file.exists(fname)) {
        return(read.table(fname, header=T, sep="\t"))
    }
    xx <- as.list(GOTERM)
    mat <- do.call(rbind, lapply(1:length(xx), function(x){return(c(GOID(xx[[x]]), Term(xx[[x]]), Ontology(xx[[x]])))}))        
    write.table(mat, 'go.all.mat.tsv', sep="\t")
    return(mat)
}

go_term <<- make.go.mat()
colnames(go_term) <- c("diseaseId", "diseaseName", "ontology")
go_term <<- subset(go_term, go_term[,"ontology"] == "BP")
go_mat <<- GO.arm[,c(2,3)]
term2gene <<- go_mat[,c(2,1)]
colnames(go_mat) <- c("gene", "go_id")
colnames(term2gene) <- c("diseaseId", "geneId")
colnames(genes.id) <- c("symbol", "entrez")
go_mat_non <<- subset(go_mat, GO.arm[,4] != 'IEA')
term2gene_non <<- subset(term2gene, GO.arm[,4] != 'IEA')

add.id.annotation <- function(gene_set) {
    gene_set <- gene_set[unlist(sapply(gene_set[,'gene_id'], function(x){grepl("^ENSDNOG", x)})),,drop=F]
    if (!any(colnames(gene_set) == 'symbol')) {
        m <- merge(gene_set, genes, by='gene_id', all.x=TRUE)
        colnames(m)[colnames(m) == 'name'] = 'symbol'
        gene_set = m
    }
    gene_set <- merge(gene_set, genes.id, all.x=TRUE, by='symbol')
    return(gene_set)
}

marker.set.enrichment <- function(pred.each, marker.set, header, category=10, height=4) {
    require(tidyr)
    require(ggplot2)
    require(dplyr)
    marker_term = data.frame(set_id=paste0('GO:', 1:dim(marker.set)[2]), name=colnames(marker.set))
    colnames(marker.set) = paste0('GO:', 1:dim(marker.set)[2])
    marker2gene = gather(add_rownames(marker.set, var='gene'), "marker", "exist", colnames(marker.set))
    marker2gene = data.frame(marker2gene[marker2gene[,3] > 0,])
    marker2gene = marker2gene[,c(2, 1)]
    print(head(marker2gene))
    marker2gene[,2] = unlist(sapply(marker2gene[,2], function(x) {
        return(which(x == rownames(pred.each))[1])
    }))
    write.table(rownames(pred.each), file=paste0(header, '_marker_gene_list.tsv'), sep="\t")
    print(head(marker.set))
    rownames(pred.each) = 1:dim(pred.each)[1]
    mclapply(1:dim(pred.each)[2], function(i) {
        print(i)
        background = rownames(pred.each)[!is.na(pred.each[,i])]
        gene_set = rownames(pred.each[background,])[pred.each[background,i] > 0]
        thres = dim(gene_set)[1]
        prefix = paste0(header, '_marker_', colnames(pred.each)[i], '_', thres)
        tryCatch({
            x = enricher(as.integer(gene_set), as.integer(background), TERM2GENE=marker2gene, TERM2NAME=marker_term)
            # x = enricher(as.character(1:length(gene_set)), as.character(1:length(background)), TERM2GENE=term2gene, TERM2NAME=go_term)
            write.table(x@result, file=paste0(prefix, ".tsv"), sep="\t")
            pdf(paste0('dotplot_', prefix, ".pdf"), width=10, height=height, useDingbats=F)
            plot(dotplot(x, showCategory=category) + scale_y_discrete(labels=function(x) str_wrap(x, width=40)))
            dev.off()
            pdf(paste0('barplot_', prefix, ".pdf"), width=10, height=height, useDingbats=F)
            plot(barplot(x, showCategory=category) + scale_y_discrete(labels=function(x) str_wrap(x, width=40)))
            dev.off()
        }, error=function(e) print(e))
    }, mc.cores=7)
}


gene.set.enrichment <- function(gene_set, background, header) {
    gene_set = add.id.annotation(gene_set)
    background = add.id.annotation(background)
    write.table(gene_set, file=paste0('gene_table_', header, '.tsv'), sep="\t")
    gene_set = gene_set[!is.na(gene_set[,'entrez']),]
    background = background[!is.na(background[,'entrez']),]
    print(head(gene_set))
    print(head(background))
    for (filt in c('', '_nonIEA')) {
        for (thres in c(dim(gene_set)[1])) {
            prefix = paste0(header, '_', thres, '_cp', filt)
            # gene <- input_hyper_bg[input_hyper_bg[,2] >= 10000-thres,1]
            if (filt == '')
                x = enricher(gene_set[,'entrez'], background[,'entrez'], TERM2GENE=term2gene, TERM2NAME=go_term)
            else 
                x = enricher(gene_set[,'entrez'], background[,'entrez'], TERM2GENE=term2gene_non, TERM2NAME=go_term)
            tryCatch({
                write.table(x@result, file=paste0(prefix, ".tsv"), sep="\t")
                x@result = x@result
                pdf(paste0('dotplot_', prefix, ".pdf"), width=10, height=4, useDingbats=F)
                plot(dotplot(x, showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=40)))
                dev.off()
                pdf(paste0('barplot_', prefix, ".pdf"), width=10, height=4, useDingbats=F)
                plot(barplot(x, showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=40)))
                dev.off()
            }, error=function(e) print(e))
        }
    }
}
    
