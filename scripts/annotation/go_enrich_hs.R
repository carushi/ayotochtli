require(GO.db)
require(GOfuncR)
require(clusterProfiler)
require(stringr)
load("../bulk_transcriptome/data/GO.arm.Rdata")
load("../bulk_transcriptome/data/homologs.Rdata")
pathways <- read.table('../data/gene/immune_pathway.txt', header=F, sep="\t", stringsAsFactor=FALSE)
gtf_name <<- "../data/genome/16-90_transcript_gb.gtf"
gtf <- read.table(gtf_name, header=F, sep="\t", stringsAsFactor=FALSE)
colnames(gtf) <- c('chr','start', 'end', 'strand', 'gene_id', 'name')
gtf <- gtf[gtf[,'chr'] != 'chrM',c('gene_id', 'name')]
print(head(gtf))



plot.KEGG.pathway <- function(map, gene.data, outfile, kegg.code, all.gene.data=NULL) {
    require(pathview)
    cpd = NULL
    kegg = T
    layer = F
    tail = "png"
    print(outfile)
    tryCatch({
        pv.out <- pathview(gene.data=gene.data, gene.idtype="entrez", pathway.id=map, species=kegg.code, out.suffix=map, keys.align="y", kegg.native=kegg, match.data=T, same.layer=layer)
        print(str(pv.out))
        plot.name <- paste0(outfile, "_", map, ".", tail)
        print(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
        system(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
        print(plot.name)
    }, error=function(cond) {
        print(cond)
        return(NA)
    })
    if (is.null(all.gene.data)) return()
    tryCatch({
        pv.out <- pathview(gene.data=all.gene.data, gene.idtype="entrez", pathway.id=map, species=kegg.code, out.suffix=map, keys.align="y", kegg.native=kegg, match.data=T, same.layer=layer)
        print(str(pv.out))
        plot.name <- paste0(outfile, "_", map, "_all.", tail)
        print(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
        system(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
        print(plot.name)
    }, error=function(cond) {
        print(cond)
        return(NA)
    })
}


make.go.mat <- function() {
    fname = "go.all.mat.tsv"
    if (file.exists(fname)) {
        return(read.table(fname, header=T, sep="\t", stringsAsFactor=FALSE))
    }
    xx <- as.list(GOTERM)
    mat <- do.call(rbind, lapply(1:length(xx), function(x){return(c(GOID(xx[[x]]), Term(xx[[x]]), Ontology(xx[[x]])))}))        
    write.table(mat, 'go.all.mat.tsv', sep="\t")
    return(mat)
}

construct.go.data <- function() {
    go_term <- make.go.mat()
    colnames(go_term) <- c("diseaseId", "diseaseName", "ontology")
    go_term <- subset(go_term, go_term[,"ontology"] == "BP")
    go_mat <- GO.arm[,c(2,3)]
    colnames(go_mat) <- c("gene", "go_id")
    term2gene <- go_mat[,c(2,1)]
    colnames(term2gene) <- c("diseaseId", "geneId")
    return(list(go_term, go_mat, term2gene))
}

ensembl.to.entrez <- function(em_genes) {
    require(biomaRt)
    fname = 'em.to.ez.tsv'
    if (file.exists(fname)) {
        return(read.table(fname, header=T, sep="\t", stringsAsFactor=FALSE))
    }
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                            dataset = "dnovemcinctus_gene_ensembl",
                            host = "http://www.ensembl.org")

    genes <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
                values = em_genes, 
                mart = mart)
    write.table(genes, 'em.to.ez.tsv', sep="\t")
    return(genes)
}

ensembl.to.entrez.model.organism <- function(em_genes, org='mm') {
    m <- merge(data.frame(id=em_genes), genes.id2ids, by.x='id', by.y='ensemblID', all.x=TRUE)
    print(c(length(em_genes), dim(m)[1]))
    stopifnot(dim(m)[1] == length(em_genes))
    print('--merge ids with homolog data')
    print(head(m))
    if (org == 'mm') return(m[,c('id', 'EntrezGene.ID', 'Symbol')])
    else return(m[,c('id', 'EntrezGene.ID.1', 'Symbol.1')])
}

obtain.background.ez <- function(all_id_file, symbol=FALSE, species='mm') {
    features <- read.table(all_id_file, sep="\t", header=F, stringsAsFactor=FALSE)
    colnames(features) <- c("id", "symbol")
    features <- features[unlist(sapply(features[,1], function(x){grepl("^ENSDNOG", x)})),]
    genes.id <- ensembl.to.entrez.model.organism(features[,1], species)
    if (symbol) {
        genes.id <- genes.id[,c(3, 2)]
    } else {
        genes.id <- genes.id[,c(1, 2)]
    }
    print(head(genes.id))
    colnames(genes.id) <- c("id", "entrez")
    return(genes.id)
}

obtain.background.ez.id <- function(all_ids, symbol=FALSE, species='mm') {
    features <- merge(data.frame(gene_id=unique(all_ids)), gtf, by='gene_id')
    colnames(features) <- c("id", "symbol")
    features <- features[!duplicated(features), ]
    features <- features[unlist(sapply(features[,1], function(x){grepl("^ENSDNOG", x)})),]
    genes.id <- ensembl.to.entrez.model.organism(features[,1], species)
    print(head(genes.id))
    if (symbol) {
        genes.id <- genes.id[,c(3, 2)]
    } else {
        genes.id <- genes.id[,c(1, 2)]
    }
    print(head(genes.id))
    colnames(genes.id) <- c("id", "entrez")
    return(genes.id)
}


read.gene.list <- function(input_file) {
    df <- read.table(input_file, sep="\t", header=T, stringsAsFactor=FALSE)
    print(head(df))
    gene_list <- list()
    for (i in unique(df[,'cluster'])) {
        temp <- subset(df, df[,'cluster'] == i)
        temp <- temp[order(temp[,'p_val_adj']),]
        temp[,'gene'] <- unlist(sapply(temp[,'gene'], function(x){return(unlist(strsplit(x, '.', fixed=TRUE))[1])}))
        gene_list[[as.character(i)]] = data.frame(id=temp[,'gene'], order=unlist(sapply(temp[,'p_val_adj'], function(x){return(-log10(max(x, 0.0000000001)))})), logFC=temp[,'avg_logFC'])
    }
    print(gene_list)
    print(names(gene_list))
    print(length(gene_list))
    return(gene_list)
}

species = 'hs'

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    all_id_file = "../data/singlecell/final_pseudo_bulk_rna.rds"
    # all_id_file = "../data/bam/305264Solo.out/Gene/filtered/features.tsv"
    if (grepl('rds', all_id_file)) {
        x.sp <- readRDS(all_id_file)
        genes <- rownames(x.sp@assays[['RNA']])[!grepl("^gSpikein-ERCC", rownames(x.sp@assays[["RNA"]]))]
        genes.id <- obtain.background.ez.id(genes, FALSE, species)
    } else {
        genes.id <- obtain.background.ez(all_id_file, FALSE, species)
    }
    genes.id <- genes.id[!is.na(genes.id[,'entrez']),]
    print(head(genes.id))
} else {
    all_id_file = args[1]
}

# column: id and symbol used for the analysis

if (length(args) < 2) {
    if (nchar(input_file) == 0) {
        input_file = "../data/gene/rank_gene_atac.tsv"
        # input_file = "../data/gene/rank_gene_rna_mcocluster.tsv"
    }
} else {
    input_file = args[2]
}

if (length(args) < 3) {
    if (nchar(input_file) == 0) {
        header = 'rna_mcocluster_'
    }
} else {
    header = args[3]
}

# read gene2go matrix + annotation
go_info <- construct.go.data()
go_term <- go_info[[1]]
go_mat <- go_info[[2]]
term2gene <- go_info[[3]]

gene_list <- read.gene.list(input_file)

pvalue_thres = -log10(0.05)
# read all gene id and symbols from scRNA-seq data features.tsv

for (filt in c('', '_nonIEA')) {
    if (filt == "_nonIEA") {
        go_mat <- subset(go_mat, GO.arm[,4] != 'IEA')
        term2gene <- subset(term2gene, GO.arm[,4] != 'IEA')
    }
    for (i in names(gene_list)) {
        genes <- gene_list[[i]]
        m <- merge(genes, genes.id, by="id")
        m <- m[order(m[,'order'], decreasing=TRUE),]
        m <- m[!is.na(m$entrez),]
        stopifnot(length(unique(m[,'id'])) == dim(m)[1])
        max_rank = max(c(max(which(!is.na(m[,'order']))), max(m[,'order'] > pvalue_thres)))
        for (thres in c(10, 50, 100, 500, 1000)) {
            gene <- m[1:min(thres, max_rank),'entrez']
            x = enricher(gene, genes.id[,'entrez'], TERM2GENE=term2gene, TERM2NAME=go_term)
            write.table(x@result, file=paste0('enrich_cluster_', header, i, "_", thres, "_cp", filt, ".tsv"), sep="\t")
            x@result = x@result
            pdf(paste0('dotplot_', header, i, "_", min(max_rank, thres), "_cp", filt, ".pdf"), width=10, height=4)
            plot(dotplot(x, showCategory=10))
            dev.off()
            if (max_rank < thres) break
        }
        if (filt == '_nonIEA') next
        all_id_data <- m[,'logFC']
        names(all_id_data) <- m[,'entrez']
        temp <- m[m[,'order'] >= pvalue_thres, ]
        sig <- temp[,'logFC']
        names(sig) <- temp[,'entrez']
        for (map in pathways[,1]) {
            if (species == 'hs') {
                plot.KEGG.pathway(str_pad(map, 5, pad='0'), sig, paste0('cluster_', header, i, filt), 'hsa', all_id_data)
            } else {
                plot.KEGG.pathway(str_pad(map, 5, pad='0'), map, sig, paste0('cluster_', header, i, filt), 'mmu', all_id_data)
            }
        }
    }
}

