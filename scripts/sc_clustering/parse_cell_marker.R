source("../scripts/sc_clustering/multi_sc_base.R")
require(stringr)

marker <- read.table("../data/gene/all_cell_markers_lu.txt", quote="\"", sep="\t", header=T, stringsAsFactor=F) #211210 cellmarker
MOUSE = FALSE
BLOOD = TRUE
if (MOUSE) {
    marker <- marker[marker[,1] == 'Mouse',]
    marker <- marker[marker[,4] == 'Normal',]
    id_index=4
    symbol_index=3
} else {
    marker <- marker[marker[,1] == 'Human',]
    marker <- marker[marker[,4] == 'Normal',]
    id_index=2
    symbol_index=1
}
if (BLOOD) {
    bm_marker <- marker[marker[,2] == 'Blood',]
} else {
    bm_marker <- marker[marker[,2] == 'Bone marrow',]
}

bm_marker = apply(bm_marker, c(1, 2), function(x) gsub('11-Sep', 'Sept11', x, fixed=TRUE))

load("../bulk_transcriptome/data/GO.arm.Rdata")
load("../bulk_transcriptome/datahomologs.Rdata")
gtf = gene.region.data()
gtf[,6] = unlist(sapply(gtf[,6], toupper))
print(dim(bm_marker))

parsed_marker = NULL
for (i in 1:dim(bm_marker)[1]) {
    entrez = bm_marker[i,]
    print(bm_marker[i,])
    print(length(bm_marker[i,]))
    gene_symbols = unlist(strsplit(bm_marker[i, 9], ', ', fixed=TRUE))
    gene_symbols = unlist(sapply(gene_symbols, toupper))
    gene_ids = unlist(strsplit(bm_marker[i,10], ', ', fixed=TRUE))
    print(gene_symbols)
    print(gene_ids)
    ids = unlist(sapply(1:length(gene_ids), function(x) {
        gene_symbol = gsub(']', '', gsub('[', '', gene_symbols[x], fixed=TRUE), fixed=TRUE)
    if (!is.na(gene_symbol) && any(gene_symbol == gtf[,6])) {
        return(gtf[which(gene_symbol == gtf[,6])[1],5])
    } else {
        id_string = gsub(']', '', gsub('[', '', gene_ids[x], fixed=TRUE), fixed=TRUE)
        ind = unlist(sapply(id_string, function(x){return(which(genes.id2ids[,id_index] == x))}))
        if (length(ind) == 0 || is.na(ind)) return("")
        # print(genes.id2ids[ind])
        return(genes.id2ids[ind,id_index])
    }}))
    ids = ids[!is.na(ids)]
    if (length(ids) == 0) next
    ids = as.character(ids)
    mouse_symbols = unlist(sapply(ids, function(x) {
        return(genes.id2ids[genes.id2ids[,5] == x,symbol_index])
    }))
    vec = c(bm_marker[i,1], bm_marker[i,2], bm_marker[i,6], length(ids)/length(gene_ids), length(ids), paste0(bm_marker[i,8]), paste0(bm_marker[i,9]), paste0(ids, collapse=','),
                paste0(mouse_symbols, collapse=','))
    
    print(vec)
    parsed_marker = rbind(parsed_marker, vec)
    rownames(parsed_marker)[dim(parsed_marker)[1]] = i
}
rownames(parsed_marker) = 1:dim(parsed_marker)[1]
colnames(parsed_marker) = c('species', 'tissue', 'celltype', 'reliability', 'detected', 'markers', 'symbol', 'gene_id', 'mouse.symbol')
write.table(parsed_marker, file='cell_marker_genes.tsv', sep="\t", quote=F)

