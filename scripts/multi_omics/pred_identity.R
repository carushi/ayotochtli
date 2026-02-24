require(parallel)
require(jaccard)

DIR <<- "../bulk_transcriptome/data/"
ASE = FALSE
if (ASE) {
    load(paste0(DIR, "exprs_all.Rdata"))
    load(paste0(DIR, "metadata.Rdata"))
    load(paste0(DIR, "ase_ratios.test_train.Rdata"))
    load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
} else {
    load(paste0(DIR, "cpm_strand_comb.Rdata"))
    load(paste0(DIR, "counts_strand_comb.Rdata"))
    load(paste0(DIR, "gene_annotations_v0_95_mod.Rdata"))
    load(paste0(DIR, "armadillo_helper.Rdata"))
    load(paste0(DIR, "armadillo_hem.Rdata"))
    load(paste0(DIR, "ref.strand.test_train.Rdata"))
}
source("../scripts/helper.r")
source("../scripts/sc_clustering/multi_sc_base.R")
source("../scripts/annotation/gene_enrichment.R")
source("../scripts/multi_omics/bulk_pred_utils.R")
source("../scripts/multi_omics/sc_pred_utils.R")
load(paste0(DIR, "Xscaffolds.Rdata"))
scaffoldsX = scaffoldsX.prop3
scaffoldsX = scaffoldsX.sub


quad_names = get.quad.names(X.cpm.all)
fname = 'data_exp_rank.rds'
if (file.exists(fname)) {
    df = readRDS(fname)
} else {
    df = compute.exp.rank(X.cpm.all)
    saveRDS(df, 'data_exp_rank.rds')
    saveRDS(X.cpm.all, 'data_exp.rds')
}

# if (TRUE) { #BULK rank
#     rank.test = do.call(cbind, lapply(1:3, function(x) {
#         mat = df[['test']][[paste0('X.tt', x)]]
#         colnames(mat) = as.vector(sapply(colnames(mat), function(sequence) {
#             prefix = substr(sequence, 1, 4)
#             ident = substr(sequence, 6, 6)
#             return(paste0(prefix, '_', ident, '_', x))
#         }))
#         return(mat)
#     }))
#     rank.train = do.call(cbind, lapply(1:3, function(x) {
#         mat = df[['training']][[paste0('X.tr', x)]]
#         colnames(mat) = as.vector(sapply(colnames(mat), function(sequence) {
#             prefix = substr(sequence, 1, 4)
#             ident = substr(sequence, 6, 6)
#             return(paste0(prefix, '_', ident, '_', x))
#         }))
#         return(mat)
#     }))
#     saveRDS(rank.test, 'bulk_exp_rank_test.rds')
#     saveRDS(rank.train, 'bulk_exp_rank_train.rds')
#     # plot.heatmap.binary(filt.mat(rank.test), 'rank', col_fun= colorspace::diverging_hsv(10))
# }

# if (TRUE) { #BULK correlation
#     rownames(cor.all)     = rownames(X.cpm.all)
#     rownames(predcor.all) = rownames(X.cpm.all)
#     rownames(pred.all)    = rownames(X.cpm.all)
#     # How many idents are accurately identified
#     pred.all = prediction.all.comb(df, quad_names)
#     pred.each = prediction.one.comb(df, quad_names)
#     write.table(pred.each, file='gene_predictability.tsv', sep="\t")
#     plot.heatmap.binary(pred.each, 'predictability', col_fun=c('darkgray', 'lightgray', 'red'))
#     analyze.mean.var.pred(pred.each)
#     saveRDS(pred.all, 'prediction_all_comb.rds')
#     saveRDS(pred.each, 'prediction_each_comb.rds')
#     print(head(pred.all[1:10, 1:10]))
#     predcor.all = rank.correlation.trio(df, quad_names)
#     cor.all = rank.correlation.pair(X.cpm.all, quad_names)
#     colnames(predcor.all) = labels
#     colnames(pred.all)    = labels

# }

# if (FALSE) {
#     mat = prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt1']], c())
#     saveRDS(mat, 'bulk_pred_mat1.rds')
#     mat = prediction.ident.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt2']], c())
#     saveRDS(mat, 'bulk_pred_mat2.rds')
#     mat = prediction.ident.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt3']], c())
#     saveRDS(mat, 'bulk_pred_mat3.rds')
#     exit()
# }

sc_fname = 'data_exp_sc_rank.rds'
if (file.exists(sc_fname)) {
    df = readRDS(sc_fname)
} else {
    df = add.sc.tr(X.cpm.all, df, c(5), '../data/singlecell/final_pseudo_bulk_rna.rds')
    saveRDS(df, sc_fname)
}


if (TRUE) {
    mat = prediction.sc.each.comb(df, c(5), 'X.tt4', c(5), -4)
    saveRDS(mat, 'prediction_each_sc.rds')
    plot.heatmap.binary(mat, 'predictability_sc', col_fun=c('darkgray', 'lightgray', 'red'))
    # mat = prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt6']], c(), c(5), -4)
    # mat = prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt7']], c(), c(5), -4)
    # mat = prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt8']], c(), c(5), -4)
}

compute.gene.list.enrichment <- function(pair1, genes, header) {
    require(viridis)
    require(jaccard)
    seta = which(!is.na(pair1) & pair1 > 0)
    background = genes[!is.na(pair1)]
    gene.set.enrichment(data.frame(gene_id=genes[seta]), data.frame(gene_id=background), paste0(header))
    print('enriched')
}

compute.gene.list.difference <- function(pair1, pair2, gene_list, header) {
    require(viridis)
    require(jaccard)
    seta = which(pair1 > 0)
    setb = which(pair2 > 0)
    common = intersect(seta, setb)
    spe_a = setdiff(seta, setb)
    spe_b = setdiff(setb, seta)
    print(c(length(common), length(spe_a), length(spe_b)))
    background = gene_list[!is.na(pair1) & !is.na(pair2)]
    print(length(background))
    if (length(common) > 0)
        gene.set.enrichment(data.frame(gene_id=gene_list[common]), data.frame(gene_id=background), paste0(header, '_common'))
    if (length(spe_a) > 0)
        gene.set.enrichment(data.frame(gene_id=gene_list[spe_a]), data.frame(gene_id=background), paste0(header, '_set1'))
    if (length(spe_b) > 0) 
        gene.set.enrichment(data.frame(gene_id=gene_list[spe_b]), data.frame(gene_id=background), paste0(header, '_set2'))
}

write.module.data <- function(mat) {
    for (i in 1:dim(mat)[2]) {
        fname = paste0('module_pred_', colnames(mat)[i], '.txt')
        genes = rownames(mat)[mat[,i] > 0]
        genes = genes[!is.na(genes)]
        cat(paste0("# ", colnames(mat)[i], "\n"), file=fname, append=FALSE)
        cat(paste0(genes, collapse="\n"), file=fname, append=TRUE)
    }    
}

if (TRUE) {
    pred.each = readRDS('prediction_each_comb.rds')
    pred.each = cbind(pred.each, readRDS('prediction_each_sc.rds'))
    step = 4
    print(pred.each)
    for (pair1 in seq(1, dim(pred.each)[2])) {
        print(pair1)
        if (!any(colnames(pred.each)[pair1] == c('16-9_2_1', '16-9_2_2', '16-9_2_3'))) next
        print(colnames(pred.each)[pair1])
        compute.gene.list.enrichment(pred.each[,pair1], rownames(pred.each), paste0(colnames(pred.each)[pair1]))
    }

    for (pair1 in seq(1, dim(pred.each)[2]-1, 1)) {
        for (pair2 in seq(pair1+1, dim(pred.each)[2], 1)) {
            print(pair1)
            print(pair2)
            print(c(colnames(pred.each)[pair1], colnames(pred.each)[pair2]))
            compute.gene.list.difference(pred.each[,pair1], pred.each[,pair2], rownames(pred.each), paste0(colnames(pred.each)[pair1], '_', paste0(colnames(pred.each)[pair2])))
        }
    }
}


if (TRUE) {
    pred.each = readRDS('prediction_each_comb.rds')
    pred.each = cbind(pred.each, readRDS('prediction_each_sc.rds'))
    dname = './'
    files = list.files(path=dname, pattern='^rank_*')
    for (file in files) {
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
        print(head(df))
        print(prefix)
    }

}

exit()
if (FALSE) {
    for (i in 1:9) {
        plot.ranking.ratio(df[['test']][[paste0('X.tt', i)]], i)
    }
    exit()
}
if (FALSE) {
    # random partial gene set 
    # prediction.ident.different.gene.set(df[['training']][['X.tr1']], df[['test']][['X.tt4']], list(), c(5), -4, step=100)
    # feature set
    source("../scripts/multi_omics/bulk_pred_utils.R")
    prediction.ident.different.gene.set(df[['training']][['X.tr1']], df[['test']][['X.tt4']], list(), c(5), -4, step=100)
    prediction.ident.different.gene.set(df[['training']][['X.tr1']], df[['test']][['X.tt1']], list(), c(1:5), 0, step=100, file='sampled.pred.score.tt1')
    prediction.ident.different.gene.set(df[['training']][['X.tr2']], df[['test']][['X.tt2']], list(), c(1:5), 0, step=100, file='sampled.pred.score.tt2')
    prediction.ident.different.gene.set(df[['training']][['X.tr3']], df[['test']][['X.tt3']], list(), c(1:5), 0, step=100, file='sampled.pred.score.tt3')
    plot.predictability.sampled.all("sampled.pred.score.tt1_*", header='sc_ave_sampled_tt1')
    plot.predictability.sampled.all("sampled.pred.score.tt2_*", header='sc_ave_sampled_tt2')
    plot.predictability.sampled.all("sampled.pred.score.tt3_*", header='sc_ave_sampled_tt3')

}


plot.predictability.sampled("sampled.pred.score.csv", header='sc_ave_sampled')
plot.predictability.sampled.all("sampled.pred.score.csv", header='sc_ave_sampled')

# plot.predictability.different.gene.set(bob=FALSE)
exit()

if (TRUE) {
    # predictability of each gene
    results = list()
    for (train in 1:3) {
        temp = list()
        for (name in names(df[['test']])) {
            if (any(name == c('X.tt1', 'X.tt2', 'X.tt3', 'X.ttNA')))
                next
            temp[[name]] = prediction.score.one.comb(df[['training']][[paste0('X.tr', train)]], df[['test']][[name]], c(), c(5), -4)
            plot.cum.predictability(temp[[name]], paste0('pred_score_', name, '_', train, '.pdf'))
            temp[[name]][is.na(temp)] = 0
            mat = prediction.ident.one.comb(df[['training']][[paste0('X.tr', train)]], df[['test']][[name]], c(), c(5), -4)
            saveRDS(mat, paste0('sc_pred_mat', train, '_', name, '.rds'))
        }
        plot.group.cum.predictability(temp, paste0('pred_score_adjusted_', train, '.pdf'), c('X.tt5', 'X.tt6', 'X.tt7', 'X.tt8'))
        plot.group.cum.predictability(temp, paste0('pred_score_celltype_', train, '.pdf'), c(paste0('X.tt', 1:9), 'X.tt_NA'), TRUE)
        # plot.group.cum.predictability(temp, paste0('pred_score_celltype_', train, '.pdf'), c('X.tt4'))
        # plot.group.cum.predictability(temp, paste0('pred_score_celltype_', train, '.pdf'), c('X.tt9'))
        for (name in names(temp)) {
            if (train == 1) {
                results[[name]] = temp[[name]]
            } else {
                results[[name]] = results[[name]] + temp[[name]]
            }
        }
    }
    for (name in names(results)) {
        results[[name]] = results[[name]]/3
    }
    # plot.group.cum.predictability(results, paste0('pred_score_adjusted_', 'ave', '.pdf'), c('X.tt5', 'X.tt6', 'X.tt7', 'X.tt8'))
    # plot.group.cum.predictability(results, paste0('pred_score_celltype_', 'ave', '.pdf'), c(paste0('X.tt', 1:9), 'X.tt_NA'), TRUE)
    plot.group.cum.predictability(results, paste0('pred_score_adjusted_', 'tt4', '.pdf'), c('X.tt4'))
    plot.group.cum.predictability(results, paste0('pred_score_adjusted_', 'tt9', '.pdf'), c('X.tt9'))
    saveRDS(results, paste0('sc_pred_score.rds'))
}
exit()


for (i in 1:3) {
    plot.pred.matrix(paste0('bulk_pred_mat', i, '.rds'), paste0('bulk', '_', i))
}
for (i in 1:3) {
    for (name in names(df[['test']])) {
        fname = paste0('sc_pred_mat', i, '_', name, '.rds')
        if (file.exists(fname)) {
            plot.pred.matrix(fname, paste0('sc_', name, '_', i))
        }
    }
}
exit()
# prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt4']], c(), c(5), -4)
# prediction.ident.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt4']], c(), c(5), -4)
# prediction.ident.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt4']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt5']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt5']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt5']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt6']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt6']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt6']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt7']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt7']], c(), c(5), -4)
prediction.ident.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt7']], c(), c(5), -4)

# score = prediction.score.one.comb(df[['training']][['X.tr1s']], df[['test']][['X.tt4']], c(), c(1), 0)
score = prediction.score.one.comb(df[['training']][['X.tr2s']], df[['test']][['X.tt1']], c(), c(1), +4)
plot.cum.predictability(score, 'pred_score_t1_2s.pdf')
score = prediction.score.one.comb(df[['training']][['X.tr1']], df[['test']][['X.tt4']], c(), c(1), 0)
score = score+prediction.score.one.comb(df[['training']][['X.tr2']], df[['test']][['X.tt4']], c(), c(1), 0)
score = score+prediction.score.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt4']], c(), c(1), 0)
score = score/3
plot.cum.predictability(score, 'pred_score_s_ave.pdf')
plot.hist.predictability(score, 'hist_score_s_ave.pdf')
prediction.score.one.comb(df[['training']][['X.tr3']], df[['test']][['X.tt1']], c(), c(5), 0)
