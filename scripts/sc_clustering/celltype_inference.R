require(plyr)
require(ggplot2)
require(Seurat)
source("../scripts/sc_clustering/multi_sc_base.R")

genes <<- read.gtf.file()

visualize.cell.type <- function(coembed, cluster='cluster') {
    require(scales)
    require(viridis)
    hex_codes1 <- rev(hue_pal()(length(unique(coembed@meta.data[,cluster]))))
    hex_codes2 <- rev(hue_pal()(length(unique(coembed@meta.data[,'rna_celltype']))))
    for (position in c('stack', 'fill')) {
        for (data in c('ATAC', 'RNA')) {
            for (filt in c('', '_filt')) {
                temp = coembed[,coembed@meta.data[,'orig.ident'] == data]
                if (filt != '') {
                    temp = temp[,temp@meta.data[,'rna_celltype'] != 'Undefined']
                }
                print(head(temp@meta.data))
                pair = temp@meta.data[,c('batch', cluster)]
                count_table = plyr::count(pair, vars=c("batch", cluster)) 
                count_table[,cluster] = factor(count_table[,cluster], levels=sort(unique(coembed@meta.data[,cluster])))
                pdf(paste0('celltype_population_', data, '_', position, filt, '.pdf'), width=10)
                g <- ggplot(count_table, aes_string(x='batch', y='freq', fill=cluster))+
                    geom_bar(stat='identity', position=position, color='black', width=0.9)+theme_classic()+scale_fill_manual(values=hex_codes1, drop=F)
                plot(g)
                dev.off()
                if (filt == '') {
                    pdf(paste0('celltype_population_', data, '_', position, filt, '_batch.pdf'), width=10)
                    g <- ggplot(count_table, aes_string(x='batch', y='freq', fill='batch'))+
                        geom_bar(stat='identity', width=0.9)+theme_classic()+scale_fill_viridis(discrete=TRUE)
                    plot(g)
                    dev.off()
                }
                print(head(count_table))
                temp = coembed[,coembed@meta.data[,'orig.ident'] == data]
                if (filt != '') {
                    temp = temp[,temp@meta.data[,'rna_celltype'] != 'Undefined']
                }
                for (labeling in c('', '_w')) {
                    for (pcluster in c('batch', cluster)) {
                        pair = temp@meta.data[,c(pcluster, 'rna_celltype')]
                        count_table = count(pair, vars=c(pcluster, "rna_celltype"))
                        count_table[,'rna_celltype'] = factor(count_table[,'rna_celltype'], levels=sort(unique(coembed@meta.data[,'rna_celltype'])))
                        pdf(paste0('celltype_population_', data, '_', position, filt, '_', pcluster, labeling, '_imputed.pdf'), width=10)
                        g <- ggplot(count_table, aes_string(x=pcluster, y='freq', fill='rna_celltype', label='rna_celltype'))+
                            geom_bar(stat='identity', position=position, color='black')+theme_classic()+scale_fill_manual(values=hex_codes2, drop=F)
                        if (labeling != '') {
                            g <- g+geom_text(position = position, vjust=+1, size=1)
                        }
                        plot(g)
                        dev.off()
                        if (position == 'fill') {
                            total = sapply(sort(unique(count_table[,1])), function(x) {
                                return(sum(count_table[count_table[,1] == x,3]))
                            })
                            count_table[,3] = sapply(1:dim(count_table)[1], function(x) {
                                return(count_table[x,3]/total[[count_table[x,1]]])
                            })
                        }
                        count_table = count_table[count_table[,'rna_celltype'] %in% c('Monocytes', 'Monocytes + T', 'CD14+ Monocytes'),]
                        print(count_table)
                        pdf(paste0('celltype_population_', data, '_', position, filt, '_', pcluster, labeling, '_imputed_mono.pdf'), width=10)
                        g <- ggplot(count_table, aes_string(x=pcluster, y='freq', fill='rna_celltype', label='rna_celltype'))+
                            geom_bar(stat='identity', position='stack', color='black')+theme_classic()+scale_fill_manual(values=hex_codes2, drop=F)
                        if (labeling != '') {
                            g <- g+geom_text(position = position, vjust=+1, size=1)
                        }
                        plot(g)
                        dev.off()
                    }
                }
            }
        }
    }
}

compute.celltype.pvalue <- function(coembed) {
    require(parallel)
    rep_times = 10000
    for (data in c('ATAC', 'RNA')) {
        for (label in c('rna_celltype', 'celltype')[1]) {
            temp = coembed[,coembed@meta.data[,'orig.ident']== data]
            temp = temp[,temp@meta.data[,label] != 'Undefined']
            for (batch in unique(coembed@meta.data[,'batch'])) {
                for (batch2 in unique(coembed@meta.data[,'batch'])) {
                    if (batch == batch2) next
                    fname = paste0('zscore_', data, '_', label, '_', batch, '_', batch2, '_', rep_times, '.tsv')
                    if (file.exists(fname)) next
                    size = dim(temp[,temp@meta.data[,'batch'] == batch2])[2]
                    print(c(data, label, batch, batch2, size))
                    simulated_count = do.call(rbind, mclapply(1:rep_times, function(x) {
                        sampled = sample(temp[,temp@meta.data[,'batch'] == batch]@meta.data[,label], size, replace=TRUE)
                        count_table = count(sampled)
                        return(count_table)
                    }, mc.cores=5))
                    print(simulated_count)
                    count_table = count(temp[,temp@meta.data[,'batch'] == batch2]@meta.data[,label])
                    print(count_table)
                    mat = NULL
                    for (celltype in count_table[,1]) {
                        bg = simulated_count[simulated_count[,1] == celltype,2]
                        value = count_table[count_table[,1] == celltype,2]
                        mat = rbind(mat, c(celltype, (value-mean(bg))/sd(bg), value, mean(bg), sd(bg)))
                    }
                    write.table(mat, file=fname, sep="\t", quote=F)
                }
            }
        }
    }
}

visualize.celltype.pvalue <- function(coembed) {
    require(ggplot2)
    require(scales)
    hex_codes1 <- rev(hue_pal()(length(unique(coembed@meta.data[,'celltype']))))
    hex_codes2 <- rev(hue_pal()(length(unique(coembed@meta.data[,'rna_celltype']))))
    rep_times = 10000

    for (data in c('ATAC', 'RNA')) {
        for (label in c('rna_celltype', 'celltype')[1]) {
            for (batch in 0:4) {
                for (batch2 in 0:4) {
                    fname = paste0('zscore_', data, '_', label, '_', batch, '_', batch2, '_', rep_times, '.tsv')
                    if (!file.exists(fname)) next
                    mat = read.table(fname, sep="\t", header=T)
                    colnames(mat) = c('celltype', 'zscore', 'cell', 'mean', 'sd')
                    mat[,'celltype'] = factor(mat[,'celltype'], levels=sort(unique(coembed@meta.data[,label])))
                    codes = hex_codes1
                    if (label == 'rna_celltype') code = hex_codes2
                    pdf(paste0('barplot_ct_', data, '_', label, '_', batch, '_', batch2, '.pdf'))
                    g <- ggplot(mat, aes(x=celltype, y=zscore, fill=celltype))+
                        geom_bar(stat='identity', position='stack')+theme_classic()+scale_fill_manual(values=codes, drop=F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                    plot(g)
                    dev.off()
                }
            }
        }
    }
}

construct.count.mat <- function(temp, label) {
    temp = temp[,temp@meta.data[,label] != 'Undefined']
    mat = NULL
    for (batch in unique(temp@meta.data[,'batch'])) {
        count_table = count(temp[,temp@meta.data[,'batch'] == batch]@meta.data[,label])
        size = sum(count_table[,2])
        count_table[,2] = count_table[,2]/size
        count_table = cbind(count_table, ind=rep(batch, dim(count_table)[1]))
        mat = rbind(mat, count_table)
    }
    colnames(mat) = c('celltype', 'ratio', 'ind')
    return(mat)
}

plot.zscore.dist <- function(mat, codes, header, count_mat) {
    require(ggsci)
    pdf(paste0(header, '.pdf'), useDingbats=FALSE)
    g <- ggplot(mat, aes(x=celltype, y=zscore, fill=celltype))+
        geom_bar(stat='identity', position='stack', color='black')+theme_classic()+scale_fill_manual(values=codes, drop=F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot(g)
    dev.off()
    mat = mat[order(mat[,'zscore']),]
    mat[,'monocyte'] = unlist(sapply(mat[,'celltype'], function(x) {
        if (any(x == c('CD14+ Monocytes', 'Monocytes', 'Monocytes + T'))) return('1 Monocyte')
        else if (any(x == 'DC')) return('5 DC')
        else if (x == 'NK') return('4 NK')
        else if (any(x == c('Plasma B', 'B'))) return('3 B cell')
        else if (x == 'HSC') return('6 HSC')
        else return('2 T cell')
    }))
    mat = cbind(mat, x=paste0(1:dim(mat)[1], ' ', mat[,'celltype']))
    mat[,'x'] = factor(mat[,'x'], levels=mat[,'x'])
    pdf(paste0(header, '_ord.pdf'), useDingbats=FALSE)
    g <- ggplot(mat, aes(x=x, y=zscore, fill=celltype))+
        geom_bar(stat='identity', position='stack', color='black')+theme_classic()+scale_fill_manual(values=codes, drop=F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        coord_flip()
    plot(g)
    dev.off()
    m <- merge(mat, count_mat, by='celltype')
    m[,'ratio'] = m[,'ratio']*100
    print(head(m))
    pdf(paste0(header, '_circ.pdf'), useDingbats=FALSE)
    g <- ggplot(m, aes(x=x, y=zscore, fill=celltype, size=ratio))+
        geom_point(shape=21)+theme_classic()+scale_fill_manual(values=codes, drop=F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        coord_flip()+scale_size_continuous(range = c(1, 13))
    plot(g)
    dev.off()
    pdf(paste0(header, '_monocyte_pop.pdf'), useDingbats=FALSE)
    g <- ggplot(m, aes(x=x, y=zscore, color=monocyte))+
        geom_point(size=8)+theme_classic()+scale_color_npg()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        coord_flip()
    plot(g)
    dev.off()
}

visualize.celltype.zscore <- function(coembed) {
    require(ggplot2)
    require(scales)
    hex_codes1 <- rev(hue_pal()(length(unique(coembed@meta.data[,'celltype']))))
    hex_codes2 <- rev(hue_pal()(length(unique(coembed@meta.data[,'rna_celltype']))))
    for (data in c('ATAC', 'RNA')) {
        for (label in c('rna_celltype', 'celltype')[1]) {
            for (filt in c('', '_filt')) {
                count_mat = construct.count.mat(coembed[,coembed@meta.data[,'orig.ident']== data], label)
                fname = paste0('zscore_', data, '_', label, filt, '_conservative.tsv')
                if (!file.exists(fname)) next
                all_mat = read.table(fname, sep="\t", header=T)
                print(all_mat)
                colnames(all_mat) = c('celltype', 'zscore', 'cell', 'mean', 'sd', 'batch')
                for (batch in 0:4) {
                    mat = all_mat[all_mat[,'batch'] == batch,]
                    if (dim(mat)[1] == 0) next
                    mat[,'celltype'] = factor(mat[,'celltype'], levels=sort(unique(coembed@meta.data[,label])))
                    codes = hex_codes1
                    if (label == 'rna_celltype') codes = hex_codes2
                    plot.zscore.dist(mat, codes, paste0('barplot_ct_cons_', data, '_', label, '_', batch, filt), count_mat[count_mat[,'ind'] == batch,])
                }
                all_mat[,'celltype'] = factor(all_mat[,'celltype'], levels=sort(unique(coembed@meta.data[,label])))
                codes = hex_codes1
                if (label == 'rna_celltype') codes = hex_codes2
                pdf(paste0('barplot_ct_cons_', data, '_', label, '_', filt, '.pdf'), width=10, height=10)
                g <- ggplot(all_mat, aes(x=batch, y=zscore, fill=celltype))+
                    geom_bar(stat='identity', position='stack', color='black')+theme_classic()+scale_fill_manual(values=codes, drop=F)+facet_wrap(vars(celltype), scales="free")+
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                plot(g)
                dev.off()
            }
        }
    }
}

compute.celltype.zscore.wosampling <- function(coembed) {
    for (data in c('ATAC', 'RNA')[1]) {
        for (label in c('rna_celltype', 'celltype')[1]) {
            for (filt in c('', '_filt')) {
                mat = construct.count.mat(coembed[,coembed@meta.data[,'orig.ident']== data], label)
                omat = NULL
                for (celltype in unique(mat[,1])) {
                    temp_mat = mat[mat[,1] == celltype,]
                    bg = temp_mat[,2]
                    if (filt != '') {
                        bg = temp_mat[temp_mat[,3] != batch,2]
                    }
                    for (batch in temp_mat[,3]) {
                        value = temp_mat[temp_mat[,3] == batch,2]
                        omat = rbind(omat, c(celltype, (value-mean(bg))/sd(bg), value, mean(bg), sd(bg), batch))
                    }
                }
                write.table(omat, file=paste0('zscore_', data, '_', label, filt, '_conservative.tsv'), sep="\t", quote=F)
            }
        }
    }

}

extract.one.omics <- function(data, coembed, label) {
    temp = coembed[,coembed@meta.data[,'orig.ident'] == data]
    if (data == 'RNA') {
        temp = temp[,temp@meta.data[,'rna_celltype'] != 'Undefined']
    } else {
        temp = temp[,temp@meta.data[,label] != 'Undefined']
    }
    return(temp)
}

compute.rna.atac.celltype.pvalue <- function(coembed) {
    require(parallel)
    rep_times = 10000
    omics_list = c('ATAC', 'RNA')
    for (i in 1:2) {
        data = omics_list[i]
        data2 = omics_list[2/i]
        for (label in c('rna_celltype', 'celltype')[1]) {
            temp = extract.one.omics(data, coembed, label)
            temp2 = extract.one.omics(data2, coembed, label)
            label1 = label
            label2 = label
            for (batch in unique(temp@meta.data[,'batch'])) {
                for (batch2 in unique(temp2@meta.data[,'batch'])) {
                    if (batch != batch2) next
                    fname = paste0('zscore_multi_', data, '_', label, '_', batch, '_', batch2, '_', rep_times, '.tsv')
                    size = dim(temp2[,temp2@meta.data[,'batch'] == batch2])[2]
                    print(c(data, label, batch, batch2, size))
                    simulated_count = do.call(rbind, mclapply(1:rep_times, function(x) {
                        sampled = sample(temp[,temp@meta.data[,'batch'] == batch]@meta.data[,label1], size, replace=TRUE)
                        count_table = count(sampled)
                        return(count_table)
                    }, mc.cores=5))
                    print(simulated_count)
                    count_table = count(temp2[,temp2@meta.data[,'batch'] == batch2]@meta.data[,label2])
                    print(count_table)
                    mat = NULL
                    for (celltype in count_table[,1]) {
                        bg = simulated_count[simulated_count[,1] == celltype,2]
                        value = count_table[count_table[,1] == celltype,2]
                        mat = rbind(mat, c(celltype, (value-mean(bg))/sd(bg), value, mean(bg), sd(bg)))
                    }
                    write.table(mat, file=fname, sep="\t", quote=F)
                }
            }
        }
    }
}


coembed <- readRDS('final_integrated.rds')

# gmat = t(coembed[,coembed@meta.data[,'orig.ident'] == 'ATAC']@assays[['RNA']]@data)
# cluster = coembed[,coembed@meta.data[,'orig.ident'] == 'ATAC']@meta.data[,'rna_celltype']
# plot.marker.mean(t(coembed[,coembed@meta.data[,'orig.ident'] == 'ATAC']@assays[['RNA']]@data), cluster, '_imputed')

# gmat = t(coembed[,coembed@meta.data[,'orig.ident'] == 'RNA']@assays[['RNA']]@data)
# cluster = coembed[,coembed@meta.data[,'orig.ident'] == 'RNA']@meta.data[,'rna_celltype']
# plot.marker.mean(t(coembed[,coembed@meta.data[,'orig.ident'] == 'RNA']@assays[['RNA']]@data), cluster, '_rna')

# cluster = coembed[,coembed@meta.data[,'orig.ident'] == 'RNA']@meta.data[,'cluster']
# plot.marker.mean(t(coembed[,coembed@meta.data[,'orig.ident'] == 'RNA']@assays[['RNA']]@data), cluster, '_rnac')

visualize.cell.type(coembed)
# compute.celltype.pvalue(coembed)

compute.celltype.zscore.wosampling(coembed)
visualize.celltype.zscore(coembed)
# visualize.celltype.pvalue(coembed)
# compute.rna.atac.celltype.pvalue(coembed)