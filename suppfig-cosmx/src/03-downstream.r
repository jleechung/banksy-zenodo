rm(list=ls())

library(Banksy) 
library(data.table)
library(harmony)
library(Seurat)
library(SeuratWrappers) 
library(cowplot)
library(ggplot2)
library(repr)
sethw = function(h,w) options(repr.plot.height=h, repr.plot.width=w)

SEED = 1000

pal = c(pals::kelly()[-c(1,3)], pals::glasbey(n=32))

clip = function(x, low=0.02, high=0.98) {
    qt = quantile(x, c(low, high))
    x[x < qt[1]] = qt[1]
    x[x > qt[2]] = qt[2]
    x
}

ulen = function(x) length(unique(x))

getNeighborhood = function(seu, idents, cluster) {
    loc = seu@meta.data[,sdims]
    knn = kNN(loc, k = 15)
    tab = table(seu@meta.data[, idents][unique(as.vector(knn$id[which(seu@meta.data[, idents] == cluster),]))])
    tab = tab[!grepl(cluster, names(tab))]
    tab = tab[order(tab, decreasing = TRUE)]
    tab / sum(tab)
    sethw(4,12)
    plot(tab / sum(tab), ylab='Proportions', main=cluster)
}

make.unique_custom <- function(x, sep = ".", start = 1) {
  if (!is.character(x)) {
    stop("Input must be a character vector.")
  } 
  counts <- table(x)
  unique_x <- unique(x)
  new_x <- character(length(x))
  for (i in seq_along(unique_x)) {
    indices <- which(x == unique_x[i])
    if (length(indices) > 1) {
      new_x[indices] <- paste(unique_x[i], sep, start + seq_along(indices) - 1, sep = "")
    } else {
      new_x[indices] <- unique_x[i]
    }
  }
  new_x
}

avgCol = function(C1, C2) {
    library(zeallot)
    C1 = as.vector(col2rgb(C1))
    C2 = as.vector(col2rgb(C2))
    c(R1, G1, B1) %<-% C1
    c(R2, G2, B2) %<-% C2    
    c(R3, G3, B3) %<-% c(
        sqrt((R1^2+R2^2)/2),
        sqrt((G1^2+G2^2)/2),
        sqrt((B1^2+B2^2)/2)
        )
    rgb(R3, G3, B3, maxColorValue = 255)
}


seu = readRDS('../out/gsm3683.rds')

# 1.4, 1.1
nsp = as.numeric(seu@meta.data[,'RNA_snn_res.1.4'])
bank = as.numeric(seu@meta.data[,'BANKSY_snn_res.1.1'])
table(nsp)
table(bank)
bank = Banksy:::mapToSeed(val = bank, seed = nsp)
nsp = factor(nsp - 1)
bank = factor(bank - 1)
seu$nonspatial = nsp
seu$banksy = bank

#' collate meta.data df
mdf_nsp = cbind(Embeddings(seu, reduction = 'umap'), clusters = seu$nonspatial, seu@meta.data[, sdims])
colnames(mdf_nsp) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

mdf_bank = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$banksy, seu@meta.data[, sdims])
colnames(mdf_bank) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

mdf_fov = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$fov, seu@meta.data[, sdims])
colnames(mdf_fov) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

#' gen plots
pout = lapply(list(mdf_nsp, mdf_bank, mdf_fov), function(mdf) {
    
    #' spatial plot
    pt_size = 0.1
    pt_alpha = 0.5
    sp = ggplot(mdf, aes(x=dimx, y=dimy, col=clusters)) +
        geom_point(size=pt_size, alpha=pt_alpha) +
        theme_minimal() + 
        theme(legend.title=element_blank()) +
        coord_equal() + 
        scale_color_manual(values = pal) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 

    #' umap
    pt_size = 0.3
    pt_alpha = 0.5
    um = ggplot(mdf, aes(x=UMAP1, y=UMAP2, col=clusters)) +
        geom_point(size=pt_size, alpha=pt_alpha) +
        theme_minimal() + 
        theme(legend.title=element_blank()) +
        # coord_equal() + 
        scale_color_manual(values = pal) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 

    #' spatial wrapped
    pt_size = 0.01
    spw = ggplot(mdf, aes(x=dimx, y=dimy, col=clusters)) +
        geom_point(size=pt_size, alpha=pt_alpha) +
        theme_minimal() + 
        theme(legend.title=element_blank()) +
        coord_equal() + 
        scale_color_manual(values = pal) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        facet_wrap(~clusters, ncol = 4)
    
    list(um, sp, spw)
})

pout = unlist(pout, recursive=FALSE)
sout = pout[c(2,5,3,6)]
uout = pout[c(1,4,7)]

sethw(7,25)
uout[[1]] = uout[[1]] + labs(title='Non-spatial') + theme(plot.title = element_text(size=40), legend.position = 'none') 
uout[[2]] = uout[[2]] + labs(title='BANKSY') + theme(plot.title = element_text(size=40))
plot_grid(plotlist = uout, nrow = 1, rel_widths = c(0.9,1,1))

Idents(seu) = seu$banksy
seu = FindSubCluster(seu, cluster = 8, graph.name = 'BANKSY_snn', resolution = 0.1)
seu$sub.cluster = factor(seu$sub.cluster, levels = c(as.character(0:7), c('8_0', '8_1'), as.character(9:20)))
Idents(seu) = seu$sub.cluster

mdf_bank = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$sub.cluster, seu@meta.data[, sdims])
colnames(mdf_bank) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

spal = pal[c(1:9,25,10:21,22)]

#' umap
pt_size = 0.1
pt_alpha = 0.2
um = ggplot(mdf_bank, aes(x=UMAP1, y=UMAP2, col=clusters)) +
    geom_point(size=pt_size, alpha=pt_alpha) +
    theme_minimal() + 
    theme(legend.title=element_blank()) +
    scale_color_manual(values = spal) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 

# Fine and coarse predictions
labels_coarse = readRDS('../out/pred_coarse.rds')$labels
labels_fine = readRDS('../out/pred_fine.rds')$labels

tab = table(labels_coarse, seu$sub.cluster)
cluster_coarse = make.unique(rownames(tab)[max.col(t(tab))], sep = '.')

tab = table(labels_fine, seu$sub.cluster)
cluster_fine = make.unique_custom(rownames(tab)[max.col(t(tab))], sep = '.')
cluster_fine_refined = make.unique_custom(c(
    'Epithelial', 'Plasma', 'T helper', 'Epithelial', 'Epithelial', 'Plasma', 'Fibroblast', 'Stromal', 'Macrophage', 'Macrophage', 'Stromal', 
    'Smooth Muscle', 'Epithelial', 'Endothelial', 'B', 'Epithelial', 'Stromal', 'Mast', 'Pericyte', 'T reticular', 'Tuft',
    'Doublet'
))
cluster_fine_refined_short = make.unique_custom(c(
    'Epi', 'Plasma', 'T helper', 'Epi', 'Epi', 'Plasma', 'Fib', 'Stromal', 'Mac', 'Mac', 'Stromal', 
    'SMC', 'Epi', 'Endo', 'B', 'Epi', 'Stromal', 'Mast', 'Pericyte', 'T reticular', 'Tuft',
    'Doub'
))

seu$celltype_coarse = factor(seu$sub.cluster, labels = cluster_coarse)
pal_coarse = spal[order(levels(seu$celltype_coarse))]
seu$celltype_coarse = factor(seu$celltype_coarse, levels = sort(levels(seu$celltype_coarse)))

seu$celltype_fine_refined = factor(seu$sub.cluster, labels = cluster_fine_refined)
pal_fine = spal[order(levels(seu$celltype_fine_refined))]
seu$celltype_fine_refined = factor(seu$celltype_fine_refined, levels = sort(levels(seu$celltype_fine_refined)))

seu$celltype_fine_refined_short = factor(seu$sub.cluster, labels = cluster_fine_refined_short)
seu$celltype_fine_refined_short = factor(seu$celltype_fine_refined_short, levels = sort(levels(seu$celltype_fine_refined_short)))

mdf_bank = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$celltype_fine_refined, seu@meta.data[, sdims])
colnames(mdf_bank) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

#' umap
pt_size = 0.3
pt_alpha = 0.5
um_anno = ggplot(mdf_bank, aes(x=UMAP1, y=UMAP2, col=clusters)) +
    geom_point(size=pt_size, alpha=pt_alpha) +
    theme_minimal() + 
    theme(
        legend.title=element_blank(), 
        plot.title = element_text(size=40), 
        legend.text = element_text(size=20),
        legend.spacing.y = unit(0.11, 'cm')
    ) +
    scale_color_manual(values = pal_fine) +
    guides(colour = guide_legend(byrow=TRUE, ncol=1, override.aes = list(size=5, alpha = 1))) +
    labs(title = 'BANKSY') 

leg = ggpubr::get_legend(um_anno)

um_anno = um_anno + theme(legend.position='none') +
    geom_segment(
        x = -7, y = 4,
        xend = -5.7, yend = 3,
        lineend = "round", # See available arrow types in example above
        linejoin = "mitre",
        linewidth = 0.8, 
        arrow = arrow(type='closed', length = unit(0.1, "inches")),
        color = 'black'
    )

sethw(7,16)
plot_grid(uout[[1]], um_anno, leg, nrow=1, rel_widths = c(1,1,0.4))

sethw(7,16)
Idents(seu) = seu$celltype_fine_refined
VlnPlot(seu, features = 'nFeature_RNA', pt.size = 0, cols = pal_fine) + theme_minimal() + 
    theme(legend.position='none', 
          text = element_text(size=20),
          axis.text.x = element_text(size=18, angle=45, hjust=1)
         ) + 
    labs(title='Number of detected genes') + xlab('BANKSY clusters')

sethw(7,16)
Idents(seu) = seu$celltype_fine_refined
VlnPlot(seu, features = 'nCount_RNA', pt.size = 0, cols = pal_fine) + theme_minimal() + 
    theme(legend.position='none', 
          text = element_text(size=20),
          axis.text.x = element_text(size=18, angle=45, hjust=1)
         ) + 
    labs(title='Total transcripts') + xlab('BANKSY clusters')

Idents(seu) = seu$banksy
seu = FindSubCluster(seu, cluster = 8, graph.name = 'BANKSY_snn', resolution = 0.1)

seu$sub.cluster = factor(seu$sub.cluster, levels = c(as.character(0:7), '8_0', '8_1', as.character(9:20)))

#' collate meta.data df
mdf_nsp = cbind(Embeddings(seu, reduction = 'umap'), clusters = seu$sub.cluster, seu@meta.data[, sdims])
colnames(mdf_nsp) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

mdf_bank = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$sub.cluster, seu@meta.data[, sdims])
colnames(mdf_bank) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

mdf_fov = cbind(Embeddings(seu, reduction = 'umap_banksy'), clusters = seu$fov, seu@meta.data[, sdims])
colnames(mdf_fov) = c('UMAP1', 'UMAP2', 'clusters', 'dimx', 'dimy')

# BANKSY subplot
sub_bank = mdf_bank
sub_bank$clusters = factor(sub_bank$clusters, labels = c(rep('Other cells', 8), 'Mac.1', 'Mac.2', rep('Other cells', 12)))
sub_bank$clusters = factor(sub_bank$clusters, levels = c('Mac.1', 'Mac.2', 'Other cells'))
sub_bank = sub_bank[c(
    which(sub_bank$clusters == 'Other cells'),
    which(sub_bank$clusters == 'Mac.1'),
    which(sub_bank$clusters == 'Mac.2')
),]
sub_bank$size = 1
sub_bank$size[grep('Mac', sub_bank$clusters)] = 2

sp8 = ggplot(sub_bank, aes(x=dimx, y=dimy, col=clusters, size=size)) +
    geom_point(alpha=0.9) +
    theme_minimal() + 
    theme(legend.title=element_blank(), legend.position=c(0.9,0.9), legend.text = element_text(size=15)) +
    coord_equal() + 
    scale_color_manual(values = c(pal[c(9,25)], 'grey80')) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)), size='none') +
    theme(axis.text = element_blank()) + xlab(NULL) + ylab(NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_size(range=c(0.3,2))
sethw(6,15)
sp8

Idents(seu) = seu$sub.cluster
DefaultAssay(seu) = 'RNA'
markers_cluster = FindMarkers(seu, ident.1 = '8_0', ident.2 = '8_1', only.pos = TRUE)
getNeighborhood(seu, 'sub.cluster', '8_0')
markers1 = FindMarkers(seu, ident.1 = '1', only.pos = TRUE)
markers2 = FindMarkers(seu, ident.1 = '2', only.pos = TRUE)
markers0 = FindMarkers(seu, ident.1 = '0', only.pos = TRUE)
markers9 = FindMarkers(seu, ident.1 = '9', only.pos = TRUE)
markers13 = FindMarkers(seu, ident.1 = '13', only.pos = TRUE)
markers_cluster = setdiff(rownames(markers_cluster), c(rownames(markers1), rownames(markers2), rownames(markers0), rownames(markers9), rownames(markers13)))

seusub = seu[, which(seu$sub.cluster %in% c('8_0', '8_1'))]
sethw(3,10)
v_8_0_all = VlnPlot(seu, features = markers_cluster[1:2], cols = pal[c(1:9,25,10:21)], pt.size = 0, ncol=2)
sethw(3,6)
v_8_0 = VlnPlot(seusub, features = markers_cluster[1:2], cols = pal[c(9,25)], pt.size = 0, ncol=2)

Idents(seu) = seu$sub.cluster
DefaultAssay(seu) = 'RNA'
markers_cluster = FindMarkers(seu, ident.1 = '8_1', ident.2 = '8_0', only.pos = TRUE)
getNeighborhood(seu, 'sub.cluster', '8_1')
markers15 = FindMarkers(seu, ident.1 = '15', only.pos = TRUE)
markers10 = FindMarkers(seu, ident.1 = '10', only.pos = TRUE)
markers12 = FindMarkers(seu, ident.1 = '12', only.pos = TRUE)
markers17 = FindMarkers(seu, ident.1 = '17', only.pos = TRUE)
markers2 = FindMarkers(seu, ident.1 = '2', only.pos = TRUE)
markers_cluster = setdiff(rownames(markers_cluster), c(rownames(markers15), rownames(markers10), rownames(markers12), rownames(markers17), rownames(markers2)))

seusub = seu[, which(seu$sub.cluster %in% c('8_0', '8_1'))]
sethw(3,20)
v_8_1_all = VlnPlot(seu, features = markers_cluster[c(1,2,4,6)], cols = pal[c(1:9,25,10:21)], pt.size = 0, ncol=6)
sethw(3,12)
v_8_1 = VlnPlot(seusub, features = markers_cluster[c(1,2,4,6)], cols = pal[c(9,25)], pt.size = 0, ncol=6)

sethw(6,7)
seusub = seu[, which(seu$sub.cluster %in% c('8_0', '8_1'))]
Idents(seusub) = seusub$celltype_fine_refined_short
vlist = VlnPlot(seusub, features = c('ITGAX', 'TYK2', 'CD163', 'MRC1', 'GAS6', 'CD14'), cols = pal[c(9,25)], pt.size = 0, ncol=3, combine=FALSE)
vlist = lapply(vlist, function(x) x + theme(legend.position='none') + ylab(NULL) + xlab(NULL))
plot_grid(plotlist = vlist, ncol=3)

