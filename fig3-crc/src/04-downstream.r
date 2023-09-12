
#' Load libs -------------------------------------------------------------------

library(Banksy)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(arrow)
library(dplyr)
library(SeuratWrappers)
library(Nebulosa)

# Plotting
library(gridExtra)
library(ggplot2)
library(ggalluvial)
library(gridExtra)
library(ggpubr)
library(repr)
library(ComplexHeatmap)
sethw = function(h,w) options(repr.plot.height=h, repr.plot.width=w)
sethw(10,12)

# Data handling
library(data.table)
library(mclust)
library(plyr)

# Colors
library(Palo)
library(randomcoloR)
ARI = adjustedRandIndex 

# Parallel
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 6000 * 1024^2)

pal = pals::kelly()[-1]
pal[14] = pals::brewer.paired(3)[3]
pal[16] = pals::glasbey()[10]
pal[19] = pals::glasbey()[15]
pal = c(pal, pals::glasbey()[1])

seu = readRDS('../out/HumanColonCancerPatient1_k30_seurat_object.rds')

rna = as.numeric(seu@meta.data$RNA_snn_res.0.7)
rna = mapvalues(
    rna,
    from = as.numeric(names(table(rna))),
    to = as.numeric(names(sort(table(rna), decreasing = TRUE)))
)

m0 = as.numeric(seu@meta.data$M0_snn_res.0.7)
m0 = mapvalues(
    m0,
    from = as.numeric(names(table(m0))),
    to = as.numeric(names(sort(table(m0), decreasing = TRUE)))
)

rna = Banksy:::mapToSeed(rna, m0)
m1 = as.numeric(seu@meta.data$M1_snn_res.0.7)
m1 = Banksy:::mapToSeed(m1, m0)
rna = Banksy:::mapToSeed(rna, m1)

seu$clust_nsp = factor(rna, levels=seq(min(rna), max(rna)))
seu$clust_m0 = factor(m0, levels=seq(min(m0), max(m0)))
seu$clust_m1 = factor(m1, levels=seq(min(m1), max(m1)))

table(seu$clust_nsp)
table(seu$clust_m0)
table(seu$clust_m1)

set.seed(100)
N_CLUST = 22
pal = Palo(
    seu@meta.data[,c('sdimx', 'sdimy')],
    seu$clust_m1,
    pal[seq_len(N_CLUST)],
    rgb_weight = c(3,4,2)  
)
order = order(as.numeric(names(pal)))
pal = pal[order]

scales::show_col(pal)
pal[c(9,20)] = pal[c(20,9)]

#' Load annotations ------------------------------------------------------------

singler_results = readRDS('../out/singleR_results_wilcox_10.rds')

tmp = as.character(singler_results$labels)

seu$singler = tmp

tab = table(seu$singler, seu$clust_m1)
singler_ct = make.unique(rownames(tab)[max.col(t(tab))], sep = '.')
seu$singler_ct = mapvalues(
    seu$clust_m1,
    from = 1:22,
    to = singler_ct
)

pal_named = pal
names(pal_named) = singler_ct

tmp = as.character(singler_results$pruned.labels)
tmp[is.na(tmp)] = 'Ambiguous'
tab = table(tmp, seu$clust_m1)

umap = Embeddings(seu, reduction = 'umap')
umap = data.frame(umap)
umap_data = cbind(umap, clust_nsp=seu$clust_nsp)
colnames(umap_data)[1:2] = c('UMAP_1', 'UMAP_2')
unsp = ggplot(umap_data, aes(x=UMAP_1, y=UMAP_2, col=clust_nsp)) + 
    geom_point(size = 0.05, alpha = 0.1, shape= '.') +
    scale_color_manual(values = pal) +
    theme_minimal() + guides(col = guide_legend(override.aes = list(size = 10))) +
    theme(text = element_text(size = 12), legend.title = element_blank(), legend.position = 'none',
          panel.grid = element_blank(), axis.line = element_line(color = 'grey60'), 
          axis.ticks = element_blank(), axis.text = element_blank(), plot.subtitle = element_text(size = 15), plot.title=element_text(size=20)) + 
        labs(title = '') +
        xlab('UMAP 1') + ylab('UMAP 2')

umap = Embeddings(seu, reduction = 'umap_m0')
umap = data.frame(umap)
umap_data = cbind(umap, clust_m0=seu$clust_m0)
colnames(umap_data)[1:2] = c('UMAP_1', 'UMAP_2')
um0 = ggplot(umap_data, aes(x=UMAP_1, y=UMAP_2, col=clust_m0)) + 
    geom_point(size = 0.05, alpha = 0.1, shape = '.') +
    scale_color_manual(values = pal) +
    theme_minimal() + guides(col = guide_legend(override.aes = list(size = 10))) +
    theme(text = element_text(size = 20)) + theme(legend.position='none')

umap = Embeddings(seu, reduction = 'umap_m1')
umap = data.frame(umap)
umap_data = cbind(umap, clust_m1=seu$clust_m1)
colnames(umap_data)[1:2] = c('UMAP_1', 'UMAP_2')
um1 = ggplot(umap_data, aes(x=UMAP_1, y=UMAP_2, col=clust_m1)) + 
    geom_point(size = 0.05, alpha = 0.1, shape = '.') +
    scale_color_manual(values = pal) +
    theme_minimal() + guides(col = guide_legend(override.aes = list(size = 10))) +
    theme(text = element_text(size = 12), legend.title = element_blank(), legend.position = 'none',
          panel.grid = element_blank(), axis.line = element_line(color = 'grey60'), 
          axis.ticks = element_blank(), axis.text = element_blank(), plot.subtitle = element_text(size = 15), plot.title=element_text(size=20)) + 
        labs(title = '') +
        xlab('UMAP 1') + ylab('UMAP 2')

sethw(12,25)
grid.arrange(unsp, um1, ncol = 2)

#' Labels ----------------------------------------------------------------------

labels = c(
    'cE01 (TA-like prolif)', 'cE02 (Stem-like)', 'cS01 (MMP3+ CAF)', 
    'cE03 (Stem-like)', 'cE04 (TA-like prolif)', 'cE05 (TA-like prolif)',
    'cE06 (Stem/TA-like)', 'cE07 (Enterocytes)', 'cE08 (Cycling)',
    'cI01 (Monocytes)', 'cE10 (Immature Goblet)', 'cS02 (GREM1+ CAF)',
    'cE11 (Enterocytes)', 'cS03 (Pericytes)', 'cI02 (Granulocytes)',
    'cS04 (Endothelial)', 'cI03 (Macrophage-like)', 'cI04 (CD4+ Treg)',
    'cS05 (Myofibro)', 'cE12 (Cycling)', 'cS06 (CAF stem niche)',
    'cE13 (Stem/TA-like)'
)

seu$celltypes = factor(seu$singler_ct, labels = labels)

labels_short = gsub(' \\(.*', '', labels)
seu$celltypes_short = factor(seu$singler_ct, labels = labels_short)

DefaultAssay(seu) = 'RNA'
Idents(seu) = seu$celltypes_short
celltype_markers = FindAllMarkers(seu)

marker_lst = split(x=celltype_markers, f=celltype_markers$cluster)
marker_lst = lapply(marker_lst, function(x) x[x$avg_log2FC>0,])
openxlsx::write.xlsx(
    marker_lst, 
    '../out/HumanColonCancerPatient1_markers.xlsx'
)
