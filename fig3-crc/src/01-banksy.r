
#' Load libs -------------------------------------------------------------------

library(Banksy)
library(gridExtra)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
library(future)
library(data.table)
library(plyr)
library(Palo)
library(ComplexHeatmap)

plan("multicore", workers = 12)
options(future.globals.maxSize = 6000 * 1024^2)

pal = c(pals::kelly()[-1], pals::glasbey())

pal[14] = pals::brewer.paired(3)[3]
pal[16] = pals::glasbey()[10]
pal[19] = pals::glasbey()[15]

#' Switches
RUN_PROCESS = FALSE
RUN_MARKERS = FALSE
RUN_UMAP = FALSE
RUN_SPATIAL = FALSE
RUN_PSEUDOBULK = FALSE

#' Read data -------------------------------------------------------------------

sample_prefix = 'HumanColonCancerPatient1'
gcm = readRDS(sprintf('../data/%s_gcm.rds', sample_prefix))
locs = readRDS(sprintf('../data/%s_locs.rds', sample_prefix))
mdata = readRDS(sprintf('../data/%s_metadata.rds', sample_prefix))
sample_prefix = paste0(sample_prefix, '_k30')

if (RUN_PROCESS) {

seu = CreateSeuratObject(counts = gcm, meta.data = cbind(locs, mdata), min.cells = 10)
seu = subset(seu, 
        nCount_RNA < quantile(nCount_RNA, 0.98) &  nCount_RNA > quantile(nCount_RNA, 0.05) 
)

#' Run clustering for 6 resolutions
RES = c(0.5, 0.6, 0.7, 0.8)

print(median(colSums(seu@assays$RNA@counts)))
seu = NormalizeData(seu, scale.factor = median(colSums(seu@assays$RNA@counts)))
seu = ScaleData(seu, features = rownames(seu))
seu = RunPCA(seu, features = rownames(seu), npcs = 50, verbose=FALSE)
seu = RunUMAP(seu, dims = 1:20, reduction.name = 'umap', verbose=FALSE, min.dist = 0.05, n.neighbors = 20)
seu = FindNeighbors(seu, dims = 1:20)
seu = FindClusters(seu, resolution = RES)

#' ## Run BANKSY with k_geom = 15, lambda = 0.2:
seu = RunBanksy(seu, M = 0, lambda = 0.2, k_geom = 15, assay = 'RNA', slot = 'data', 
		dimx = 'sdimx', dimy = 'sdimy', features = 'all', assay_name = 'M0')
seu = RunPCA(seu, features = rownames(seu), npcs = 50, verbose=FALSE)
seu = RunUMAP(seu, dims = 1:20, reduction.name = 'umap_m0', verbose=FALSE, min.dist = 0.05, n.neighbors = 20)
seu = FindNeighbors(seu, dims = 1:20)
seu = FindClusters(seu, resolution = RES)

#' ## Run BANKSY with k_geom = 15, 60 for M = 0, 1, lambda = 0.2
seu = RunBanksy(seu, M = 1, lambda = 0.2, k_geom = c(15,30), assay = 'RNA', slot = 'data',
                dimx = 'sdimx', dimy = 'sdimy', features = 'all', assay_name = 'M1')
seu = RunPCA(seu, features = rownames(seu), npcs = 50, verbose=FALSE)
seu = RunUMAP(seu, dims = 1:20, reduction.name = 'umap_m1', verbose=FALSE, min.dist = 0.05, n.neighbors = 20)
seu = FindNeighbors(seu, dims = 1:20)
seu = FindClusters(seu, resolution = RES)

#' Save:
saveRDS(seu, sprintf("../out/%s_seurat_object.rds", sample_prefix))
saveRDS(seu@reductions, sprintf("../out/%s_seurat_reductions.rds", sample_prefix))
saveRDS(seu@meta.data, sprintf("../out/%s_seurat_metadata.rds", sample_prefix))
saveRDS(seu@commands, sprintf("../out/%s_seurat_commands.rds", sample_prefix))
    
} else {

seu = readRDS(sprintf("../out/%s_seurat_object.rds", sample_prefix))

}

apply(seu@meta.data, 2, function(x) length(unique(x)))

sample_prefix = paste0(sample_prefix, '_res0.7')

#' Map -------------------------------------------------------------------------

rna = as.numeric(seu$RNA_snn_res.0.7)
rna = mapvalues(
    rna,
    from = as.numeric(names(table(rna))),
    to = as.numeric(names(sort(table(rna), decreasing = TRUE)))
)

m0 = as.numeric(seu$M0_snn_res.0.7)
m0 = mapvalues(
    m0,
    from = as.numeric(names(table(m0))),
    to = as.numeric(names(sort(table(m0), decreasing = TRUE)))
)

rna = Banksy:::mapToSeed(rna, m0)
m1 = as.numeric(seu$M1_snn_res.0.7)
m1 = Banksy:::mapToSeed(m1, m0)

N_CLUST = max(m1)

seu$RNA_snn_res.0.7 = factor(rna, levels=seq(min(rna), max(rna)))
seu$M0_snn_res.0.7 = factor(m0, levels=seq(min(m0), max(m0)))
seu$M1_snn_res.0.7 = factor(m1, levels=seq(min(m1), max(m1)))

set.seed(100)
pal = Palo(
	seu@meta.data[,c('sdimx', 'sdimy')],
	seu$M1_snn_res.0.7,
	pal[seq_len(N_CLUST)],
	rgb_weight = c(3,4,2)
)
order = order(as.numeric(names(pal)))
pal = pal[order]

print(pal)

#' Markers ---------------------------------------------------------------------

DefaultAssay(seu) = 'RNA'
n_markers = 10

Idents(seu) = seu$RNA_snn_res.0.7
if (RUN_MARKERS) {
markers_rna = FindAllMarkers(seu, only.pos = TRUE)
saveRDS(markers_rna, sprintf('../out/%s_markers_rna.rds', sample_prefix))
} else {
markers_rna = readRDS(sprintf('../out/%s_markers_rna.rds', sample_prefix))
}

markers_rna = data.table(markers_rna)
markers_rna = markers_rna[,head(.SD, n_markers),by=cluster]$gene

Idents(seu) = seu$M0_snn_res.0.7
if (RUN_MARKERS) {
markers_m0 = FindAllMarkers(seu, only.pos = TRUE)
saveRDS(markers_m0, sprintf('../out/%s_markers_m0.rds', sample_prefix))
} else {
markers_m0 = readRDS(sprintf('../out/%s_markers_m0.rds', sample_prefix))
}
markers_m0 = data.table(markers_m0)
markers_m0 = markers_m0[,head(.SD, n_markers),by=cluster]$gene

Idents(seu) = seu$M1_snn_res.0.7
if (RUN_MARKERS) {
markers_m1 = FindAllMarkers(seu, only.pos = TRUE)
saveRDS(markers_m1, sprintf('../out/%s_markers_m1.rds', sample_prefix))
} else {
markers_m1 = readRDS(sprintf('../out/%s_markers_m1.rds', sample_prefix))
}
markers_m1 = data.table(markers_m1)
markers_m1 = markers_m1[,head(.SD, n_markers),by=cluster]$gene


#' Run UMAP --------------------------------------------------------------------
if (RUN_UMAP) {

if (!file.exists('../plots')) dir.create('../plots', recursive = TRUE)

udata = data.frame(cbind(Embeddings(seu, reduction='umap'), Cluster=factor(seu$RNA_snn_res.0.7)))
print(head(udata))
colnames(udata) = c('UMAP_1', 'UMAP_2', 'Clusters')
udata$Clusters = factor(udata$Clusters, levels=seq(min(udata$Clusters), max(udata$Clusters)))
u1 = ggplot(udata, aes(x = UMAP_1, y = UMAP_2, col=Clusters)) +
        geom_point(size = 0.2, alpha = 0.3, shape='.') +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('Non-spatial clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

udata = data.frame(cbind(Embeddings(seu, reduction='umap_m0'), Cluster=factor(seu$M0_snn_res.0.7)))
colnames(udata) = c('UMAP_1', 'UMAP_2', 'Clusters')
udata$Clusters = factor(udata$Clusters, levels=seq(min(udata$Clusters), max(udata$Clusters)))
u2 = ggplot(udata, aes(x = UMAP_1, y = UMAP_2, col=Clusters)) +
        geom_point(size = 0.2, alpha = 0.3, shape='.') +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M0 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')


udata = data.frame(cbind(Embeddings(seu, reduction='umap_m1'), Cluster=factor(seu$M1_snn_res.0.7)))
colnames(udata) = c('UMAP_1', 'UMAP_2', 'Clusters')
udata$Clusters = factor(udata$Clusters, levels=seq(min(udata$Clusters), max(udata$Clusters)))
u3 = ggplot(udata, aes(x = UMAP_1, y = UMAP_2, col=Clusters)) +
       geom_point(size = 0.2, alpha = 0.3, shape='.') +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M1 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

pdf(sprintf('../plots/%s_umap.pdf', sample_prefix), height = 60, width = 185)
grid.arrange(u1, u2, u3, ncol=3)
dev.off()

}

#' Spatial plots ---------------------------------------------------------------

if (RUN_SPATIAL) {

plot.df <- seu@meta.data[, grep('snn_res.0.7|sdim', names(seu@meta.data))]
head(plot.df)

nid <- grep('RNA', names(plot.df))
kid <- which(table(plot.df[,nid])>10)
plot.df.nsp <- plot.df[plot.df[,nid] %in% kid,]
head(plot.df.nsp)
splot_nonspatial <- ggplot(plot.df.nsp, aes(x = sdimx, y = sdimy, col=as.factor(RNA_snn_res.0.7))) + 
	geom_point(size = 0.01, alpha = 0.3, shape='.') + 
	scale_color_manual(values = pal) +
	guides(color = guide_legend(override.aes = list(size = 3))) + 
	xlab('x coordinates') + ylab('y coordinates') + ggtitle('Non-spatial clustering') + 
	Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

splot_banksy <- ggplot(plot.df, aes(x = sdimx, y = sdimy, col=as.factor(M0_snn_res.0.7))) +
	geom_point(size = 0.01, alpha = 0.3, shape='.') +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M0 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')


splot_banksy1 <- ggplot(plot.df, aes(x = sdimx, y = sdimy, col=as.factor(M1_snn_res.0.7))) +
        geom_point(size = 0.01, alpha = 0.3, shape='.') +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M1 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')


splot_nonspatial_wrap <- splot_nonspatial + facet_wrap(~ RNA_snn_res.0.7)
splot_banksy_wrap <- splot_banksy + facet_wrap(~ M0_snn_res.0.7)
splot_banksy1_wrap <- splot_banksy1 + facet_wrap(~ M1_snn_res.0.7)

pdf(sprintf('../plots/%s_spatial.pdf', sample_prefix), height = 20, width = 65)
grid.arrange(splot_nonspatial, splot_banksy, splot_banksy1, ncol = 3)
dev.off()

pdf(sprintf('../plots/%s_spatial_wrap.pdf', sample_prefix), height = 60, width = 200)
grid.arrange(splot_nonspatial_wrap, splot_banksy_wrap, splot_banksy1_wrap, ncol = 3)
dev.off()

plot.df.nsp = plot.df.nsp[plot.df.nsp$sdimx > 2000 & plot.df.nsp$sdimx < 6000 & plot.df.nsp$sdimy > 3000 & plot.df.nsp$sdimy < 5500,]
plot.df = plot.df[plot.df$sdimx > 2000 & plot.df$sdimx < 6000 & plot.df$sdimy > 3000 & plot.df$sdimy < 5500,]

splot_nonspatial <- ggplot(plot.df.nsp, aes(x = sdimx, y = sdimy, col=as.factor(RNA_snn_res.0.7))) +
        geom_point(size = 0.2, alpha = 0.6) +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('Non-spatial clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

splot_banksy <- ggplot(plot.df, aes(x = sdimx, y = sdimy, col=as.factor(M0_snn_res.0.7))) +
        geom_point(size = 0.2, alpha = 0.6) +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M0 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')


splot_banksy1 <- ggplot(plot.df, aes(x = sdimx, y = sdimy, col=as.factor(M1_snn_res.0.7))) +
        geom_point(size = 0.2, alpha = 0.6) +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M1 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

splot_nonspatial_wrap <- splot_nonspatial + facet_wrap(~ RNA_snn_res.0.7)
splot_banksy_wrap <- splot_banksy + facet_wrap(~ M0_snn_res.0.7)
splot_banksy1_wrap <- splot_banksy1 + facet_wrap(~ M1_snn_res.0.7)

pdf(sprintf('../plots/%s_spatial_subset.pdf', sample_prefix), height = 10, width = 35)
grid.arrange(splot_nonspatial, splot_banksy, splot_banksy1, ncol = 3)
dev.off()

pdf(sprintf('../plots/%s_spatial_subset_wrap.pdf', sample_prefix), height = 20, width = 70)
grid.arrange(splot_nonspatial_wrap, splot_banksy_wrap, splot_banksy1_wrap, ncol = 3)
dev.off()

leg <- ggplot(plot.df, aes(x = sdimx, y = sdimy, col=as.factor(M1_snn_res.0.7))) +
        geom_point(size = 0.2, alpha = 1) +
        scale_color_manual(values = pal) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        xlab('x coordinates') + ylab('y coordinates') + ggtitle('M1 clustering') +
        Banksy:::theme_blank(legend.text.size = 10, main.size = 20, legend.pos = 'right')

leg <- cowplot::get_legend(leg)
pdf(sprintf('../plots/%s_legend.pdf', sample_prefix), height=8, width=5)
grid::grid.draw(leg)
dev.off()

}

#' Pseudobulk ------------------------------------------------------------------

if (RUN_PSEUDOBULK) {

DefaultAssay(seu) = 'RNA'
pbulk = AverageExpression(seu, group.by = 'orig.ident')
saveRDS(pbulk$RNA, '../out/HumanColonCancerPatient1_pseudobulk.rds')
    
}