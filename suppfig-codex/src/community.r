rm(list=ls())
# remotes::install_github('prabhakarlab/Banksy@6d76280')
# remotes::install_github('jleechung/seurat-wrappers@feat-aft')

library(Banksy) 
library(data.table)
library(harmony)
library(Seurat)
library(SeuratWrappers) 
library(future)
plan("multicore", workers = 20)
options(future.globals.maxSize = 6000 * 1024^2)

SEED = 1000

pal = c(pals::kelly()[-c(1,3)], pals::glasbey(n=32))

clip = function(x, low=0.02, high=0.98) {
    qt = quantile(x, c(low, high))
    x[x < qt[1]] = qt[1]
    x[x > qt[2]] = qt[2]
    x
}

ulen = function(x) length(unique(x))

path = '../data/hickey-intestine/CODEX_HuBMAP_alldata_Dryad.csv.gz'
df = data.frame(fread(path))
table(df$unique_region)
region = unique(df$unique_region)
region = region[grep('B006_Ascending', region)]
print(region)

for (i in seq(length(region))) {

    curr_region = region[i]
    message(curr_region)
    dfr = df[df$unique_region == curr_region,]
    gcm_sid1 = grep('MUC2', colnames(dfr))
    gcm_eid1 = grep('CD161', colnames(dfr))
    gcm_sid2 = grep('OLFM4', colnames(dfr))
    gcm_eid2 = grep('MUC6', colnames(dfr))
    gene_cols = c(seq(gcm_sid1, gcm_eid1), seq(gcm_sid2, gcm_eid2))
    gcm = t(as.matrix(dfr[, gene_cols]))
    loc = dfr[, setdiff(seq(ncol(dfr)), gene_cols)]
    colnames(gcm) = paste0('cell_', 1:ncol(gcm))
    rownames(loc) = colnames(gcm)

    print(dim(gcm))
    print(dim(loc))

    #' Init. seurat object 
    seu = CreateSeuratObject(counts = gcm, meta.data = loc)
    res_grid = seq(0.01, 0.5, by = 0.01)

    #' Non-spatial
    seu = SetAssayData(seu, assay.type = 'RNA', new.data = gcm, slot = 'scale.data')
    seu = RunPCA(seu, features = rownames(seu), npcs = 20, verbose = FALSE)
    seu = RunUMAP(seu, reduction = 'pca', reduction.name = 'umap', dims = 1:20, verbose = FALSE, 
                  umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
    message('Clustering')
    seu = FindNeighbors(seu, reduction = 'pca', dims = 1:20, k.param = 50, verbose = FALSE)
    seu = FindClusters(seu, resolution = res_grid, verbose = FALSE, random.seed = SEED)

    #' banksy
    seu = RunBanksy(seu, lambda = 0.8, M = 1, dimx = 'x', dimy = 'y', features = 'all', 
                    k_geom = c(15,30), assay = 'RNA', slot = 'data', verbose = FALSE)
    seu = RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 20, reduction.name = 'pca_banksy', verbose = FALSE)
    seu = RunUMAP(seu, reduction = 'pca_banksy', reduction.name = 'umap_banksy', dims = 1:20, verbose = FALSE, 
                  umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
    message('Clustering')
    seu = FindNeighbors(seu, reduction = 'pca_banksy', dims = 1:20, k.param = 50, verbose = FALSE) 
    seu = FindClusters(seu, resolution = res_grid, verbose = FALSE, random.seed = SEED)

    message('Saving to ', sprintf('../out/seu-%s-community.rds', curr_region))
    saveRDS(seu, sprintf('../out/seu-%s-community.rds', curr_region))

}
