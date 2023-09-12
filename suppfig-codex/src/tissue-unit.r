
rm(list=ls())
# remotes::install_github('prabhakarlab/Banksy@6d76280')
# remotes::install_github('jleechung/seurat-wrappers@5f19815')

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
donor = 'B012'

for (i in seq(length(donor))) {

    curr_region = donor[i]
    message(curr_region)
    dfr = df[df$donor == curr_region,]
    fovs = c('B012_Ileum', 'B012_Right', 'B012_Trans')
    dfr = dfr[dfr$unique_region %in% fovs, ]
    print(dfr$unique_region)
    gcm_sid1 = grep('MUC2', colnames(dfr))
    gcm_eid1 = grep('CD161', colnames(dfr))
    gene_cols = c(seq(gcm_sid1, gcm_eid1))
    gcm = t(as.matrix(dfr[, gene_cols]))
    loc = dfr[, setdiff(seq(ncol(dfr)), gene_cols)]
    colnames(gcm) = paste0('cell_', 1:ncol(gcm))
    rownames(loc) = colnames(gcm)
    
    tmp = loc[,c('x','y','unique_region')]
    loc_dt = data.table(tmp)
    colnames(loc_dt) = c('x', 'y', 'group')
    loc_dt$group = as.numeric(factor(loc_dt$group))
    loc_dt[, x := x - min(x), by = group]
    global_max <- max(loc_dt$x) * 1.2
    loc_dt[, x := x + (group - 1) * global_max]
    loc$x = loc_dt$x
    loc$y = loc_dt$y

    print(dim(gcm))
    print(dim(loc))

    #' Init. seurat object 
    seu = CreateSeuratObject(counts = gcm, meta.data = loc)
    seu$group = factor(seu$unique_region)
    res_grid = seq(0.02, 0.06, by = 0.002)

    #' Non-spatial
    seu = SetAssayData(seu, assay.type = 'RNA', new.data = gcm, slot = 'scale.data')
    seu = RunPCA(seu, features = rownames(seu), npcs = 20, verbose = FALSE)
    seu = RunHarmony(seu, group.by.vars = 'group',  do_pca=FALSE, reduction = 'pca', 
                     max.iter.harmony = 20, reduction.save = 'harmony')
    seu = RunUMAP(seu, reduction = 'harmony', reduction.name = 'umap', dims = 1:20, verbose = FALSE, 
                  umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
    message('Clustering')
    seu = FindNeighbors(seu, reduction = 'harmony', dims = 1:20, k.param = 50, verbose = FALSE)
    seu = FindClusters(seu, resolution = res_grid, verbose = FALSE, random.seed = SEED)

    #' banksy
    seu = RunBanksy(seu, lambda = 0.8, M = 1, dimx = 'x', dimy = 'y', features = 'all', 
                    k_geom = c(15, 30), assay = 'RNA', slot = 'data', verbose = FALSE)
    seu = RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 20, reduction.name = 'pca_banksy', verbose = FALSE)
    seu = RunHarmony(seu, group.by.vars = 'group',  do_pca=FALSE, reduction = 'pca_banksy', 
                     max.iter.harmony = 20, reduction.save = 'harmony_banksy')
    seu = RunUMAP(seu, reduction = 'harmony_banksy', reduction.name = 'umap_banksy', dims = 1:20, verbose = FALSE, 
                  umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
    message('Clustering')
    seu = FindNeighbors(seu, reduction = 'harmony_banksy', dims = 1:20, k.param = 50, verbose = FALSE) 
    seu = FindClusters(seu, resolution = res_grid, verbose = FALSE, random.seed = SEED)

    message('Saving to ', sprintf('../out/seu-%s-tissue-unit.rds', curr_region))
    saveRDS(seu, sprintf('../out/seu-%s-tissue-unit.rds', curr_region))

}
