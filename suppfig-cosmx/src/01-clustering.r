rm(list=ls())

library(Banksy) 
library(data.table)
library(harmony)
library(Seurat)
library(SeuratWrappers) 

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


#' set paths
gcm_path = '../data/GSM7473683_HC_b_exprMat_file.csv.gz'
meta_path = '../data/GSM7473683_HC_b_metadata_file.csv.gz'

cmeta = fread(meta_path)
cmeta$ID = paste0(cmeta$fov, '_', cmeta$cell_ID)

#' Gene cell matrix
gcm = fread(gcm_path)
gcm$ID = paste0(gcm$fov, '_', gcm$cell_ID)
gcm = gcm[ID %in% cmeta$ID]
cnm = gcm$ID
gcm = t(as.matrix(gcm[,3:962]))
colnames(gcm) = cnm
cmeta = cmeta[match(cmeta$ID, colnames(gcm))]
cmeta$nCount = clip(colSums(gcm))
cmeta$NODG = clip(colSums(gcm > 0))
cmeta$CPG = cmeta$nCount / cmeta$NODG
cmetadf = data.frame(cmeta)
rownames(cmetadf) = colnames(gcm)
cmetadf$fov = factor(cmetadf$fov)

#' Init. seurat object 
rownames(cmetadf) = colnames(gcm)
seu = CreateSeuratObject(counts = gcm, meta.data = cmetadf)
seu = subset(seu, nCount_RNA > quantile(nCount_RNA, 0.05) & 
                 nCount_RNA < quantile(nCount_RNA, 0.98))
seu = subset(seu, fov %in% c(1:7, 13:16))
seu$fov = droplevels(seu$fov)

sdims = c('x', 'y')
seu@meta.data[, sdims] = seu@meta.data[, c('CenterX_global_px', 'CenterY_global_px')]

#' Non-spatial
seu = NormalizeData(seu, scale.factor = median(colSums(seu@assays$RNA@counts)))
seu = ScaleData(seu, verbose = FALSE)
seu = RunPCA(seu, features = rownames(seu), npcs = 20, verbose = FALSE)
seu = RunHarmony(seu, group.by.vars = 'fov', do_pca = FALSE, verbose = FALSE)
seu = RunUMAP(seu, reduction = 'harmony', reduction.name = 'umap', dims = 1:20, verbose = FALSE, umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
message('Clustering')
seu = FindNeighbors(seu, reduction = 'harmony', dims = 1:20, k.param = 50, verbose = FALSE)
seu = FindClusters(seu, resolution = seq(0.5,2.5,by=0.1), verbose = FALSE, random.seed = SEED)

#' banksy
seu = RunBanksy(seu, lambda = 0.2, M = 1, dimx = sdims[1], dimy = sdims[2], features = 'all', k_geom = 15, assay = 'RNA', slot = 'data', verbose = FALSE)
seu = RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 20, reduction.name = 'pca_banksy', verbose = FALSE)
seu = RunHarmony(seu, group.by.vars = 'fov',  do_pca=FALSE, reduction = 'pca_banksy', reduction.save = 'harmony_banksy', verbose = FALSE)
seu = RunUMAP(seu, reduction = 'harmony_banksy', reduction.name = 'umap_banksy', dims = 1:20, verbose = FALSE, umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
message('Clustering')
seu = FindNeighbors(seu, reduction = 'harmony_banksy', dims = 1:20, k.param = 50, verbose = FALSE) 
seu = FindClusters(seu, resolution = seq(0.5,2.5,by=0.1), verbose = FALSE, random.seed = SEED)

saveRDS(seu, '../out/gsm3683.rds')
