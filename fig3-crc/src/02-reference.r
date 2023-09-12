
library(Seurat)
library(data.table)
library(repr)
library(ggplot2)
library(gridExtra)
library(future)
library(SingleCellExperiment)
plan("multicore", workers = 14)
plan()
options(future.globals.maxSize = 10000 * 1024^2)
sethw = function(h,w) options(repr.plot.height=h, repr.plot.width=w)
pal = pals::kelly()[-c(1:3)]
pall = Giotto:::getDistinctColors(100)

gc()
expr = Read10X_h5('../data/GSE178341_crc10x_full_c295v4_submit.h5')
gc()

clust = fread('../data/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz')
mdata = fread('../data/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz')
mdata = cbind(mdata, clust[,-1])
gc()

mdata = data.frame(mdata)
rownames(mdata) = colnames(expr)

## Subset the data 
kid = grepl('Regev', mdata$TISSUE_PROCESSING_TEAM)
expr = expr[,kid]
mdata = mdata[kid,]

## Create Seurat
seu = CreateSeuratObject(counts = expr, meta.data = mdata)
gc()

seu = NormalizeData(seu, scale.factor = median(colSums(expr)))
gc()

pbulk = AverageExpression(seu, group.by = 'PID')
pbulk = pbulk$RNA
gc()

## Compute correlation to Vizgen
vizgen = readRDS('../out/HumanColonCancerPatient1_pseudobulk.rds')
colnames(vizgen) = 'vizgen'

common = intersect(rownames(pbulk), rownames(vizgen))
length(common)

pbulk = pbulk[common,]
vizgen = vizgen[common,,drop=FALSE]

comb = cbind(vizgen=vizgen, pbulk)
cor.mat = cor(comb, method = 'spearman')
high_cor_samples = names(which(cor.mat[-1,1] > 0.55))
keep_cells = rownames(mdata)[mdata$PID %in% high_cor_samples]
expr = expr[,keep_cells]
mdata = mdata[keep_cells,]

## Create Seurat
seu = CreateSeuratObject(counts = expr, meta.data = mdata)
seu = NormalizeData(seu, scale.factor = median(colSums(expr)))
seu = FindVariableFeatures(seu, nfeatures = 2000)
seu = ScaleData(seu)
seu = RunPCA(seu, features = VariableFeatures(seu), verbose=FALSE)
gc()

seu$cluster = paste0(seu$SPECIMEN_TYPE, seu$cl295v11SubFull, ' ')

sce = as.SingleCellExperiment(seu)
common = intersect(rownames(sce), rownames(vizgen))
sce = sce[common,]
saveRDS(sce, '../out/singleR_ref.rds')
