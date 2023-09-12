
rm(list=ls())

library(data.table)
library(Seurat)
library(SeuratWrappers) 
library(Matrix)
library(SingleR)
library(BiocParallel)

SEED = 1000

pal = c(pals::kelly()[-c(1,3)], pals::glasbey(n=32))

clip = function(x, low=0.02, high=0.98) {
    qt = quantile(x, c(low, high))
    x[x < qt[1]] = qt[1]
    x[x > qt[2]] = qt[2]
    x
}

ulen = function(x) length(unique(x))

seu = readRDS('../out/gsm3683.rds')

query = GetAssayData(seu, slot = 'data', assay = 'RNA')
dim(query)

# Load 
ref = Matrix::readMM('../data/Colon-Healthy-adult.mtx')
genes = fread('../data/Colon-Healthy-adult-genes.csv.gz')
metadata = fread('../data/Colon-Healthy-adult-metadata.csv.gz')

# Get commmon genes
ref = t(ref)
genes = unlist(genes)
names(genes) = NULL
common_genes = sort(intersect(genes, rownames(query)))
rownames(ref) = genes
ref = ref[common_genes,]
query = query[common_genes,]

head(ref[,1:5])
head(query[,1:5])
print('---------------------')
tail(ref[,1:5])
tail(query[,1:5])

labels_coarse = metadata$category
table(labels_coarse)
labels_fine = metadata$Integrated_05
table(labels_fine)

pred_coarse = SingleR(
    query,
    ref = ref,
    assay.type.test = 1,
    labels = labels_coarse,
    BPPARAM = BiocParallel::MulticoreParam(30),
    de.method = 'wilcox',
    de.n = 30
)

pred_fine = SingleR(
    query,
    ref = ref,
    assay.type.test = 1,
    labels = labels_fine,
    BPPARAM = BiocParallel::MulticoreParam(30),
    de.method = 'wilcox',
    de.n = 15
)

saveRDS(pred_coarse, '../out/pred_coarse.rds')
saveRDS(pred_fine, '../out/pred_fine.rds')
