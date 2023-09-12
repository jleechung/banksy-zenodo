
#' Load libs -------------------------------------------------------------------

library(SingleR)
library(SingleCellExperiment)
library(Seurat)
library(SeuratWrappers)
library(future)
library(data.table)
library(plyr)
library(matrixStats)
library(BiocParallel)
plan("multicore", workers = 12)
options(future.globals.maxSize = 6000 * 1024^2)

sample_prefix = 'HumanColonCancerPatient1'
sample_prefix = paste0(sample_prefix, '_k30')
seu = readRDS(sprintf("../out/%s_seurat_object.rds", sample_prefix))
sample_prefix = paste0(sample_prefix, '_res0.7')

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

#' Run SingleR -----------------------------------------------------------------

ref = readRDS('../out/singleR_ref.rds')
print(ref)

message('\n Convert to sce')
sce = as.SingleCellExperiment(seu)

message('\n Run singler')

method = 'wilcox'
n = 10
save_path = sprintf('../out/singleR_results_%s_%s.rds', method, n)

message(sprintf('\n Method: %s', method))
message(sprintf('\n Num. genes: %s', n))
message(sprintf('\n Save path: %s', save_path))

pred = SingleR(
    sce,
    ref = ref,
    assay.type.test = 1,
    labels = ref$cl295v11SubFull,
    BPPARAM = BiocParallel::MulticoreParam(28),
    de.method = method,
    de.n = n
)

print(
    table(pred$labels)
)

message('\n Saving')
saveRDS(pred, save_path)
