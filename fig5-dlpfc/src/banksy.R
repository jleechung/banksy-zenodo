rm(list=ls())
graphics.off()

#' Parse args ------------------------------------------------------------------

library(optparse)

option_list = list(
	make_option(c("--sampleid"), type="numeric"),
    make_option(c("--seed"), type="numeric", default = 1000)
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$sampleid)) {
	print_help(opt_parser)
	stop()
} else {
	print(opt)
}

i = opt$sampleid
SEED = opt$seed

library(repr)
sethw = function(h,w) options(repr.plot.height=h, repr.plot.width=w)

#' Load libs -------------------------------------------------------------------
suppressPackageStartupMessages({

    library(Seurat)
    library(Banksy)
    library(data.table)
    library(aricode)
    library(mltools)

})
    
#' Params
K_GEOM = c(18, 18)
LAM = 0.2
all_samples = as.character(c(151507:151510, 151669:151676))
all_domains = c(rep(7, 4), rep(5, 4), rep(7, 4))
sample_size = length(all_samples)

#' Paths
out = sprintf('../out/seed%s', SEED)
message('Output directory:', out)
if (!file.exists(out)) dir.create(out, recursive = TRUE)
data_path = '../data/'

#' Load data -------------------------------------------------------------------

#' Read data
curr_sample = all_samples[i]
curr_domain = all_domains[i]
message(curr_sample)
curr_sample_files = list.files(data_path, full.names = TRUE, pattern = as.character(curr_sample))
gcm = readRDS(curr_sample_files[grep('exp', curr_sample_files)])
locs = readRDS(curr_sample_files[grep('loc', curr_sample_files)])
anno = readRDS(curr_sample_files[grep('man', curr_sample_files)])
rid = which(is.na(anno))
anno = anno[-rid]
gcm = gcm[,-rid]
locs = locs[-rid,]
lvl = c(sort(unique(anno)))
anno = factor(anno, levels = sort(unique(anno)))
anno = as.numeric(anno)

scale.factor = median(colSums(gcm))

#' Run BANKSY HVG --------------------------------------------------------------
    
seu = CreateSeuratObject(counts = gcm)
seu = NormalizeData(seu, normalization.method = 'RC', scale.factor = scale.factor, verbose = FALSE)
seu = FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)

bank = BanksyObject(own.expr = gcm, cell.locs = locs, meta.data = data.frame(clust_anno = anno))
bank = NormalizeBanksy(bank, norm_factor = scale.factor)
#' Subset after normalization
bank = SubsetBanksy(bank, features = VariableFeatures(seu))
bank = ComputeBanksy(bank, k_geom = K_GEOM, verbose=FALSE)
bank = ScaleBanksy(bank)
bank = RunBanksyPCA(bank, lambda = LAM)
bank = ClusterBanksy(bank, lambda = LAM, resolution = seq(0.2,1.2,0.02), verbose = FALSE, seed = SEED, num.cores = 10)

hits = names(which(
    apply(bank@meta.data[,clust.names(bank)], 2, function(x) length(unique(x))) == curr_domain
))            
if (length(hits) == 1) {
    diff = apply(bank@meta.data[,clust.names(bank)[-1]], 2, function(x) abs(length(unique(x)) - curr_domain))
    hits = c('clust_anno', names(which(diff == min(diff))))
}
                 
bank@meta.data = bank@meta.data[, c('cell_ID', 'nCount', 'NODG', hits)]
bank = SmoothLabels(bank, cluster_names = clust.names(bank)[-1], k=6, verbose=FALSE)
bank = ConnectClusters(bank, map.to = 'clust_anno')
hits = c('clust_anno', clust.names(bank)[grepl('smooth', clust.names(bank))])

#' Compute metrics --------------------------------------------------------------

getMetric <- function(bank, digits, func) {
    clust <- bank@meta.data[, clust.names(bank), drop = FALSE]
    n.clust <- ncol(clust)
    if (n.clust < 2) {
        stop('ARI will only be calculated for at least 2 clustering runs.')
    }
    comb <- combn(names(clust), 2)
    n.comb <- ncol(comb)
    ari <- numeric(length = ncol(comb))
    for (i in seq_len(n.comb)) {
        ari[i] <- func(clust[[comb[1, i]]], clust[[comb[2, i]]])
    }

    ari.mat <- diag(nrow = n.clust, ncol = n.clust)
    ari.mat[lower.tri(ari.mat)] <- ari
    rownames(ari.mat) <- colnames(ari.mat) <- colnames(clust)
    ari.mat <-
        ari.mat + t(ari.mat) - diag(nrow = n.clust, ncol = n.clust)
    ari.mat <- round(ari.mat, digits)
    return(ari.mat)
}
                 
ari = median(getMetric(bank, digits = 5, func = aricode::ARI)[,1][hits][-1])
message(sprintf('Seed=%s, Sample=%s, %s=%s', SEED, curr_sample, 'ARI', ari))
saveRDS(ari, sprintf('../out/%s/%s_%s.rds', out, 'ari', curr_sample))

nmi = median(getMetric(bank, digits = 5, func = aricode::NMI)[,1][hits][-1])
message(sprintf('Seed=%s, Sample=%s, %s=%s', SEED, curr_sample, 'NMI', nmi))
saveRDS(nmi, sprintf('../out/%s/%s_%s.rds', out, 'nmi', curr_sample))

mcc = median(getMetric(bank, digits = 5, func = mltools::mcc)[,1][hits][-1])
message(sprintf('Seed=%s, Sample=%s, %s=%s', SEED, curr_sample, 'MCC', mcc))
saveRDS(mcc, sprintf('../out/%s/%s_%s.rds', out, 'mcc', curr_sample))

bank@meta.data = bank@meta.data[, hits]
saveRDS(bank@meta.data, sprintf('../out/%s/metadata_%s.rds', out, curr_sample))
