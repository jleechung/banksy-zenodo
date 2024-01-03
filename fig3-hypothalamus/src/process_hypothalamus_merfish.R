# process_merfish_main.R
rm(list=ls())
graphics.off()
out.dir = 'fig3-hypothalamus/out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'fig3-hypothalamus/out/merfish'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)
# You need the (legacy) version of BANKSY (v 0.1.5) for recreating this analysis: 
# remotes::install_github("prabhakarlab/Banksy@main") # to install 0.1.5. 

################## # User inputs # ################## 
# You can either use our provided banksy object (see download instructions below) 
# or cluster the data yourself. Even if you wish to cluster the data yourself, we 
# recommend running the script once with the provided banksy object. 
# 1. We provide a Banksy object (banksyObj_provided.rds) that has had clustering 
# already performed on all naive animals. 
# Downloading this object, and placing it in the correct directory (see below), 
# along with setting the flag:
USE_PROVIDED_BANKSY_OBJ = TRUE
# to TRUE allows for exact reproduction of the results in the paper. 
# Download instructions: 
# The Banksy Object can be downloaded from 
 # the Dropbox link at
# https://www.dropbox.com/scl/fi/eq05c2gip0g61vc0i6n1e/banksyObj_provided.rds?rlkey=z7w2tywtn4h8jiapeymv0l54e&dl=0
# Save this object in the banksy-zenodo/fig3-hypothalamus/data directory. 
# In the Zenodo version of this repo, the banksy object is provided, and does not need to be downloaded. 

# 2. If performing clustering yourself (i.e., not using the Banksy object we provide), 
# then you may want to use fewer than the 11 naive animals in your initial run, 
# as the full 11 animals comprise ~0.5 million cells, and can some time and RAM to 
# run the full UMAP and  clustering steps. 
# As a start, we recommend using just 
# 2 or 4 animals, which lead to qualitatively similar results as the full set of naive animals. 
number_of_animals = 2 # set this to 11 to use all naive animals. 

################# # If running the clustering yourself # ################
# It is possible that the numbers assigned to the clusters are slightly different than what were obtained 
# by us. In some of the plotting and DE gene calling analysis below, we have manually specified these cluster numbers, 
# so to reproduce our results, you will need to determine which cluster corresponds to the mature oligodendrocytes. 
# See the note at line 272. 
# 
################## # Analysis Code # ################## 

library(Banksy)
library(gridExtra) # grid.arrange
library(ggplot2) # facet_wrap 
# library(Seurat)
library(scales) # show_col
library(data.table) # fread
# library(irlba)
# library(Matrix)
# library(tidyverse)
library(plyr) # mapvalues
# library(scran)
library(tictoc)
library(ComplexHeatmap)
# library(dbscan)
# library(circlize)
library(peakRAM)
# library(qvalue) # check how many of these are needed.

results.dir = 'fig3-hypothalamus/out/merfish/'
data.dir = 'fig3-hypothalamus/data/'

k_geom = c(15, 30);lambda = c(0.2);npcs = 20;k_expr = 50;res = seq(0.5, 5, 0.25)
list_of_animal_IDs = 1:number_of_animals 

if (USE_PROVIDED_BANKSY_OBJ){
  bank <- readRDS(file = paste0(data.dir, 'banksyObj_provided.rds'))
  
} else {
  # # Load data
  all_mfish = fread('fig3-hypothalamus/data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv') # see the readme file in the data dir
  all_mfish <- all_mfish[,-c('Fos')]# remove Fos gene per Moffitt manuscript
  all_mfish = cbind(cell_ids = paste0('cell_', 1:nrow(all_mfish)), all_mfish)
  m.list = lapply(list_of_animal_IDs, function(x) all_mfish[all_mfish$Animal_ID==x,])
  expr.list <- lapply(m.list, function(x){
    expr = t(as.matrix(x[,-c(1:10)]))
    colnames(expr) <- x$cell_ids
    expr})
  locs.list <- lapply(m.list, function(x){
    df = as.data.frame(cbind(sdimx = x$Centroid_X, sdimy = x$Centroid_Y,
                             sdimz = 1000*x$Bregma))
    rownames(df) <- x$cell_ids
    df})
  names(expr.list) <- names(locs.list) <- paste0('Animal_', list_of_animal_IDs)
  
  # ------ start clustering block ---------
  bank <- BanksyObject(own.expr = expr.list, cell.locs = locs.list)
  all_mfish_by_animal_list = lapply(list_of_animal_IDs, function(x) all_mfish[Animal_ID %in% x,
                                                                              c('cell_ids',
                                                                                'Animal_ID',
                                                                                'Animal_sex',
                                                                                'Behavior',
                                                                                'Bregma',
                                                                                'Centroid_X',
                                                                                'Centroid_Y',
                                                                                'Cell_class',
                                                                                'Neuron_cluster_ID')])
  bank@meta.data = cbind(bank@meta.data, do.call(rbind, all_mfish_by_animal_list))
  
  # bank <- ComputeBanksy(bank, k_geom = k_geom)
  ram_compute_banksy = peakRAM(bank <- ComputeBanksy(bank, k_geom = k_geom))
  
  # bank <- ScaleBanksy(bank)
  ram_scale_banksy = peakRAM(
    bank <- ScaleBanksy(bank)
  )
  
  lambdas<-c(0, lambda)
  
  # bank <- Banksy:::RunBanksyPCA(bank, lambda = lambdas, npcs = npcs)
  ram_pca_banksy = peakRAM(
    bank <- Banksy:::RunBanksyPCA(bank, lambda = lambdas, npcs = npcs)
  )
  
  # bank <- Banksy:::RunBanksyUMAP(bank, lambda = lambdas, npcs = npcs, nneighbors = k_expr)
  
  ram_umap_banksy = peakRAM(
    bank <- Banksy:::RunBanksyUMAP(bank, lambda = lambdas, npcs = npcs, nneighbors = k_expr)
  )
  
  clust_major_numeric = as.numeric(factor(bank@meta.data$Cell_class))
  clust_major_names = as.character(factor(bank@meta.data$Cell_class))
  bank@meta.data = cbind(bank@meta.data, clust_major_num = clust_major_numeric)
  set.seed(42)
  ptm <- proc.time()
  # bank <- ClusterBanksy(bank, lambda = lambdas, pca = TRUE, npcs = npcs,
  #                       method = 'leiden', k.neighbors = k_expr, resolution = res)
  ram_cluster_banksy = peakRAM(
    bank <- ClusterBanksy(bank, lambda = lambdas, pca = TRUE, npcs = npcs,
                          method = 'leiden', num.cores = 8,
                          k.neighbors = k_expr, resolution = res)
  )
  time_taken = proc.time() - ptm
  print(time_taken)
  saveRDS(bank, file = paste0(data.dir, 'banksyObj_naive_', number_of_animals, '.rds'))
  # bank <- readRDS(file = paste0(data.dir, 'banksyObj_naive_run.rds'))
  # ------ end clustering block ---------
}



reorder_genes <- function(bank){
  if (is.list(bank@own.expr)) {
    x<-bank@own.expr[[1]]
    d_gene <- dist(as.matrix(x))
    hc_gene <- hclust(d_gene)
    bank@own.expr <- lapply(bank@own.expr, function(x) x[hc_gene$order,])
    bank@nbr.expr <- lapply(bank@nbr.expr, function(x) x[hc_gene$order,])
    m1.list <- lapply(bank@harmonics, function(x) list(m1 = x$m1[hc_gene$order,]))
  } else {
    x<-bank@own.expr
    d_gene <- dist(as.matrix(x))
    hc_gene <- hclust(d_gene)
    bank@own.expr <- bank@own.expr[hc_gene$order,]    
    bank@nbr.expr <- bank@nbr.expr[hc_gene$order,] 
    bank@harmonics$m1 <- bank@harmonics$m1[hc_gene$order,]
  }
  return(bank)
}

non_ambig = bank@meta.data$cell_ID[which(bank@meta.data$Cell_class != 'Ambiguous')]
bank = SubsetBanksy(bank, cells = non_ambig)
nonspatial.main.run = 'clust_M1_lam0_k50_res0.5'
spatial.main.run = 'clust_M1_lam0.2_k50_res0.5'
moffitt.labels = 'clust_major_num'
bank = ConnectClusters(bank = bank, map.to = spatial.main.run)
num_clusters<-max(bank@meta.data[,clust.names(bank)])
hypo.cols<-Banksy:::getPalette(num_clusters)
names(hypo.cols)<-1:num_clusters
hypo.cols.old = hypo.cols
scales::show_col(hypo.cols)

bank<- reorder_genes(bank)
set.seed(42);
sampled.cells = bank@meta.data$cell_ID[sample.int(nrow(bank@meta.data),
                                                       min(80000, nrow(bank@meta.data)))]
bank.sampled = SubsetBanksy(bank, cells = sampled.cells)



hypo.cols['3'] = '#848482'

hypo.cols['5'] = '#008856'

hypo.cols['7'] = '#BE0032'
hypo.cols['8'] = '#F38400'

hypo.cols['10'] = '#F99379'
hypo.cols['11'] = '#5F6B6D' 

hypo.cols['13'] = '#0067A5'

hypo.cols['15'] = '#F6A6A0'



bank.runs = c('clust_M1_lam0_k50_res0.5', 
  'clust_M1_lam0.2_k50_res0.5') 
bank.reductions<- paste0('umap_M1_lam', as.numeric(gsub('.*lam|_.*', '', bank.runs)))

png(paste0(results.dir, '/umap_banksy.png'), units = 'in',  res = 200,  
    height = 9*1.25, width = 18*1.25)
umapdims <- mapply(FUN = function(reduction, by, title) plotReduction(bank.sampled, 
                                                                      reduction = reduction, 
                                                                      main = title, legend = TRUE,
                                                                      by = by, col.discrete = hypo.cols,
                                                                      type = 'discrete', 
                                                                      pt.size = 0.2, main.size = 35), 
                   reduction = as.list(bank.reductions), 
                   by = as.list(c(bank.runs[2], bank.runs[2])), 
                   title = as.list(c('Non-spatial UMAP space', 'BANKSY UMAP space')),
                   SIMPLIFY = FALSE)
do.call("grid.arrange", c(umapdims, ncol = length(bank.runs)))
dev.off()



# spatial plot for animal 1. 
bank.animal1 = copy(bank)
animal_id = 'Animal_1'
bank.animal1@own.expr = bank.animal1@own.expr[[animal_id]]
bank.animal1@nbr.expr = bank.animal1@nbr.expr[[animal_id]]
bank.animal1@harmonics = bank.animal1@harmonics[[animal_id]]
bank.animal1@cell.locs = bank.animal1@cell.locs[[animal_id]]
bank.animal1@meta.data = bank.animal1@meta.data[bank.animal1@meta.data$dataset %in% c(animal_id), ]
bank.animal1@reduction$pca_M1_lam0$x <- bank.animal1@reduction$pca_M1_lam0$x[grep(paste0(animal_id, '_'), 
                                                                                  rownames(bank.animal1@reduction$pca_M1_lam0$x) ),
]
bank.animal1@reduction$pca_M1_lam0.2$x <- bank.animal1@reduction$pca_M1_lam0.2$x[grep(paste0(animal_id, '_'), 
                                                                                      rownames(bank.animal1@reduction$pca_M1_lam0.2$x) ),
]
bank.animal1@reduction$umap_M1_lam0 <- bank.animal1@reduction$umap_M1_lam0[grep(paste0(animal_id, '_'), 
                                                                                rownames(bank.animal1@reduction$umap_M1_lam0)),
]
bank.animal1@reduction$umap_M1_lam0.2 <- bank.animal1@reduction$umap_M1_lam0.2[grep(paste0(animal_id, '_'), 
                                                                                    rownames(bank.animal1@reduction$umap_M1_lam0.2)),
]


titles = c('Non-spatial', 'BANKSY')
names(titles) = bank.runs

spatial_plots2 <- function(x, colorscheme = hypo.cols) {
  for (runid in bank.runs){
    png(paste0(results.dir, '/spatial_',runid,'.png'), height = 7, width = 7, 
        units = 'in', res = 300)
    layer_IDs <- unique(bank.animal1@meta.data$Bregma)
    # print(titles[runid])
    spatdims<-vector(mode = "list", length = length(layer_IDs))
    for (jj in 1:length(layer_IDs)){
      layer_id<-layer_IDs[jj]
      bank_layer<- SubsetBanksy(x,
                                cells = x@meta.data$cell_ID[x@meta.data$Bregma %in% layer_id]) # this works...
      spatdims[[jj]]<-plotSpatial(bank_layer, type = 'discrete',
                                  by = runid,
                                  col.discrete = colorscheme,
                                  pt.size = 0.2,legend = FALSE, 
                                  # main = paste0(layer_id, 'mm'),
                                  main.size = 15)+#facet_wrap(~feature, ncol = 5)+ 
        scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
        scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = NULL) +
        theme(axis.title.y=element_text(size = 50))
    }
    do.call("grid.arrange", c(spatdims, ncol = 4))
    dev.off()
  }
}
spatial_plots2(bank.animal1, hypo.cols)




################# Identify DE genes in the oligodendrocyte subpopulations #####################
# In the BANKSY clustering (clust_M1_lam0.2_k50_res0.5), clusters 7 and 8 correspond to the 
# white and grey matter oligodendrocytes. In the nonspatial clustering (clust_M1_lam0_k50_res0.5),
# cluster 7 merges these two clusters. The the numerical cluster IDs could be different if you 
# are not using the provided banksy object, and are instead clustering the data from scratch. 
# If this is the case, you will need to determine the numbers yourself, by either looking 
# at the UMAP, heatmap, or the spatial distributions of the clusters. 
# 
# Nonpspatial: 
# ODM: cluster 7 (red)
#
# Spatial: 
# ODM: clusters 7, 8 (7=red=dense, 8=orange=sparse)

# banksy objects with just ODs
clusname <- c(7) # red cluster (mixed ACA and general hypothalamic region
OD_cells_nonsp <- bank@meta.data[[nonspatial.main.run]] %in% clusname
clusname <- c(7,8) # full od cluster
OD_cells_sp <- bank@meta.data[[spatial.main.run]] %in% clusname
OD_cells<- bank@meta.data$cell_ID[OD_cells_nonsp  & OD_cells_sp]
OD_cells_sp<-NULL
OD_cells_nonsp<- NULL
bank_OD <- SubsetBanksy(bank,cells = OD_cells) #, features = OD.mk
bank_OD_unscaled <- copy(bank_OD)
bank_OD <- ScaleBanksy(bank_OD)

num_de_genes<- 20 
mk_sp_w_od = scran::findMarkers(do.call(cbind, bank_OD@own.expr),
                                groups = bank_OD@meta.data[[spatial.main.run]], 
                                pval.type = 'all', 
                                test.type = 'wilcox')
genes_sp_w_od<-as.data.frame(lapply(mk_sp_w_od, function(x) rownames(x)[1:num_de_genes]))
names(genes_sp_w_od)<-paste0('cluster_', names(mk_sp_w_od))

mk_sp_t_od = scran::findMarkers(do.call(cbind, bank_OD@own.expr),
                                groups = bank_OD@meta.data[[spatial.main.run]],
                                pval.type = 'all',
                                test.type = 't')
genes_sp_t_od<-as.data.frame(lapply(mk_sp_t_od, function(x) rownames(x)[1:num_de_genes]))
names(genes_sp_t_od)<-paste0('cluster_', names(mk_sp_t_od))

# Subset out and plot OD markers and cells
OD.mk.de<-union(unique(unlist(genes_sp_w_od[1:10,])), 
                       unique(unlist(genes_sp_t_od[1:10,])))
# [1] "Mbp"   "Mlc1"  "Cck"   "Gad1"  "Sln"   "Lpar1" "Trh"   "Gjc3"  "Gnrh1" "Ucn3"  "Plin3" "Cbln2" "Syt4"  "Dgkk" 

#### - redo with the top 10 de for wilcox and t test genes: 
bank_OD_mk <- SubsetBanksy(bank, cells = OD_cells, features = OD.mk.de) #
bank_OD_mk_scaled <- ScaleBanksy(bank_OD_mk)
bank_OD_mk_scaled<-reorder_genes(bank_OD_mk_scaled)

pdf(paste0(results.dir, '/hm_OD' ,'.pdf'), 
    height = 4, width = 7)
p.bank.od.scaled_10only <- Banksy::plotHeatmap(bank_OD_mk_scaled,
                                               # features = genes_ord,
                                               assay = 'banksy',
                                               M = 1,
                                               max.cols = 4000,
                                               lambda = 0.2,
                                               cex.row = 6,
                                               annotate = TRUE,
                                               annotate.by = c(nonspatial.main.run, spatial.main.run),
                                               order.by = spatial.main.run,
                                               col.discrete = hypo.cols,
                                               annotation.name = TRUE,
                                               rasterize = TRUE,
                                               cluster_row_slices = FALSE,
                                               cluster_rows = FALSE,
                                               cluster_column_slices = FALSE,
                                               cluster_columns = FALSE)
print(p.bank.od.scaled_10only)
dev.off()
#

# next, plot the spatial plots of the relevant clusters and z slices.
bank.animal1_OD = copy(bank_OD)
animal_id = 'Animal_1'
bank.animal1_OD@own.expr = bank.animal1_OD@own.expr[[animal_id]]
bank.animal1_OD@nbr.expr = bank.animal1_OD@nbr.expr[[animal_id]]
bank.animal1_OD@harmonics = bank.animal1_OD@harmonics[[animal_id]]
bank.animal1_OD@cell.locs = bank.animal1_OD@cell.locs[[animal_id]]
bank.animal1_OD@meta.data = bank.animal1_OD@meta.data[bank.animal1_OD@meta.data$dataset %in% c(animal_id), ]
bank.animal1_OD@reduction$pca_M1_lam0$x <- bank.animal1_OD@reduction$pca_M1_lam0$x[grep(paste0(animal_id, '_'), 
                                                                                        rownames(bank.animal1_OD@reduction$pca_M1_lam0$x) ),
]
bank.animal1_OD@reduction$pca_M1_lam0.2$x <- bank.animal1_OD@reduction$pca_M1_lam0.2$x[grep(paste0(animal_id, '_'), 
                                                                                            rownames(bank.animal1_OD@reduction$pca_M1_lam0.2$x) ),
]

bank.animal1_OD@reduction$umap_M1_lam0 <- bank.animal1_OD@reduction$umap_M1_lam0[grep(paste0(animal_id, '_'), 
                                                                                      rownames(bank.animal1_OD@reduction$umap_M1_lam0)),
]
bank.animal1_OD@reduction$umap_M1_lam0.2 <- bank.animal1_OD@reduction$umap_M1_lam0.2[grep(paste0(animal_id, '_'), 
                                                                                          rownames(bank.animal1_OD@reduction$umap_M1_lam0.2)),
]


runid<-nonspatial.main.run
pdf(paste0(results.dir, '/spatial_OD_',runid, '.pdf'), height = 2.4, width = 3.5)
layer_IDs <- c(0.16, 0.26)
spatdims<-vector(mode = "list", length = length(layer_IDs))
for (i in 1:length(layer_IDs)){
  layer_id<-layer_IDs[i]
  bank_layer<- SubsetBanksy(bank.animal1_OD, 
                            metadata = Bregma %in% layer_id)
  spatdims[[i]]<-plotSpatial(bank_layer, type = 'discrete',
                             by = runid, 
                             # main = runid,
                             legend = FALSE, 
                             col.discrete = hypo.cols,
                             pt.size = 0.35, 
                             main.size = 10)+facet_wrap(~feature)+
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
    scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = NULL) +
    theme(axis.title.y=element_text(size = 50))
}
do.call("grid.arrange", c(spatdims, ncol = 2))
dev.off()

runid<-spatial.main.run
pdf(paste0(results.dir, '/spatial_OD_',runid, '.pdf'),height = 2.4, width = 7)
layer_IDs <- c(0.16, 0.26)
spatdims<-vector(mode = "list", length = length(layer_IDs))
for (i in 1:length(layer_IDs)){
  layer_id<-layer_IDs[i]
  bank_layer<- SubsetBanksy(bank.animal1_OD, 
                            metadata = Bregma %in% layer_id)
  spatdims[[i]]<-plotSpatial(bank_layer, type = 'discrete',
                             by = runid, legend = FALSE, 
                             # main = runid,
                             col.discrete = hypo.cols,
                             pt.size = .35, 
                             main.size = 10)+facet_wrap(~feature)+
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
    scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = NULL) +
    theme(axis.title.y=element_text(size = 50))
}
do.call("grid.arrange", c(spatdims, ncol = 2))
dev.off()

# ------------------------------------------------------------------------------
# Plot metagene plots (ODM)
# Dense subcluster markers: Mbp, Lpar1
# Sparse subclustr markers: Gad1, Mlc1, Cbln2, Syt4
OD.mk.dense<- rownames(bank_OD_mk_scaled@own.expr[[1]])[10:14]#cluster 7
metagene_own1 <- colMeans(do.call(cbind, bank_OD_mk_scaled@own.expr)[OD.mk.dense,])
metagene_nbr1 <- colMeans(do.call(cbind, bank_OD_mk_scaled@nbr.expr)[paste0(OD.mk.dense, '.nbr'),])

OD.mk.sparse<- rownames(bank_OD_mk_scaled@own.expr[[1]])[1:9]#cluster 8
metagene_own2 <- colMeans(do.call(cbind, bank_OD_mk_scaled@own.expr)[OD.mk.sparse,])
metagene_nbr2 <- colMeans(do.call(cbind, bank_OD_mk_scaled@nbr.expr)[paste0(OD.mk.sparse, '.nbr'),])

metagene_own <- metagene_own1-metagene_own2
metagene_nbr <- metagene_nbr1-metagene_nbr2
clust_labels <- bank_OD_mk_scaled@meta.data[[spatial.main.run]]
metagene_df <- data.frame(metagene_own, metagene_nbr, clust_labels )
colormaps = hypo.cols[c('7', '8')]                
cols = colormaps[as.character(as.vector(unlist(clust_labels)))]
OD_means<-aggregate(.~clust_labels, metagene_df, mean)
# 


df <- data.frame(x = sqrt(0.8)*metagene_own, 
                 y = sqrt(0.2)*metagene_nbr, 
                 colors=factor(cols))
p1 <- ggplot(df, aes(x, y, color=colors)) + geom_point(alpha=0.04) +
  scale_color_manual(values=unname(hypo.cols[c( '7', '8' )]))+
  theme_bw() + theme(legend.position="none",
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+
  xlim(-3, 3)+ylim(-1.5, 1.5)
p2<- p1+ geom_point(aes(x=sqrt(0.8)*OD_means$metagene_own[1],
                        y=sqrt(0.2)*OD_means$metagene_nbr[1]),
                    colour="black",
                    pch=4,
                    cex=2)
p3<-p2+geom_point(aes(x=sqrt(0.8)*OD_means$metagene_own[2],
                      y=sqrt(0.2)*OD_means$metagene_nbr[2]),
                  colour="black", pch=4,
                  cex=2)
p4<-ggExtra::ggMarginal(
  p3,
  type = 'density',
  margins = 'both',
  size = 5,
  groupFill = TRUE
)
pdf( paste0(results.dir, '/metagene_OD_',
            spatial.main.run ,
            '.pdf'), height = 6, width = 10)
print(p4)
dev.off()


