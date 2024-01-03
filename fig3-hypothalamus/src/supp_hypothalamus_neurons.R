# supp_hypothalamus_neurons.R

rm(list=ls())
graphics.off()

################## # Load Data # ################## 
# Load the clustering result from the file process_hypothalamus_merfish.R
# In that script, users had the option for using our clustering solution 
# (using 'fig3-hypothalamus/data/banksyObj_provided.rds' (see README.md in that directory),
# or of performing the 
# clustering themselves (which takes ~24 hours on a single core). Here, we will 
# use our clustering solution for simplicity, but if users generated their own 
# clustering result, they may use that BANKSY object instead. 
# If using your own saved banksy object, change line 32 accordingly. 



out.dir = 'fig3-hypothalamus/out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'fig3-hypothalamus/out/merfish_supp_neuron'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)
data.dir = 'fig3-hypothalamus/data/'


USE_PROVIDED_BANKSY_OBJ = TRUE # if true, download 
# bank_exc_rd2.rds and bank_inh_rd2.rds using the links in 
# fig3-hypothalamus/data/README.md

bank <- readRDS(file = paste0(data.dir, 'banksyObj_provided.rds'))

################## # Load libraries # ################## 
# 
library(Banksy)
library(gridExtra) # grid.arrange
library(ggplot2) # facet_wrap 
library(scales) # show_col
library(data.table) # fread
library(plyr) # mapvalues
# 
# 
ambig.cell.ids = bank@meta.data$cell_ID[which(bank@meta.data$Cell_class != 'Ambiguous')]
bank = SubsetBanksy(bank, cells = ambig.cell.ids)

npcs = 20; 
k_expr = 50; 

nonspatial.main.run = 'clust_M1_lam0_k50_res2'
spatial.main.run = 'clust_M1_lam0.2_k50_res2'
moffitt.labels = 'clust_major_num'
bank = ConnectClusters(bank = bank, map.to = spatial.main.run)

num_clusters<-max(bank@meta.data[,clust.names(bank)])
hypo.cols<-Banksy:::getPalette(num_clusters)
names(hypo.cols)<-1:num_clusters
hypo.cols.old = hypo.cols
show_col(hypo.cols)


subset_genes <- function(bank, genes){
  # order all datasets' genes based on first dataset
  # gene_names = rownames(bank@own.expr[[1]])
  x<-bank@own.expr[[1]]
  selected = rownames(x) %in% genes 
  bank@own.expr <- lapply(bank@own.expr, function(x) x[selected,])
  bank@nbr.expr <- lapply(bank@nbr.expr, function(x) x[selected,])
  m1.list <- lapply(bank@harmonics, function(x) list(m1 = x$m1[selected,]))
  return(bank)
}
bank.genesubsetted = subset_genes(bank, c('Gad1', 'Slc17a6'))

png(paste0(results.dir, '/','hm_', 'res4' , '.png'), 
    height = 2, width = 102, units = 'in', res = 100)
p.all.mk<-Banksy::plotHeatmap(bank.genesubsetted,
                              assay = 'own.expr',
                              M = 1,
                              max.cols = 12000,
                              # lambda = 0.2,
                              cex.row = 25,
                              annotate = TRUE,
                              annotate.by = c('clust_major_num', 
                                              'clust_M1_lam0_k50_res4', 
                                              'clust_M1_lam0.2_k50_res4'),
                              order.by = 'clust_M1_lam0.2_k50_res4',
                              col.discrete = hypo.cols,
                              annotation.name = TRUE,
                              rasterize = TRUE,
                              cluster_row_slices = FALSE,
                              cluster_rows = FALSE,
                              cluster_column_slices = TRUE,
                              cluster_columns = FALSE)
print(p.all.mk)
dev.off()

exc.rd1.clust = c(7, 22, 23, 35, 37, 40, 45, 46, 48, 54, 56, 61, 63, 68) 
inh.rd1.clust = c(8, 14, 17:21, 26, 27, 30:32, 34, 36, 38, 41,47 ,53, 65)
res.save.str = '_rd1res4_'
clustering.name = 'clust_M1_lam0.2_k50_res4'
bank.exc.banksy = SubsetBanksy(bank, metadata = clust_M1_lam0.2_k50_res4 %in% exc.rd1.clust)
bank.inh.banksy = SubsetBanksy(bank, metadata = clust_M1_lam0.2_k50_res4 %in% inh.rd1.clust)

cn.to.del = clust.names(bank.exc.banksy)[!(clust.names(bank.exc.banksy)%in% c('clust_M1_lam0.2_k50_res4'))]
bank.exc.banksy@meta.data[cn.to.del]<-NULL
bank.inh.banksy@meta.data[cn.to.del]<-NULL


res.neurons.exc = seq(2.1, 2.8, 0.1) 
res.neurons.inh = seq(2.5, 3, 0.1) 

if (USE_PROVIDED_BANKSY_OBJ){
  bank.exc.banksy <- readRDS(file = paste0(data.dir, 'bank_exc_rd2.rds'))
  bank.inh.banksy <- readRDS(file = paste0(data.dir, 'bank_inh_rd2.rds'))
  
} else {
  bank.exc.banksy <- ScaleBanksy(bank.exc.banksy)
  bank.exc.banksy <- Banksy:::RunBanksyPCA(bank.exc.banksy, lambda = c(0, 0.2), npcs = npcs)
  bank.exc.banksy <- Banksy:::RunBanksyUMAP(bank.exc.banksy, lambda = c(0, 0.2),
                                            npcs = npcs, nneighbors = k_expr)
  ptm <- proc.time()
  set.seed(42)
  bank.exc.banksy <- ClusterBanksy(bank.exc.banksy, lambda = c(0, 0.2), pca = TRUE, npcs = npcs,
                                   method = 'leiden', k.neighbors = k_expr, resolution = res.neurons.exc)
  print(proc.time() - ptm)
  saveRDS(bank.exc.banksy, file = paste0(data.dir, '/', 'bank_exc_rd2.rds'))
  # bank.exc.banksy = readRDS(file = paste0(data.dir, '/', 'bank_exc_rd2.rds'))
  
  bank.inh.banksy <- ScaleBanksy(bank.inh.banksy)
  bank.inh.banksy <- Banksy:::RunBanksyPCA(bank.inh.banksy, lambda = c(0, 0.2), npcs = npcs)
  bank.inh.banksy <- Banksy:::RunBanksyUMAP(bank.inh.banksy, lambda = c(0, 0.2), npcs = npcs, nneighbors = k_expr)
  ptm <- proc.time()
  set.seed(42)
  bank.inh.banksy <- ClusterBanksy(bank.inh.banksy, lambda = c(0, 0.2), pca = TRUE, npcs = npcs,
                                   method = 'leiden', k.neighbors = k_expr, resolution = res.neurons.inh)
  print(proc.time() - ptm)
  saveRDS(bank.inh.banksy, file = paste0(data.dir, '/', 'bank_inh_rd2.rds'))
  # bank.inh.banksy = readRDS(file = paste0(data.dir, '/', 'bank_inh_rd2.rds'))
}


exc.labels.ord = paste0('E-', 1:31)
inh.labels.ord = c(paste0('I-', 1:27), 'H-1', paste0('I-', 29:39))  


bank.exc.banksy@meta.data$Neuron_cluster_ID[bank.exc.banksy@meta.data$Neuron_cluster_ID==""] = 'non_neuron'
clust_neuron_numeric = as.numeric(factor(bank.exc.banksy@meta.data$Neuron_cluster_ID,
                                         levels = c(exc.labels.ord, inh.labels.ord, "non_neuron")))

clust_neuron_names = as.character(factor(bank.exc.banksy@meta.data$Neuron_cluster_ID,
                                         levels = c(exc.labels.ord, inh.labels.ord, "non_neuron")))
clust_neuron_named_vec = clust_neuron_numeric
names(clust_neuron_named_vec) = clust_neuron_names
bank.exc.banksy@meta.data$clust_neuron_num = clust_neuron_numeric

bank.exc.banksy = ConnectClusters(bank = bank.exc.banksy, map.to = 'clust_neuron_num')
num_clusters<-max(bank.exc.banksy@meta.data[,clust.names(bank.exc.banksy)])
Exc.cols<-Banksy:::getPalette(num_clusters)
names(Exc.cols)<-1:num_clusters
Exc.cols.old = Exc.cols
show_col(Exc.cols)

bank.inh.banksy@meta.data$Neuron_cluster_ID[bank.inh.banksy@meta.data$Neuron_cluster_ID==""] = 'non_neuron'
clust_neuron_numeric = as.numeric(factor(bank.inh.banksy@meta.data$Neuron_cluster_ID, 
                                         levels = c(inh.labels.ord, exc.labels.ord, "non_neuron")))

clust_neuron_names = as.character(factor(bank.inh.banksy@meta.data$Neuron_cluster_ID,  
                                         levels = c(inh.labels.ord, exc.labels.ord, "non_neuron")))
clust_neuron_named_vec = clust_neuron_numeric
names(clust_neuron_named_vec) = clust_neuron_names
bank.inh.banksy@meta.data$clust_neuron_num = clust_neuron_numeric

bank.inh.banksy = ConnectClusters(bank = bank.inh.banksy, map.to = 'clust_neuron_num')
num_clusters<-max(bank.inh.banksy@meta.data[,clust.names(bank.inh.banksy)])
Inh.cols<-Banksy:::getPalette(num_clusters)
names(Inh.cols)<-1:num_clusters
Inh.cols.old = Inh.cols
show_col(Inh.cols)

# !! change these to the appropriate names:
exc.nonsp.clustering = 'clust_M1_lam0_k50_res2.2'
exc.bank.clustering = 'clust_M1_lam0.2_k50_res2.2'
inh.nonsp.clustering = 'clust_M1_lam0_k50_res2.8'
inh.bank.clustering = 'clust_M1_lam0.2_k50_res3'
ANIMAL_IDS = c('Animal_1'
  ,'Animal_1'
  ,'Animal_1'
  ,'Animal_1'
)
CLUSTERINGS = c(
  exc.nonsp.clustering 
  ,exc.bank.clustering
  ,inh.nonsp.clustering
  ,inh.bank.clustering
)
# 
CLUST.LABELS.UNIQUE = list(
  sort(unique(bank.exc.banksy@meta.data[[exc.nonsp.clustering]]))
  ,sort(unique(bank.exc.banksy@meta.data[[exc.bank.clustering]]))
  ,sort(unique(bank.inh.banksy@meta.data[[inh.nonsp.clustering]]))
  ,sort(unique(bank.inh.banksy@meta.data[[inh.bank.clustering]]))
)

ROW.LABELS = lapply(CLUST.LABELS.UNIQUE, function(x) 1:length(x))
ROW.LABELS = mapply(function(x, y) {
  names(y) = as.character(x)
  return(y)
  }, CLUST.LABELS.UNIQUE, ROW.LABELS)
names(ROW.LABELS)<-CLUSTERINGS
SAVE.STRINGS = c(
  'exc_nonsp'
  ,
  'exc_banksy'
  ,
  'inh_nonsp'
  ,
  'inh_banksy'
)


# Plot the different clusters for each layer, for animal 1: 
PLOT.OBJS = mapply(function(animal_id, clustering.name, clust.labels.ord, save.string){
  if (isTRUE(grepl('inh', save.string))){
    bank_layer.neur = copy(bank.inh.banksy)
  } else if (isTRUE(grepl('exc', save.string))) {
    bank_layer.neur = copy(bank.exc.banksy)
  }
  bank_layer.neur@own.expr = bank_layer.neur@own.expr[[animal_id]]
  bank_layer.neur@nbr.expr = bank_layer.neur@nbr.expr[[animal_id]]
  bank_layer.neur@harmonics = bank_layer.neur@harmonics[[animal_id]]
  bank_layer.neur@cell.locs = bank_layer.neur@cell.locs[[animal_id]]
  bank_layer.neur@meta.data = bank_layer.neur@meta.data[bank_layer.neur@meta.data$dataset %in% c(animal_id), ]
  bank_layer.neur@reduction$pca_M1_lam0$x <- bank_layer.neur@reduction$pca_M1_lam0$x[grep(animal_id, 
                                                                                          rownames(bank_layer.neur@reduction$pca_M1_lam0$x) ),
  ]
  bank_layer.neur@reduction$pca_M1_lam0.2$x <- bank_layer.neur@reduction$pca_M1_lam0.2$x[grep(animal_id, 
                                                                                              rownames(bank_layer.neur@reduction$pca_M1_lam0.2$x) ),
  ]
  
  bank_layer.neur@reduction$umap_M1_lam0 <- bank_layer.neur@reduction$umap_M1_lam0[grep(animal_id, 
                                                                                        rownames(bank_layer.neur@reduction$umap_M1_lam0)),
  ]
  bank_layer.neur@reduction$umap_M1_lam0.2 <- bank_layer.neur@reduction$umap_M1_lam0.2[grep(animal_id, 
                                                                                            rownames(bank_layer.neur@reduction$umap_M1_lam0.2)),
  ]
  nested.plot.list = lapply(unique(bank_layer.neur@meta.data$Bregma), function(y){
    bank_layer.neur1 = SubsetBanksy(bank_layer.neur, metadata = Bregma %in% !!y)
    print(y)
    inner.plot.list = lapply(clust.labels.ord, function(x){
      print(x)
      selected.cells = rep('unselected', times = nrow(meta.data(bank_layer.neur1)))
      selected.cells[bank_layer.neur1@meta.data[[clustering.name]] %in% x]='selected'
      bank_layer.neur2 = copy(bank_layer.neur1)
      bank_layer.neur2@meta.data$clust_selected = selected.cells
      pal <- c('grey', 'red')
      names(pal) <- c('unselected', 'selected')
      rowlabel = ROW.LABELS[[clustering.name]][as.character(x)]
      if (y == unique(bank_layer.neur@meta.data$Bregma)[1]){
        
        plotSpatial(bank_layer.neur2, type = 'discrete', 
                    by = 'clust_selected',legend = FALSE,
                    col.discrete = pal, 
                    pt.size = 2
        ) + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
          scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = paste0('cl_', rowlabel)) +
          theme(axis.title.y=element_text(size = 50))
      } else {
        plotSpatial(bank_layer.neur2, type = 'discrete', 
                    by = 'clust_selected',legend = FALSE,
                    col.discrete = pal, 
                    pt.size = 1.5
        ) + 
          scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
          scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = NULL)
      }
    })
  })
  unnested.plot.list = unlist(nested.plot.list, recursive = FALSE)
  png(paste0(results.dir, '/','spatial_', '_', save.string, '_', clustering.name , '.png'), 
      res = 20, units = 'in', height = length(unnested.plot.list)/length(nested.plot.list)*2.75, 
      width = length(nested.plot.list)*2.5)
  plot.obj = gridExtra::grid.arrange(grobs = unnested.plot.list, 
                                     layout_matrix = matrix(data = 1:length(unnested.plot.list), 
                                                            byrow = FALSE, 
                                                            ncol = length(nested.plot.list))
  )
  print(plot.obj)
  dev.off()
  return(plot.obj)
}, ANIMAL_IDS, CLUSTERINGS, CLUST.LABELS.UNIQUE, SAVE.STRINGS)

animal_id = 'Animal_1'
inh.fig5.clustering = inh.bank.clustering
exc.fig5.clustering = exc.bank.clustering

############################################


bank.inh.banksy.relabeled = copy(bank.inh.banksy)
bank.inh.banksy.relabeled@own.expr = bank.inh.banksy.relabeled@own.expr[[animal_id]]
bank.inh.banksy.relabeled@nbr.expr = bank.inh.banksy.relabeled@nbr.expr[[animal_id]]
bank.inh.banksy.relabeled@harmonics = bank.inh.banksy.relabeled@harmonics[[animal_id]]
bank.inh.banksy.relabeled@cell.locs = bank.inh.banksy.relabeled@cell.locs[[animal_id]]
bank.inh.banksy.relabeled@meta.data = bank.inh.banksy.relabeled@meta.data[bank.inh.banksy.relabeled@meta.data$dataset %in% c(animal_id), ]
bank.inh.banksy.relabeled@reduction$pca_M1_lam0$x <- bank.inh.banksy.relabeled@reduction$pca_M1_lam0$x[grep(animal_id, 
                                                                                                            rownames(bank.inh.banksy.relabeled@reduction$pca_M1_lam0$x) ),
]
bank.inh.banksy.relabeled@reduction$pca_M1_lam0.2$x <- bank.inh.banksy.relabeled@reduction$pca_M1_lam0.2$x[grep(animal_id, 
                                                                                                                rownames(bank.inh.banksy.relabeled@reduction$pca_M1_lam0.2$x) ),
]

bank.inh.banksy.relabeled@reduction$umap_M1_lam0 <- bank.inh.banksy.relabeled@reduction$umap_M1_lam0[grep(animal_id, 
                                                                                                          rownames(bank.inh.banksy.relabeled@reduction$umap_M1_lam0)),
]
bank.inh.banksy.relabeled@reduction$umap_M1_lam0.2 <- bank.inh.banksy.relabeled@reduction$umap_M1_lam0.2[grep(animal_id, 
                                                                                                              rownames(bank.inh.banksy.relabeled@reduction$umap_M1_lam0.2)),
]


list.of.labels = c(3, 5, 12, 13)*1000
ii = inh.fig5.clustering
bank.inh.banksy.relabeled@meta.data[[ii]] = bank.inh.banksy.relabeled@meta.data[[ii]]*1000
inh.cols.fig5 = c('#0067A5', '#F38400' , '#BE0032', '#875692')
names(inh.cols.fig5) = list.of.labels
plot.cols.inh = rep('#E5E7E9', length(Inh.cols))
names(plot.cols.inh) = paste0(names(Inh.cols), '000')
plot.cols.inh[names(plot.cols.inh) %in% names(inh.cols.fig5)] = inh.cols.fig5
show_col(plot.cols.inh)


bank.exc.banksy.relabeled = copy(bank.exc.banksy)
bank.exc.banksy.relabeled@own.expr = bank.exc.banksy.relabeled@own.expr[[animal_id]]
bank.exc.banksy.relabeled@nbr.expr = bank.exc.banksy.relabeled@nbr.expr[[animal_id]]
bank.exc.banksy.relabeled@harmonics = bank.exc.banksy.relabeled@harmonics[[animal_id]]
bank.exc.banksy.relabeled@cell.locs = bank.exc.banksy.relabeled@cell.locs[[animal_id]]
bank.exc.banksy.relabeled@meta.data = bank.exc.banksy.relabeled@meta.data[bank.exc.banksy.relabeled@meta.data$dataset %in% c(animal_id), ]
bank.exc.banksy.relabeled@reduction$pca_M1_lam0$x <- bank.exc.banksy.relabeled@reduction$pca_M1_lam0$x[grep(animal_id, 
                                                                                                            rownames(bank.exc.banksy.relabeled@reduction$pca_M1_lam0$x) ),
]
bank.exc.banksy.relabeled@reduction$pca_M1_lam0.2$x <- bank.exc.banksy.relabeled@reduction$pca_M1_lam0.2$x[grep(animal_id, 
                                                                                                                rownames(bank.exc.banksy.relabeled@reduction$pca_M1_lam0.2$x) ),
]

bank.exc.banksy.relabeled@reduction$umap_M1_lam0 <- bank.exc.banksy.relabeled@reduction$umap_M1_lam0[grep(animal_id, 
                                                                                                          rownames(bank.exc.banksy.relabeled@reduction$umap_M1_lam0)),
]
bank.exc.banksy.relabeled@reduction$umap_M1_lam0.2 <- bank.exc.banksy.relabeled@reduction$umap_M1_lam0.2[grep(animal_id, 
                                                                                                              rownames(bank.exc.banksy.relabeled@reduction$umap_M1_lam0.2)),
]



options(scipen=999)
list.of.labels = c(3,9)*1000000
ii = exc.fig5.clustering
bank.exc.banksy.relabeled@meta.data[[ii]] = bank.exc.banksy.relabeled@meta.data[[ii]]*1000000
exc.cols.fig5 = c('#008856', '#882D17' )
names(exc.cols.fig5) = list.of.labels
plot.cols.exc = rep('#E5E7E9', length(Exc.cols))
names(plot.cols.exc) = paste0(names(Exc.cols), '000000')
plot.cols.exc[names(plot.cols.exc) %in% names(exc.cols.fig5)] = exc.cols.fig5
show_col(plot.cols.exc)
plot.cols = c(plot.cols.inh, plot.cols.exc)
# now merge the two objects. 
own.merged = cbind(bank.inh.banksy.relabeled@own.expr,
                   bank.exc.banksy.relabeled@own.expr)
nbr.merged = cbind(bank.inh.banksy.relabeled@nbr.expr,
                   bank.exc.banksy.relabeled@nbr.expr)
har.merged = cbind(bank.inh.banksy.relabeled@harmonics$m1, 
                   bank.exc.banksy.relabeled@harmonics$m1)
cell.locs.merged = rbind(bank.inh.banksy.relabeled@cell.locs, 
                         bank.exc.banksy.relabeled@cell.locs)
metadata.merged.temp = rbind(bank.inh.banksy.relabeled@meta.data[,1:12], 
                             bank.exc.banksy.relabeled@meta.data[,1:12])
ii = inh.fig5.clustering
jj = exc.fig5.clustering
metadata.merged = cbind(metadata.merged.temp, 
                        clust_banksy = as.character(c(bank.inh.banksy.relabeled@meta.data[[ii]], 
                                                      bank.exc.banksy.relabeled@meta.data[[jj]])))
reduction.nonsp.merged = rbind(bank.inh.banksy.relabeled@reduction$umap_M1_lam0, 
                               bank.exc.banksy.relabeled@reduction$umap_M1_lam0)
reduction.bank.merged = rbind(bank.inh.banksy.relabeled@reduction$umap_M1_lam0.2, 
                              bank.exc.banksy.relabeled@reduction$umap_M1_lam0.2)

bank.merged = copy(bank.inh.banksy.relabeled)
bank.merged@own.expr = own.merged
bank.merged@nbr.expr = nbr.merged
bank.merged@harmonics$m1 = har.merged
bank.merged@cell.locs = cell.locs.merged
bank.merged@meta.data = metadata.merged
bank.merged@reduction$umap_M1_lam0 = reduction.nonsp.merged
bank.merged@reduction$umap_M1_lam0.2 = reduction.bank.merged
bank.merged@reduction$pca_M1_lam0 <-NULL
bank.merged@reduction$pca_M1_lam0.2 <-NULL

labels.to.plot.top = c('3000', '5000', '3000000')
labels.to.plot.bottom = c('12000', '13000', '9000000')

plot.cols.5a.top = plot.cols
plot.cols.5a.top[labels.to.plot.bottom] = c('#E5E7E9','#E5E7E9','#E5E7E9')

plot.cols.5a.bottom = plot.cols
plot.cols.5a.bottom[labels.to.plot.top] = c('#E5E7E9','#E5E7E9','#E5E7E9')

ptsize = 1.2

nested.plot.list = lapply(unique(bank.merged@meta.data$Bregma)[c(1,3,5,7,9,11)], function(y){
  bank_layer.neur1 = SubsetBanksy(bank.merged,
                                  metadata = Bregma %in% !!y)
  
  data.test1 <- Banksy:::getLocations(bank_layer.neur1)
  sdimx <- sdimy <- NULL
  feature1 <- Banksy:::getFeature(bank_layer.neur1,
                                  by = 'clust_banksy', dataset = NULL)
  
  data.test1 <- cbind(data.test1, feature = as.factor(feature1), 
                      tosort = as.numeric(feature1 %in% labels.to.plot.top))
  sortorder.inh = order(data.test1$tosort)
  
  data.test1 <- data.test1[sortorder.inh,]
  plot <- ggplot(data.test1, aes(x = sdimx, y = sdimy, col = as.factor(feature1[sortorder.inh]))) +
    scale_color_manual(values = Banksy:::getDiscretePalette(feature1[sortorder.inh],
                                                            plot.cols.5a.top)) #+ 
  plot <- plot + geom_point(size = 1, alpha = 0.8) +theme_bw()+ggtitle(paste0(y, 'mm'))
  
  plot<-plot +theme_bw()+theme(legend.position="none")+    scale_x_discrete(labels = NULL, 
                                                                            breaks = NULL) + labs(x = NULL) +
    scale_y_discrete(labels = NULL, 
                     breaks = NULL) + labs(y = NULL)
  plot
})
pdf(paste0(results.dir, '/','spatial_fig5_toprow','.pdf'),
    height = 2.5,
    width = length(nested.plot.list)*2.25)
plot.obj = gridExtra::grid.arrange(grobs = nested.plot.list,  ncol = length(nested.plot.list))
print(plot.obj)
dev.off()

nested.plot.list = lapply(unique(bank.merged@meta.data$Bregma)[c(1,3,5,7,9,11)], function(y){
  bank_layer.neur1 = SubsetBanksy(bank.merged,
                                  metadata = Bregma %in% !!y)

  data.test1 <- Banksy:::getLocations(bank_layer.neur1)
  sdimx <- sdimy <- NULL
  feature1 <- Banksy:::getFeature(bank_layer.neur1,
                                  by = 'clust_banksy', dataset = NULL)
  
  data.test1 <- cbind(data.test1, feature = as.factor(feature1), 
                      tosort = as.numeric(feature1 %in% labels.to.plot.bottom))
  sortorder.inh = order(data.test1$tosort)
  
  data.test1 <- data.test1[sortorder.inh,]
  plot <- ggplot(data.test1, aes(x = sdimx, y = sdimy, col = as.factor(feature1[sortorder.inh]))) +
    scale_color_manual(values = Banksy:::getDiscretePalette(feature1[sortorder.inh],
                                                            plot.cols.5a.bottom)) #+ 
  plot <- plot + geom_point(size = 1, alpha = 0.8) +theme_bw()+ggtitle(paste0(y, 'mm'))
  
  plot<-plot +theme_bw()+theme(legend.position="none")+    scale_x_discrete(labels = NULL, 
                                                                            breaks = NULL) + labs(x = NULL) +
    scale_y_discrete(labels = NULL, 
                     breaks = NULL) + labs(y = NULL)
  plot
})
pdf(paste0(results.dir, '/','spatial_fig5_bottom', '.pdf'),
    height = 2.5,
    width = length(nested.plot.list)*2.25)
plot.obj = gridExtra::grid.arrange(grobs = nested.plot.list,  ncol = length(nested.plot.list))
print(plot.obj)
dev.off()

