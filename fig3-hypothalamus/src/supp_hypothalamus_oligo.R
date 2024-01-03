# supp_hypothalamus_oligo.R

rm(list=ls())
graphics.off()

################## # Load Data # ################## 
# Load the clustering result from the file process_hypothalamus_merfish.R
# In that script, users had the option for using our clustering solution 
# (using 'fig3-hypothalamus/data/banksyObj_provided.rds', or of performing the 
# clustering themselves). Here, we will 
# use our clustering solution for simplicity, but if users generated their own 
# clustering result, they may use that BANKSY object instead, by changing line 23 
# to load the the correct saved object. 

out.dir = 'fig3-hypothalamus/out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'fig3-hypothalamus/out/merfish_supp_oligo'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)
data.dir = 'fig3-hypothalamus/data/'
bank.conn <- readRDS(file = paste0(data.dir, 'banksyObj_provided.rds'))

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


# with ambiguous cells removed: (maybe want to rerun without the ambigs in the first place...)
ambig.cell.ids = bank.conn@meta.data$cell_ID[which(bank.conn@meta.data$Cell_class != 'Ambiguous')]
bank.conn = SubsetBanksy(bank.conn, cells = ambig.cell.ids)


nonspatial.main.run = 'clust_M1_lam0_k50_res2'
spatial.main.run = 'clust_M1_lam0.2_k50_res2'
moffitt.labels = 'clust_major_num'
bank.conn = ConnectClusters(bank = bank.conn, map.to = spatial.main.run)

num_clusters<-max(bank.conn@meta.data[,clust.names(bank.conn)])
mOD.cols<-Banksy:::getPalette(num_clusters)
names(mOD.cols)<-1:num_clusters
mOD.cols.old = mOD.cols
show_col(mOD.cols)


bank.animal1 = copy(bank.conn)
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

num_of_clust = unlist(lapply(bank.animal1@meta.data[clust.names(bank.animal1)], function(x) length(unique(x))))



png(paste0(results.dir, '/','UMAPs_conn_anim1_fewres' , '.png'), res = 150, units = 'in', height = 1.6*10, width = 1.6*20)
grid.arrange(
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0', by = 'clust_major_num', type = 'discrete',
                main = paste0('Moffitt et al. annots, \nnum. clust. = ', num_of_clust[names(num_of_clust)=='clust_major_num']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res0.5', type = 'discrete',
                main = paste0('Nonspatial clusters, \nres = 0.5, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0_k50_res0.5']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res2', type = 'discrete',
                main = paste0('Nonspatial clusters, \nres = 2, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0_k50_res2']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res5', type = 'discrete',
                main = paste0('Nonspatial clusters, \nres = 5, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0_k50_res5']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0.2', by = 'clust_major_num', type = 'discrete',
                main = paste0('Moffitt et al. annots, \nnum. clust. = ', num_of_clust[names(num_of_clust)=='clust_major_num']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res0.5', type = 'discrete',
                main = paste0('BANKSY clusters, \nres = 0.5, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0.2_k50_res0.5']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res2', type = 'discrete',
                main = paste0('BANKSY clusters, \nres = 2, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0.2_k50_res2']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.animal1, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res5', type = 'discrete',
                main = paste0('BANKSY clusters, \nres = 5, num. clust. = ', num_of_clust[names(num_of_clust)=='clust_M1_lam0.2_k50_res5']), 
                main.size = 40, pt.size = 0.15, col.discrete = mOD.cols, legend = FALSE),
  ncol = 4
)
dev.off()


odm.names = grep('OD Mature', unique(bank.animal1@meta.data$Cell_class), value = TRUE)
set.seed(42)
# 
# 
bank.mOD = SubsetBanksy(bank.animal1, metadata = clust_M1_lam0.2_k50_res2 %in% c(8, 9, 16, 32, 38))
bank.mOD = SubsetBanksy(bank.mOD, metadata = clust_M1_lam0.2_k50_res0.5 %in% c(8,9))
bank.mOD = SubsetBanksy(bank.mOD, metadata = Cell_class %in% odm.names)


png(paste0(results.dir, '/','ODM_UMAPs_conne_anim1' , '.png'), res = 150, units = 'in', height = 1.6*10, width = 1.6*20)
grid.arrange(
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0', by = 'clust_major_num', type = 'discrete',
                main = 'Moffitt et al. annots', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res0.5', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 0.5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res2', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 2', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0', by = 'clust_M1_lam0_k50_res5', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_major_num', type = 'discrete',
                main = 'Moffitt et al. annots', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res0.5', type = 'discrete',
                main = 'BANKSY clusters, \nres = 0.5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res2', type = 'discrete',
                main = 'BANKSY clusters, \nres = 2', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0.2_k50_res5', type = 'discrete',
                main = 'BANKSY clusters, \nres = 5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  ncol = 4
)
dev.off()


png(paste0(results.dir, '/','ODM_UMAPs_swap_conne_anim1' , '.png'), res = 150, units = 'in', height = 1.6*5, width = 1.6*15)
grid.arrange(
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0_k50_res0.5', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 0.5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0_k50_res2', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 2', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  plotReduction(bank.mOD, reduction = 'umap_M1_lam0.2', by = 'clust_M1_lam0_k50_res5', type = 'discrete',
                main = 'Nonspatial clusters, \nres = 5', main.size = 40, pt.size = 6*0.15, col.discrete = mOD.cols, legend = FALSE),
  ncol = 3
)
dev.off()
