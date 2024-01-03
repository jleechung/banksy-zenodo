# supp_hypothalamus_oligo.R

rm(list=ls())
graphics.off()

################## # Load Data # ################## 
# Load the clustering result from the file process_hypothalamus_merfish.R
# In that script, users had the option for using our clustering solution 
# (using 'fig3-hypothalamus/data/banksyObj_provided.rds', or of performing the 
# clustering themselves (which takes ~24 hours on a single core). Here, we will 
# use our clustering solution for simplicity, but if users generated their own 
# clustering result, they may use that BANKSY object instead.

out.dir = 'fig3-hypothalamus/out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'fig3-hypothalamus/out/merfish_supp_rarecell'
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
hypo.cols<-Banksy:::getPalette(num_clusters)
names(hypo.cols)<-1:num_clusters
hypo.cols.old = hypo.cols
show_col(hypo.cols)


unlist(lapply(bank.conn@meta.data[clust.names(bank.conn)], function(x) length(unique(x))))
####

# plot the res = 0.5 run, and also plot the corresponding spatial plots for animal 1. 

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



idlist = list(c(1, 1), c(1,2), c(2,1), c(2,2))
clustering.run = c('clust_major_num', 'clust_M1_lam0.2_k50_res0.5' )
names(clustering.run) = c('Moffitt', 'BANKSY')
labels.to.plot = c(20, 26)
names(labels.to.plot) = c('Ependymal', 'Microglia')
hypo.cols.rare = rep('#E5E7E9', length(hypo.cols))
names(hypo.cols.rare) = names(hypo.cols)


bank.conn.raremeta = bank.conn@meta.data[clustering.run]
for (x in idlist){
  celltypeid = x[1] # 1 is microglia, 2 is ependymal
  clusteringrunid = x[2]
  curr.celltype.name = names(labels.to.plot)[celltypeid]
  bool.vect = as.numeric(bank.conn.raremeta[[clustering.run[clusteringrunid]]]==labels.to.plot[celltypeid])
  
  bank.conn.raremeta[[paste0(clustering.run[clusteringrunid], '_', curr.celltype.name)]] = bool.vect
}



lapply(idlist, function(x){
  celltypeid = x[1] # 1 is microglia, 2 is ependymal
  clusteringrunid = x[2]
  hypo.cols.rare[labels.to.plot[celltypeid]] = hypo.cols[labels.to.plot[celltypeid]]
  nested.plot.list = lapply(unique(bank.animal1@meta.data$Bregma)[c(8, 10, 12)], function(y){#[c(1,3,5,7,9,11)]
    bank_layer = SubsetBanksy(bank.animal1,
                              metadata = Bregma %in% !!y)
    # how to do this with a variable containing 'clust_M1_lam0_k50_res2.8'?
    
    # bank_layer.neur2 = SubsetBanksy(bank_layer,
    #                                 metadata = clust_banksy %in% !!labels.to.plot)
    data.test1 <- Banksy:::getLocations(bank_layer)
    sdimx <- sdimy <- NULL
    feature1 <- Banksy:::getFeature(bank_layer,
                                    by = clustering.run[clusteringrunid], dataset = NULL)
    
    data.test1 <- cbind(data.test1, feature = as.factor(feature1), 
                        tosort = as.numeric(feature1 %in% labels.to.plot[celltypeid]))
    sortorder.inh = order(data.test1$tosort)
    
    data.test1 <- data.test1[sortorder.inh,]
    plot <- ggplot(data.test1, aes(x = sdimx, y = sdimy, col = as.factor(feature1[sortorder.inh]))) +
      scale_color_manual(values = Banksy:::getDiscretePalette(feature1[sortorder.inh],
                                                              hypo.cols.rare)) #+ 
    plot <- plot + geom_point(size = 0.8, alpha = 0.8) +theme_bw()+ggtitle(paste0(y, 'mm'))
    
    plot<-plot +theme_bw()+theme(legend.position="none")+    scale_x_discrete(labels = NULL, 
                                                                              breaks = NULL) + labs(x = NULL) +
      scale_y_discrete(labels = NULL, 
                       breaks = NULL) + labs(y = NULL)
    plot
  })
  pdf(paste0(results.dir, '/','spatial_rare_3breg_', names(labels.to.plot)[celltypeid], 
             '_', names(clustering.run)[clusteringrunid] , '.pdf'),
      # res = 75, units = 'in',
      height = 2.5,
      width = length(nested.plot.list)*2.25)
  plot.obj = gridExtra::grid.arrange(grobs = nested.plot.list,  ncol = length(nested.plot.list))
  print(plot.obj)
  dev.off()
  
})


# compute ARIs// 
ependymal_ari = aricode::ARI(bank.conn.raremeta$clust_major_num_Ependymal, bank.conn.raremeta$clust_M1_lam0.2_k50_res0.5_Ependymal)
microglia_ari = aricode::ARI(bank.conn.raremeta$clust_major_num_Microglia, bank.conn.raremeta$clust_M1_lam0.2_k50_res0.5_Microglia)
print(ependymal_ari)
print(microglia_ari)
# > print(ependymal_ari)
# [1] 0.9578389
# > print(microglia_ari)
# [1] 0.8071589