# supp_hypothalamus_lambda.R
# Lambda sweep on Moffitt et al. MERFISH mouse hypothalamus data
# Here, we use the mouse hypothalamus data to show what happens at large lambdas. 
# Load libraries

rm(list=ls())
graphics.off()

# set this to FALSE to perform clustering yourself. 
USE_PROVIDED_BANKSY_OBJ = TRUE
# If setting this to true, download the required banksy object from 
# https://www.dropbox.com/scl/fi/72skj0ms4akw1s3szfghh/bank_anim1_lamsweep.rds?rlkey=9rd7cubivtaonaf7nwc1ipnsl&dl=0

# params
animal_ID = 1 
k_geom = c(15, 30); 
lambda = seq(0, 1, 0.05)
npcs = 20; 
k_expr = 50; 
res = c(0.75)

library(Banksy)
library(gridExtra) # grid.arrange
library(ggplot2) # facet_wrap 
library(scales) # show_col
library(data.table) # fread
library(plyr) # mapvalues

out.dir = 'fig3-hypothalamus/out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'fig3-hypothalamus/out/merfish_supp_lambda'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)
data.dir = 'fig3-hypothalamus/data/'

if (USE_PROVIDED_BANKSY_OBJ){
  bank = readRDS( file = paste0(data.dir, '/bank_anim1_lamsweep.rds'))
} else {
  all_mfish = fread('fig3-hypothalamus/data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv')
  all_mfish <- all_mfish[,-c('Fos')]# remove Fos gene per Moffitt manuscript
  all_mfish = cbind(cell_ids = paste0('cell_', 1:nrow(all_mfish)), all_mfish)
  m = all_mfish[all_mfish$Animal_ID==animal_ID,]
  m = m[(!(m$Cell_class=='Ambiguous')),]
  expr = t(as.matrix(m[,-c(1:10)]))
  colnames(expr) <- m$cell_ids
  locs = as.data.frame(cbind(sdimx = m$Centroid_X, 
                             sdimy = m$Centroid_Y, 
                             sdimz = 1000*m$Bregma))
  rownames(locs) <- m$cell_ids
  bank <- BanksyObject(own.expr = expr,
                       cell.locs = locs
                       # )
                       ,
                       meta.data = m[, c('cell_ids',
                                         'Animal_ID',
                                         'Animal_sex',
                                         'Behavior',
                                         'Bregma',
                                         'Centroid_X',
                                         'Centroid_Y',
                                         'Cell_class',
                                         'Neuron_cluster_ID')])
  bank <- ComputeBanksy(bank, k_geom = k_geom)
  bank <- ScaleBanksy(bank)
  lambdas<-lambda
  #
  bank <- Banksy:::RunBanksyPCA(bank, lambda = lambdas, npcs = npcs)
  #
  print(names(bank@reduction))
  for (lam in lambdas){
    pdf(paste0(results.dir, '/nn-pca_-lam-',lam, '.pdf'), height = 6, width = 12)
    p1 <- plotReduction(bank, reduction = paste0('pca_M1_lam', lam))
    p2 <- plotScree(bank, lambda = lam)
    gridExtra::grid.arrange(p1, p2, ncol = 2)
    dev.off()
  }
  #
  bank <- Banksy:::RunBanksyUMAP(bank, lambda = lambdas, npcs = npcs, nneighbors = k_expr)
  
  clust_major = as.numeric(factor(unique(bank@meta.data$Cell_class)))
  names(clust_major) = as.character(factor(unique(bank@meta.data$Cell_class)))
  all(names(sort(clust_major))==levels(factor(unique(bank@meta.data$Cell_class))))
  
  clust_major_numeric = as.numeric(factor(bank@meta.data$Cell_class))
  clust_major_names = as.character(factor(bank@meta.data$Cell_class))
  clust_major_named_vec = clust_major_numeric
  names(clust_major_named_vec) = clust_major_names
  
  bank@meta.data = cbind(bank@meta.data, clust_major_num = clust_major_numeric)
  bank.backup = copy(bank)
  
  set.seed(42)
  ptm <- proc.time()
  bank <- ClusterBanksy(bank, lambda = lambdas, pca = TRUE, npcs = npcs,
                        method = 'leiden', k.neighbors = k_expr, resolution = res)
  clust.time = proc.time() - ptm
  print(clust.time)
  saveRDS(bank, file = paste0(data.dir, '/bank_anim1_lamsweep.rds'))
}

reorder_genes <- function(bank){
  x<-bank@own.expr
  d_gene <- dist(as.matrix(x))
  hc_gene <- hclust(d_gene)
  bank@own.expr <- bank@own.expr[hc_gene$order,]
  bank@nbr.expr <- bank@nbr.expr[hc_gene$order,]
  bank@harmonics$m1 = bank@harmonics$m1[hc_gene$order,]
  return(bank)
}
bank<- reorder_genes(bank)


##################### 
# Since we swept over all lambdas (from 0 to 1), we have both cell typing and domain segmentation 
# clustering runs. It does not make sense to do cluster consensus between domain segmentation
# and cell typing runs. Thus, to get a cluster color harmonized banksy object, we split the 
# object into cell typing runs, and domain segmentation runs, perform connectClusters, and then merge back. 
# strategy:
# map moffitt, lam 0 to 0.45, and 0.8 to lambda = 0.2. 
# map lam 0.5 to 1 to lambda = 0.8 

bank.1 = copy(bank)
bank.1@meta.data[clust.names(bank.1)[c(12:17, 19:22)]]<-NULL
bank.1 <- ConnectClusters(bank.1, map.to = 'clust_M1_lam0.2_k50_res0.75')

bank.2 = copy(bank)
bank.2@meta.data$clust_M1_lam0.8_k50_res0.75 = bank.1@meta.data$clust_M1_lam0.8_k50_res0.75
bank.2@meta.data[clust.names(bank.2)[c(1:11)]]<-NULL
bank.2 <- ConnectClusters(bank.2, map.to = 'clust_M1_lam0.8_k50_res0.75')

bank.1@meta.data$clust_M1_lam0.8_k50_res0.75<-NULL
bank.1@meta.data = cbind(bank.1@meta.data, bank.2@meta.data[clust.names(bank.2)])

bank.mapped = bank.1
bank.1 = NULL
bank.2 = NULL

num_clusters<-max(bank.mapped@meta.data[,clust.names(bank.mapped)])
hypo.cols<-Banksy:::getPalette(num_clusters)
names(hypo.cols)<-1:num_clusters

nested.plot.list = lapply(unique(bank.mapped@meta.data$Bregma)[c(7,9,11)], function(y){ #[c(1,3,5,7,9,11)]
  bank_layer = SubsetBanksy(bank.mapped,
                            metadata = Bregma %in% !!y)
  
  data.test1 <- Banksy:::getLocations(bank_layer)
  sdimx <- sdimy <- NULL
  
  inner.plot.list = lapply(clust.names(bank_layer)[round(seq(2, 22, 4))], function(x){
    feature1 <- Banksy:::getFeature(bank_layer,
                                    by = x, dataset = NULL)
    data.test1 <- cbind(data.test1, feature = as.factor(feature1))
    plot <- ggplot(data.test1, aes(x = sdimx, y = sdimy, col = as.factor(feature1))) +
      scale_color_manual(values = Banksy:::getDiscretePalette(feature1,
                                                              hypo.cols)) #+ 
    plot <- plot + geom_point(size = 1, alpha = 0.8) +theme_bw() 
    
    plot<-plot +theme_bw()+theme(legend.position="none") + 
      scale_x_discrete(labels = NULL,
                       breaks = NULL) + 
      labs(x = NULL)
    
    if (x == clust.names(bank_layer)[2]){
      plot = plot +ggtitle(paste0( y, 'mm') ) + 
        theme(plot.title = element_text(size = 36, face = "bold")) # , face = "bold"
    }
    # +
    # scale_y_discrete(labels = NULL, 
    #                  breaks = NULL) + labs(y = NULL)
    if (y == '-0.04'){
      clustering.run = paste0('lam = ' , gsub('clust_M1_lam', '', gsub('_k50_res0.75', '', x)))
      plot <- plot+scale_y_discrete(labels = NULL, 
                                    breaks = NULL) + labs(y = clustering.run) + 
        # theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) # ,axis.title=element_text(size=14,face="bold")
        theme(axis.title=element_text(size=36)) # ,axis.title=element_text(size=14,face="bold")
    } else {
      plot <- plot+scale_y_discrete(labels = NULL, 
                                    breaks = NULL) + labs(y = NULL)
    }
    plot
  })
})


unnested.plot.list = unlist(nested.plot.list, recursive = FALSE)



png(paste0(results.dir, '/','lamsweep_0to1', '.png'), 
    res = 50, units = 'in', height = length(unnested.plot.list)/length(nested.plot.list)*3, 
    width = length(nested.plot.list)*3)
plot.obj = gridExtra::grid.arrange(grobs = unnested.plot.list, 
                                   layout_matrix = matrix(data = 1:length(unnested.plot.list), 
                                                          byrow = FALSE, 
                                                          ncol = length(nested.plot.list))
)
print(plot.obj)
dev.off()




nested.plot.list = lapply(unique(bank.mapped@meta.data$Bregma)[c(7,9,11)], function(y){ #[c(1,3,5,7,9,11)]
  bank_layer = SubsetBanksy(bank.mapped,
                            metadata = Bregma %in% !!y)
  
  data.test1 <- Banksy:::getLocations(bank_layer)
  sdimx <- sdimy <- NULL
  
  inner.plot.list = lapply(clust.names(bank_layer)[round(c(17:22))], function(x){
    feature1 <- Banksy:::getFeature(bank_layer,
                                    by = x, dataset = NULL)
    data.test1 <- cbind(data.test1, feature = as.factor(feature1))
    plot <- ggplot(data.test1, aes(x = sdimx, y = sdimy, col = as.factor(feature1))) +
      scale_color_manual(values = Banksy:::getDiscretePalette(feature1,
                                                              hypo.cols)) #+ 
    plot <- plot + geom_point(size = 1, alpha = 0.8) +theme_bw() 
    
    plot<-plot +theme_bw()+theme(legend.position="none") + 
      scale_x_discrete(labels = NULL,
                       breaks = NULL) + 
      labs(x = NULL)
    
    if (x == clust.names(bank_layer)[2]){
      plot = plot +ggtitle(paste0( y, 'mm') ) + 
        theme(plot.title = element_text(size = 36, face = "bold")) # , face = "bold"
    }
    # +
    # scale_y_discrete(labels = NULL, 
    #                  breaks = NULL) + labs(y = NULL)
    if (y == '-0.04'){
      clustering.run = paste0('lam = ' , gsub('clust_M1_lam', '', gsub('_k50_res0.75', '', x)))
      plot <- plot+scale_y_discrete(labels = NULL, 
                                    breaks = NULL) + labs(y = clustering.run) + 
        # theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) # ,axis.title=element_text(size=14,face="bold")
        theme(axis.title=element_text(size=36)) # ,axis.title=element_text(size=14,face="bold")
    } else {
      plot <- plot+scale_y_discrete(labels = NULL, 
                                    breaks = NULL) + labs(y = NULL)
    }
    plot
  })
})


unnested.plot.list = unlist(nested.plot.list, recursive = FALSE)


png(paste0(results.dir, '/','lamsweep_075to1', '.png'), 
    res = 50, units = 'in', height = length(unnested.plot.list)/length(nested.plot.list)*3, 
    width = length(nested.plot.list)*3)
plot.obj = gridExtra::grid.arrange(grobs = unnested.plot.list, 
                                   layout_matrix = matrix(data = 1:length(unnested.plot.list), 
                                                          byrow = FALSE, 
                                                          ncol = length(nested.plot.list))
)
print(plot.obj)
dev.off()


