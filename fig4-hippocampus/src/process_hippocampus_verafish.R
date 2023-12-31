# process_hippocampus_verafish.R
# You need the (legacy) version of BANKSY (v 0.1.5) for recreating this analysis: 
# remotes::install_github("prabhakarlab/Banksy@main") # to install 0.1.5. 
# See also the instructions in lines 133-151 below. 
# The data for this analysis is provided with the package. You can load it as follows:
library(Banksy)
library(ggplot2)
library(gridExtra)
library(scales)
library(ComplexHeatmap)
library(data.table)
library(irlba)

library(Matrix)

results.dir <- 'fig4-hippocampus/out/verafish'
data.dir <- 'fig4-hippocampus/data'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)
# 

k_geom = c(15, 30)
m.val = 1
lambda = 0.2
npcs = 20
k_expr = 50
res = 1.5
data(hippocampus)
expr <- hippocampus$expression
locs <- hippocampus$locations

# -------------------------------------------------------- #
# Load custom functions 
plotSpatialFeatures_ <- function (bank, by, type, dataset = NULL, nrow = NULL, ncol = NULL, 
                                  return.plot = FALSE, ...) 
{
  valid <- length(by) == length(type)
  if (!valid) 
    stop("Ensure a 1-1 correspondence between features to plot (by)\n                   and type of features (type)")
  plots <- Map(f = function(feature, feature.type, ...) {
    plotSpatial_(bank, dataset = dataset, by = feature, type = feature.type, 
                 ...)
  }, by, type, ...)
  do.call(grid.arrange, c(plots, nrow = nrow, ncol = ncol))
  if (return.plot) 
    return(plots)
}
plotSpatial_ <- function (bank, dataset = NULL, by = NA, type = c("discrete", 
                                                                  "continuous"), pt.size = 0.5, 
                          pt.alpha = 0.7, col.midpoint = NULL, 
                          col.lowpoint = NULL, col.highpoint = NULL, col.low = "blue", 
                          col.mid = "gray95", col.high = "red", na.value = "gray", 
                          col.discrete = NULL, main = NULL, main.size = 5, legend = TRUE, 
                          legend.text.size = 6, legend.pt.size = 3, wrap = FALSE) 
{
  data <- Banksy:::getLocations(bank, dataset = dataset)
  sdimx <- sdimy <- NULL
  if (is.na(by)) 
    plot <- ggplot(data, aes(x = sdimx, y = sdimy))
  else {
    feature <- Banksy:::getFeature(bank, by = by, dataset = dataset)
    Banksy:::checkType(type)
    if (type == "continuous") {
      if (is.null(col.midpoint)) 
        col.midpoint <- median(feature)
      if ((is.null(col.highpoint)) | (is.null(col.lowpoint))) {
        plot <- ggplot(data, aes(x = sdimx, y = sdimy, 
                                 col = feature)) + scale_color_gradient2(midpoint = col.midpoint, 
                                                                         low = col.low, mid = col.mid, high = col.high)
      }
      else {
        feature[feature > col.highpoint] <- col.highpoint
        feature[feature < col.lowpoint] <- col.lowpoint
        plot <- ggplot(data, aes(x = sdimx, y = sdimy, 
                                 col = feature)) + scale_color_gradient2(midpoint = col.midpoint, 
                                                                         limits = c(col.lowpoint, col.highpoint), low = col.low, 
                                                                         mid = col.mid, high = col.high, na.value = na.value)
      }
    }
    if (type == "discrete") {
      data <- cbind(data, feature = as.factor(feature))
      plot <- ggplot(data, aes(x = sdimx, y = sdimy, col = as.factor(feature))) + 
        scale_color_manual(values = Banksy:::getDiscretePalette(feature, 
                                                                col.discrete)) + 
        guides(color = guide_legend(override.aes = list(size = legend.pt.size)))
    }
  }
  legend.pos <- ifelse(legend, "right", "none")
  plot <- plot + geom_point(size = pt.size, alpha = pt.alpha) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    ggtitle(main) + xlab("") + ylab("") + 
    Banksy:::theme_blank(legend.text.size = legend.text.size, main.size = main.size, 
                         legend.pos = legend.pos)
  wrap.z <- "sdimz" %in% names(data)
  if (wrap.z & wrap) 
    plot <- plot + facet_grid(sdimz ~ feature)
  else if (wrap) 
    plot <- plot + facet_wrap(~feature)
  else if (wrap.z) 
    plot <- plot + facet_wrap(~sdimz)
  return(plot)
}
reorder_genes <- function(bank, normalize_cols = FALSE, scale_rows = FALSE){
  # the data in bank must be unnormalized and unscaled. 
  gene_names = rownames(bank@own.expr)
  x<-bank@own.expr
  if (normalize_cols){
    x<- t(t(x)/colSums(x))*100
  }
  if (scale_rows){
    x<-t(scale(t(x), center = TRUE, scale = TRUE))
  }
  
  d_gene <- dist(as.matrix(x))
  hc_gene <- hclust(d_gene)
  bank@own.expr <- bank@own.expr[hc_gene$order,]    
  bank@nbr.expr <- bank@nbr.expr[hc_gene$order,]  
  return(bank)
}

# -------------------------------------------------------- #
# /end functions


total_count <- colSums(expr)
num_genes <- colSums(expr > 0)
meta <- data.frame(total_count = total_count, num_genes = num_genes)

# You may load banksy object from our run, or uncomment the code below 
# (start and end Banksy clustering) and generate your own. To download our BANKSY object, 
# go, use the following link and place it in the 'fig4-hippocampus/data' directory. 
# https://www.dropbox.com/scl/fi/7wwb3oge18f4dsd09t269/bank.rds?rlkey=kln41aetg3louna430r7l50ue&dl=0
# (If you downloaded the code from Zenodo, then the banksy object is already in the directory; 
# you only need download it if you are using the code from the GitHub repo.)
bank = readRDS(file = paste0(data.dir, '/bank.rds'))

# If you do run the code yourself, see the notes below in lines 432-438 and on lines 
# 500-508 on correctly identifying the correct pairs of clusters to generate DE genes.
# In brief, you must identify which clusters correspond to the fornix vs thalamic oligos, 
# and the CA3 vs SSC neurons. For the oligos, you can do so by looking at the heatmap 
# generated below, and looking for marker genes which we specify below (since 
# the markers mark both clusters of interest). For the neurons, the markes mark both 
# the neuronal clusters of interest, and other neuronal clusters too (with subtle 
# differences in expression across different neuronal clusters). Thus, we take the alternate 
# approach of plotting the spatial distributions of the clusters, and look at the 
# paper to see what the distribution of the  pair of clusters 
# looks like to identify their cluster ID numbers. 
# ---------------- start banksy clustering -----------------
# bank <- BanksyObject(own.expr = expr, cell.locs = locs, meta.data = meta)
# bank <- SubsetBanksy(bank, metadata = total_count > quantile(total_count, 0.05) &
#                        total_count < quantile(total_count, 0.98))
# bank <- NormalizeBanksy(bank)
# bank <- ComputeBanksy(bank, k_geom = k_geom, M = m.val)
# bank <- ScaleBanksy(bank)
# lam <- c(0, lambda)
# bank <- Banksy:::RunBanksyPCA(bank, lambda = lam, npcs = npcs, M = m.val)
# bank <- Banksy:::RunBanksyUMAP(bank, lambda = lam, npcs = npcs, nneighbors = k_expr, M = m.val)
# set.seed(42)
# bank <- ClusterBanksy(bank, lambda = lam, pca = TRUE, npcs = npcs, M = m.val,
#                       method = 'leiden', k.neighbors = k_expr, resolution = res)
# ---------------- end banksy clustering -----------------





# You can uncomment the lines below and run MERINGUE. Or you can use the provided results (lines 226-233)
# library(MERINGUE)
# filterDists<-c(750)
# w2 <- getSpatialNeighbors(bank@cell.locs, filterDist = filterDists[1], verbose=TRUE)
# 
# pca <- prcomp_irlba(t(as.matrix(bank@own.expr)), n = 20)
# #
# set.seed(42)
# ptm <- proc.time()
# com5_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 5, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com5_42_w2, file = paste0(results.dir, '/com5_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com6_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 6, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com6_42_w2, file = paste0(results.dir, '/com6_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com7_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 7, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com7_42_w2, file = paste0(results.dir, '/com7_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com8_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 8, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com8_42_w2, file = paste0(results.dir, '/com8_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com9_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 9, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com9_42_w2, file = paste0(results.dir, '/com9_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com10_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 10, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com10_42_w2, file = paste0(results.dir, '/com10_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com12_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 12, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com12_42_w2, file = paste0(results.dir, '/com12_w2_seed42', '.rds'))
# 
# set.seed(42)
# ptm <- proc.time()
# com15_42_w2 <- getSpatiallyInformedClusters(pca$x[,1:20], W = w2, k = 15, verbose = TRUE)
# print(proc.time() - ptm)
# saveRDS(com15_42_w2, file = paste0(results.dir, '/com15_w2_seed42', '.rds'))
# 

com5_w2 = readRDS(paste0(data.dir, '/com5_w2_seed42', '.rds'))
com6_w2 = readRDS(paste0(data.dir, '/com6_w2_seed42', '.rds'))
com7_w2 = readRDS(paste0(data.dir, '/com7_w2_seed42', '.rds'))
com8_w2 = readRDS(paste0(data.dir, '/com8_w2_seed42', '.rds'))
com9_w2 = readRDS(paste0(data.dir, '/com9_w2_seed42', '.rds'))
com10_w2 = readRDS(paste0(data.dir, '/com10_w2_seed42', '.rds'))
com12_w2 = readRDS(paste0(data.dir, '/com12_w2_seed42', '.rds'))
com15_w2 = readRDS(paste0(data.dir, '/com15_w2_seed42', '.rds'))

bank@meta.data<-cbind(bank@meta.data,
                      clust_dist750_k5 = as.numeric(com5_w2), 
                      clust_dist750_k6 = as.numeric(com6_w2), 
                      clust_dist750_k7 = as.numeric(com7_w2), 
                      clust_dist750_k8 = as.numeric(com8_w2), 
                      clust_dist750_k9 = as.numeric(com9_w2), 
                      clust_dist750_k10 = as.numeric(com10_w2), 
                      clust_dist750_k12 = as.numeric(com12_w2), 
                      clust_dist750_k15 = as.numeric(com15_w2))

# ## load HMRF and other methods' results.
# alt_methods_bank <- readRDS('Jun1_bank_vfish_hip_0728_multi_domain.rds')
bayes.giotto.fict = readRDS(file = paste0(data.dir, '/hippocampus_0728_bayespace_giotto_fict.rds'))
bayes.giotto.fict.num = lapply(bayes.giotto.fict, as.numeric)

spagcn = read.csv(file = paste0(data.dir, '/spaGCN_clusters_0728.csv'))
names(spagcn) = gsub('X0728_pred_', 'clust_spagcn_', names(spagcn))
spagcn.num = lapply(spagcn[,-c(1:3)], as.numeric)
spagcn.num = lapply(spagcn.num, FUN = function(x) x+1)

# Mar 30, 2023, Now update with the new version of spicemix 
new.spicemix = readRDS(file = paste0(data.dir, 
                                     '/bank_vfish_hip_0728_multi_domain_revision.rds'))
spicedf = new.spicemix@meta.data[,-c(1:13)]
# all(old.giotto.fict.spice@meta.data$cell_ID == new.spicemix@meta.data$cell_ID)

all(spagcn$X==bank@meta.data$cell_ID)
bank@meta.data = cbind(bank@meta.data, spagcn.num, bayes.giotto.fict.num, spicedf)

# # check cell ordering before appending results across objects.
# stopifnot(all(bank@meta.data$cell_ID==alt_methods_bank@meta.data$cell_ID))

bank <- ConnectClusters(bank, map.to = 'clust_M1_lam0.2_k50_res1.5')
cnms <- clust.names(bank)

num_clusters<-max(bank@meta.data[,cnms])
hippo.cols.initial<-Banksy:::getPalette(num_clusters)
names(hippo.cols.initial)<-1:num_clusters
# 
cluster.colors = c('14' = 'darkgray',
                   '13' = 'bisque4',
                   '3' = 'goldenrod1',
                   '11' = 'chocolate4',
                   '15' = 'darkorange',
                   '12' = 'deeppink2',
                   '7' = 'darkolivegreen4',
                   '10' = 'darkblue',
                   '16' = 'aquamarine2',
                   '6' = 'cornflowerblue')
cluster.colors.names <- names(cluster.colors)
cluster.colors <- gplots::col2hex(cluster.colors)
names(cluster.colors) <- cluster.colors.names
show_col(cluster.colors)
hippo.cols <- hippo.cols.initial
hippo.cols[cluster.colors.names] <- cluster.colors
# hippo.cols = c(hippo.cols, '17' = "#8DB600")
show_col(hippo.cols)


all.runs <- clust.names(bank)

cluster.colors = c('14' = 'coral3', 
                   '13' = '#8C7E6C',
                   '3' = 'goldenrod1', 
                   '11' = 'chocolate4', 
                   '15' = 'darkorange',
                   '12' = '#ED1B89',
                   '7' = '#8E55A2', 
                   '10' = '#008856', 
                   '16' = 'aquamarine', 
                   '6' = 'cornflowerblue', 
                   '17' = '#8DB73F', 
                   '8' = '#243E90', 
                   '2' = 'darkolivegreen4' 
)
cluster.colors.names <- names(cluster.colors)
cluster.colors <- gplots::col2hex(cluster.colors)
names(cluster.colors) <- cluster.colors.names
show_col(cluster.colors)
hippo.cols.mid <- hippo.cols
hippo.cols[cluster.colors.names] <- cluster.colors

par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
show_col(hippo.cols, ncol = 4)
show_col(hippo.cols.mid, ncol = 4)

for (i in 1:length(all.runs)){
  subplot.nrow <- ceiling(length(unique(bank@meta.data[[all.runs[i]]])))
  print(subplot.nrow)
  png(paste0(results.dir, '/spatial_tall_smaller_',all.runs[i], '.png'),
      height = subplot.nrow*1.5, width = 2.5, res = 80, units = 'in')
  pp<-plotSpatial(bank, 
                  by = all.runs[i], 
                  main = gsub('clust_', '', all.runs[i]),
                  main.size = 12,
                  col.discrete = hippo.cols,
                  pt.size = .4, 
                  legend = FALSE,
                  type = 'discrete') + facet_wrap(~feature, ncol = 1)
  print(pp)
  dev.off()
}

clust.names(bank)
all.runs <- clust.names(bank)
gene_names = rownames(bank@own.expr)
x<-bank@own.expr

d_gene <- dist(as.matrix(x))
hc_gene <- hclust(d_gene)

print(gene_names[hc_gene$order][1:10])
genes_ord<-gene_names[hc_gene$order]
print(rownames(bank@own.expr)[1:10])
print(rownames(bank@nbr.expr)[1:10])
print(rownames(bank@harmonics$k1)[1:10])

bank@own.expr <- bank@own.expr[hc_gene$order,]    
bank@nbr.expr <- bank@nbr.expr[hc_gene$order,] 
bank@harmonics$k1 <- bank@harmonics$k1[hc_gene$order,]  


print(rownames(bank@own.expr)[1:10])
print(rownames(bank@nbr.expr)[1:10])
print(rownames(bank@harmonics$m1)[1:10])


## now extract out the key runs and regenerate the figures
nonspatial.annot <- 'clust_M1_lam0_k50_res1.5'
banksy.annot <- 'clust_M1_lam0.2_k50_res1.5'
meringue.annot <- "clust_dist750_k10" 
hmrf.annot <- "clust_giotto18_beta_12"
bspace.annot <- 'clust_bspace18'
spicemix.annot <- 'clust_spicemix_N_18'
# spagcn.annot0p5 <- 'clust_spagcn_0.5_18'
spagcn.annot1 <- 'clust_spagcn_1_18'
# fict.annot <- 'clust_fict18'
fict.annot <- 'clust_fict20' 

multiple.annots<-c(nonspatial = nonspatial.annot,
                   banksy = banksy.annot,
                   meringue = meringue.annot,
                   hmrf = hmrf.annot,
                   bspace = bspace.annot, 
                   spicemix = spicemix.annot, 
                   spagcn1 = spagcn.annot1,
                   fict = fict.annot)
# 


pdf(paste0(results.dir, '/','hm_', 'all_main_methods' , '.pdf'), height = 20/2, width = 32/2)
p.all.mk<-Banksy::plotHeatmap(bank,
                              assay = 'own.expr',
                              # features = rownames(bank2@own.expr),
                              # M = m.val,
                              # lambda = 0.2,
                              cex.row = 5,
                              annotate = TRUE,
                              annotate.by = multiple.annots,
                              annotation.size = 10,
                              order.by = banksy.annot,
                              col.discrete = hippo.cols,
                              annotation.name = TRUE,
                              rasterize = TRUE,
                              barplot.by = c('NODG', 'nCount'),
                              cluster_row_slices = FALSE,
                              cluster_rows = FALSE,
                              cluster_column_slices = TRUE,
                              cluster_columns = TRUE)
print(p.all.mk)
dev.off()

# de gene analysis
mk_banksy_w_pair <- scran::pairwiseWilcox(bank@own.expr, 
                                          groups = bank@meta.data$clust_M1_lam0.2_k50_res1.5)
mk_banksy_t_pair <- scran::pairwiseTTests(bank@own.expr, 
                                          groups = bank@meta.data$clust_M1_lam0.2_k50_res1.5)
topmk_w<-scran::getTopMarkers(mk_banksy_w_pair$statistics, 
                              pairs = mk_banksy_w_pair$pairs)
topmk_t<-scran::getTopMarkers(mk_banksy_t_pair$statistics, 
                              pairs = mk_banksy_t_pair$pairs)


topmk_w_combined_any<-scran::getTopMarkers(mk_banksy_w_pair$statistics, 
                                           pairs = mk_banksy_w_pair$pairs, 
                                           pairwise = FALSE)

topmk_t_combined_any<-scran::getTopMarkers(mk_banksy_t_pair$statistics, 
                                           pairs = mk_banksy_t_pair$pairs, 
                                           pairwise = FALSE)

topmk_w_combined_all<-scran::getTopMarkers(mk_banksy_w_pair$statistics, 
                                           pairs = mk_banksy_w_pair$pairs, 
                                           pairwise = FALSE, pval.type = 'all' )

saveRDS(topmk_w_combined_all, file = paste0(results.dir, '/topmk_w_combined_all.rds'))

# If you did not use our banksy object and clustering result, and instead ran the code in lines 140-153
# yourself, you might need to manually identify which two clusters numbers correspond to the 
# two oligodendrocyte clusters, and then populate lines 440-441 below with those two numbers 
# (in the provided banksy clustering, these are clusters 7 and 10). The best way to do this 
# is by looking at the heatmap just generated, which can be found in fig4-hippocampus/out/verafish/hm_all_main_methods.pdf
# The second annotation row in this heatmap corresponds to the BANKSY clustering (lam0.2). The 
# pair of clusters you are looking for are marked by high expression of Mbp, Enpp2, Plp1 and Sgk1.

cl1.oligo = as.character(7)
cl2.oligo = as.character(10)

set2w1 = topmk_w[[cl1.oligo]][[cl2.oligo]]

set2t1 = topmk_t[[cl1.oligo]][[cl2.oligo]]

fimbria_de <- unique(c(set2w1, set2t1)) 

bank.cells<-bank@meta.data$cell_ID[ bank@meta.data[[banksy.annot]] %in% c(cl1.oligo, cl2.oligo)] 
bank.subset<-SubsetBanksy(bank, cells = bank.cells, features = fimbria_de)

gene_names = rownames(bank.subset@own.expr)
x<-bank.subset@own.expr
d_gene <- dist(as.matrix(x))
hc_gene <- hclust(d_gene)
genes_ord<-gene_names[hc_gene$order]
genes_ord = c('Mobp',
              'Plp1',
              'Gfap',
              'Bcas1',
              'Mapk9',
              'Sparcl1',
              'Atp1a2',
              'Atp1b2',
              'Nefl',
              'Map',
              'Clstn2')
union(setdiff(fimbria_de, genes_ord) , setdiff(genes_ord, fimbria_de))

saveRDS(genes_ord, file = paste0(results.dir, '/fimbria_de_Oct25.rds'))


bank.subset@own.expr <- bank.subset@own.expr[order(match(rownames(bank.subset@own.expr),genes_ord)),]    
bank.subset@nbr.expr <- bank.subset@nbr.expr[order(match(rownames(bank.subset@own.expr),genes_ord)),] 
bank.subset@harmonics$k1 <- bank.subset@harmonics$k1[order(match(rownames(bank.subset@own.expr),genes_ord)),]  

bank.subset <- ScaleBanksy(bank.subset)
png(paste0(results.dir, '/','hm_', 'fimbria_scaled' , '.png'), 
    height = 5, width = 12, res = 200,
    units = 'in')
plotHeatmap(bank.subset, 
            # cells = cells,
            # features = fimbria_de, 
            assay = 'own.expr', 
            cex.row = 12,
            annotate = TRUE,
            annotate.by = multiple.annots[c(1,2)],
            annotation.size = 12,
            order.by = banksy.annot,
            col.discrete = hippo.cols,
            annotation.name = TRUE,
            rasterize = TRUE,
            barplot.by = c('NODG', 'nCount'),
            cluster_row_slices = FALSE,
            cluster_rows = FALSE,
            cluster_column_slices = TRUE,
            cluster_columns = TRUE)
dev.off()

# If you ran the banksy clustering yourself (instead of using our provided banksy object),
# you will need to identify the cluster number IDs corresponding to the CA3 and somatosensory
# cortex neurons and populate lines 511-512. In the oligodendrocytes case above, we were able to do this using 
# marker genes. But here, there are several neuronal clusters, with only subtle differences 
# in marker gene expression (the markers being Rtn1, Tbr1, Parp1, Clstn2, Cadm3, Cacna1a, among many others),
# so to identify the two neuronal clusters of interest, we must follow look at the spatial distributions
# of the clusters. Open the file out/spatial_tall_smaller_clust_M1_lam0.2_k50_res1.5.png generated above. 
# Here, you will find the clusters corresponding to the CA3 neurons and the SSC neurons (see main Fig 4 for 
# the shape / location of these clusters.). In our banksy object, their numerical IDs are 14, 18. 

# For identifying the cluster number of the two neuronal clusters 
cl1.neuron = as.character(14)
cl2.neuron = as.character(18)

set3w1 = topmk_w[[cl1.neuron]][[cl2.neuron]]
set3t1 = topmk_t[[cl1.neuron]][[cl2.neuron]]


ca3_de <- unique(c(set3w1, set3t1))
bank.cells<-bank@meta.data$cell_ID[ bank@meta.data[[banksy.annot]] %in% c(cl1.neuron,cl2.neuron)]
bank.subset<-SubsetBanksy(bank, cells = bank.cells, features = ca3_de)
genes_ord = c('Sparcl1' ,'Cd34' , 'Grm3'   ,'Egr1'  ,  'Tbr1' , 
              'Nefl',
              'Cacna1a', 'Cadm3' ,    'Clstn2', 
              'Gnaq'  )

saveRDS(genes_ord, file = paste0(results.dir, '/CA3_de_Oct25.rds'))

bank.subset@own.expr <- bank.subset@own.expr[order(match(rownames(bank.subset@own.expr),genes_ord)),]    
bank.subset@nbr.expr <- bank.subset@nbr.expr[order(match(rownames(bank.subset@own.expr),genes_ord)),] 
bank.subset@harmonics$k1 <- bank.subset@harmonics$k1[order(match(rownames(bank.subset@own.expr),genes_ord)),]  


bank.subset <- ScaleBanksy(bank.subset)
png(paste0(results.dir, '/','hm_', 'CA3' , '.png'), height = 6, width = 12, res = 200, units = 'in')
plotHeatmap(bank.subset, 
            assay = 'own.expr', 
            cex.row = 12,
            annotate = TRUE,
            annotate.by = multiple.annots[c(1,2)],
            annotation.size = 12,
            order.by = banksy.annot,
            col.discrete = hippo.cols,
            annotation.name = TRUE,
            rasterize = TRUE,
            barplot.by = c('NODG', 'nCount'),
            cluster_row_slices = FALSE,
            cluster_rows = FALSE,
            cluster_column_slices = TRUE,
            cluster_columns = TRUE)
dev.off()



all.mk1 <- Reduce(union, lapply(FUN = head, X = topmk_w_combined_all, 
                                n = 5))


all.mk3 <- Reduce(union, lapply(FUN = head, X = topmk_w_combined_all, 
                                n = 3))

all.mk2 = Reduce(union, list(fimbria_de, ca3_de))
all.mk = union(all.mk1, all.mk2)
bank.subset.all.mk<-SubsetBanksy(bank, features = all.mk)
d_gene.all.mk <- dist(bank.subset.all.mk@own.expr)
hc_gene.all.mk <- hclust(d_gene.all.mk)
bank.subset.all.mk@own.expr <- bank.subset.all.mk@own.expr[hc_gene.all.mk$order,] 
bank.subset.all.mk@nbr.expr <- bank.subset.all.mk@nbr.expr[hc_gene.all.mk$order,] 

# ---------- spatial and umap plots ----------- #
all.runs <- (multiple.annots) #names(bank@meta.data)[4:9]
reductions <- rep('umap_M1_lam0', 8) # march 2023 changed from 9 because we removed spagcn 0.5

png(paste0(results.dir, '/umap_0_all.png'), height = 8, 
    width = 16, res = 100, units = 'in')#, units = 'in',  res = 200)
umapdims <- mapply(FUN = function(reduction, by) plotReduction(bank, 
                                                               reduction = reduction, 
                                                               by = by, 
                                                               main = by,
                                                               col.discrete = hippo.cols,
                                                               type = 'discrete', 
                                                               pt.size = 0.2, main.size = 10), 
                   reduction = as.list(reductions), 
                   by = (all.runs), 
                   SIMPLIFY = FALSE)
do.call("grid.arrange", c(umapdims, ncol = 4))
dev.off()

reductions <- rep('umap_M1_lam0.2', 8)
png(paste0(results.dir, '/umap_0.2_all.png'), height = 8, width = 16, res = 100, units = 'in')#, units = 'in',  res = 200)
umapdims <- mapply(FUN = function(reduction, by) plotReduction(bank, 
                                                               reduction = reduction, 
                                                               by = by, 
                                                               main = by,
                                                               col.discrete = hippo.cols,
                                                               type = 'discrete', 
                                                               pt.size = 0.2, main.size = 10), 
                   reduction = as.list(reductions), 
                   by = all.runs, 
                   SIMPLIFY = FALSE)
do.call("grid.arrange", c(umapdims, ncol = 4))
dev.off()

#### spatial
# all.runs <- c(bank.runs, meringue.runs)
# for (i in 1:length(all.runs)){
#   subplot.nrow <- ceiling(length(unique(bank@meta.data[[all.runs[i]]]))/4)
#   print(subplot.nrow)
#   pdf(paste0(results.dir, '/spatial_final_',all.runs[i], '.pdf'),
#       height = subplot.nrow*3, width = 15)
#   pp<-plotSpatial(bank, 
#                   by = all.runs[i], 
#                   col.discrete = hippo.cols,
#                   pt.size = 1, 
#                   type = 'discrete') + facet_wrap(~feature, ncol = 4)
#   print(pp)
#   dev.off()
# }

for (i in 1:length(all.runs)){
  subplot.nrow <- ceiling(length(unique(bank@meta.data[[all.runs[i]]]))/4)
  print(subplot.nrow)
  png(paste0(results.dir, '/spatial_final_',all.runs[i], '.png'),
      height = subplot.nrow*3, width = 15, res=300, units = 'in')
  pp<-plotSpatial(bank, 
                  by = all.runs[i], 
                  col.discrete = hippo.cols,
                  pt.size = 1, 
                  type = 'discrete') + facet_wrap(~feature, ncol = 4)
  print(pp)
  dev.off()
}


# ---------- marker gene plots ----------#

mk_comb<-function(x, N, topmk){
  return(topmk[[as.character(x[1])]][[as.character(x[2])]][1:N])
}

png(paste0(results.dir, '/all_mk_rescaled', '.png'), height = 35, 
    width = 21, res = 300, units = 'in')
plotSpatialFeatures_(bank, by = all.mk, main = all.mk, main.size = 20,
                        col.highpoint=apply(bank@own.expr[all.mk,], 1, quantile, probs = 0.995),
                        col.lowpoint=apply(bank@own.expr[all.mk,], 1, quantile, probs = 0),
                        col.midpoint=apply(bank@own.expr[all.mk,], 1, quantile, probs = 0.5),
                        type = rep('continuous', length(all.mk)), 
                        nrow = 13, ncol = 5, pt.size = 0.1, legend.text.size = 12)
dev.off()

already.plotted.genes = c('Cadm3', 'Clstn2', 'Grm3', 'Egr1', 'Mobp', 'Bcas1', 'Sparcl1', 'Nefl')
fimbria_to_plot <- fimbria_de[!(fimbria_de %in% already.plotted.genes)]
ca3_to_plot <- ca3_de[!(ca3_de %in% union(already.plotted.genes, fimbria_to_plot))]

fimbria_and_ca3_clusters = c(cl1.oligo, cl2.oligo, cl1.neuron, cl2.neuron) 

already.plotted2 <- unique(c(already.plotted.genes, fimbria_to_plot, ca3_to_plot))
all.mk3 <- Reduce(union, lapply(FUN = head, 
                                X = (topmk_w_combined_all)[setdiff(names(topmk_w_combined_all), 
                                                                   fimbria_and_ca3_clusters)], 
                                n = 3))


mk.to.plot <- all.mk3[!(all.mk3 %in% already.plotted2)]

# all three sets together
all.three.mk.sets <- c( ca3_to_plot, fimbria_to_plot, mk.to.plot)
png(paste0(results.dir, '/remaining_mk_rescaled', '.png'), 
    height = ceiling(length(all.three.mk.sets)/6)*3, width = 20, res = 200, units = 'in')
plotSpatialFeatures_(bank, by = all.three.mk.sets, main = all.three.mk.sets, 
                        main.size = 20,
                        col.highpoint=apply(bank@own.expr[all.three.mk.sets,], 
                                            1, quantile, probs = 0.99),
                        col.lowpoint=apply(bank@own.expr[all.three.mk.sets,], 
                                           1, quantile, probs = 0),
                        col.midpoint=apply(bank@own.expr[all.three.mk.sets,], 
                                           1, quantile, probs = 0.5),
                        type = rep('continuous', 
                                   length(all.three.mk.sets)), 
                        nrow = ceiling(length(all.three.mk.sets)/6), ncol = 6, 
                        pt.size = 0.1, legend= FALSE)
dev.off()

# 