# process_hippocampus_scrna.R
# # In this file, we perform the verafish result validation using scrna seq data. 

cutoff = 0.01
library(Seurat)
library(Banksy)
library(data.table)
library(scater)
library(scran)

data(hippocampus)
expr <- hippocampus$expression
locs <- hippocampus$locations
num_top_genes_to_use = 25
results.dir <- paste0('fig4-hippocampus/out/scrna_neurons')

check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)

# download from: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq
# the following line in the table: 
# Gene expression matrix (Seurat)	5.4 GB	.ss.rda	This Seurat object contains the cell-by-gene expression matrix, with introns and exons combined.
# specific link: https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/Seurat.ss.rda

# -------------------------------- NOTE -------------------------------- #
# Uncomment the lines below between "START PREPROCESSING BLOCK" and "END PREPROCESSING BLOCK" to generate the 
# individual seurat objects (ss.seurat.SSC.CA3_unprocessed for neurons, ss.seurat.oligo_unprocessed for 
# oligodendrocytes.)
# You will need a computer with 32 GB RAM to run this block, which is why we provide the subsetted seurat objects, 
# which require much less RAM to work with. 
# -------------------------------- ---- -------------------------------- #

# START PREPROCESSING BLOCK (peak usage 17-20 GB RAM , so you need a 32 GB machine, probably)
# load(file = 'fig4-hippocampus/data/Seurat.ss.rda') # 16.12 GB
# # # subset out the CA3 and SSc cells. 
# ss.seurat.SSC.CA3_unprocessed = subset(x = ss.seurat,
#                                 subset = region_label%in% c("SSp", "HIP")&subclass_label %in% c("L5/6 NP CTX",
#                                                                                                 "L6 CT CTX",
#                                                                                                 "L6 IT CTX",
#                                                                                                 "L6b CTX",
#                                                                                                 "CA3"))
# saveRDS(ss.seurat.SSC.CA3_unprocessed, file = 'fig4-hippocampus/data/ss_seurat_SSC_CA3_unprocessed.rds')
# ss.seurat.oligo_unprocessed = subset(x = ss.seurat, subset = subclass_label=="Oligo")
# saveRDS(ss.seurat.oligo_unprocessed, file = 'fig4-hippocampus/data/ss_seurat_oligo_unprocessed.rds')
# ss.seurat <- NULL
# END PREPROCESSING BLOCK

ss.seurat.SSC.CA3_unprocessed<-readRDS(file = 'fig4-hippocampus/data/ss_seurat_SSC_CA3_unprocessed.rds')
data.counts <- GetAssayData(ss.seurat.SSC.CA3_unprocessed, slot = 'counts')

ca3_de <- readRDS('fig4-hippocampus/data/CA3_de.rds')
fimbria_de <- readRDS('fig4-hippocampus/data/fimbria_de.rds')

ss.seurat.SSC.CA3.ln <- NormalizeData(object = ss.seurat.SSC.CA3_unprocessed, 
                                      normalization.method = 'LogNormalize')
ss.seurat.SSC.CA3.ln.scaled <- ScaleData(object = ss.seurat.SSC.CA3.ln)

data.ln.scaled <- GetAssayData(ss.seurat.SSC.CA3.ln.scaled, slot = 'scale.data')
ca3_de<-ca3_de[ca3_de %in% rownames(data.ln.scaled)]

find_corr_genes <- function(gene, data.mat) {
  message(gene)
  cor <- rep(NA, nrow(data.mat))
  names(cor) <- rownames(data.mat)
  x <- data.mat[gene,]
  for (i in 1:nrow(data.mat)) {
    if (i %% 1000 == 0) message('Iter. ', i)
    if (grepl('^Blank', rownames(data.mat)[i])) next
    cor[i] <- cor(x, data.mat[i,])
  }
  cor <- sort(cor, decreasing = TRUE)
  return(cor)
}

cor.genes <- lapply(ca3_de, find_corr_genes, data.mat = data.ln.scaled)
top.genes.list <- lapply(cor.genes, function(x) x[1:num_top_genes_to_use])

remove.lowly.expressed.genes <- function(x){
  Ax = as.matrix(data.counts[names(x),])
  Bx = Ax>0
  cx = rowSums(Bx)
  x.to.keep = names(which(cx>cutoff*ncol(Ax)))
}
top.genes.list.filtered = lapply(top.genes.list, remove.lowly.expressed.genes) 
# no genes removed, as seen here:
mapply(function(x, y) all(names(x)==y),top.genes.list, top.genes.list.filtered)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

top.genes <- unique(unlist(top.genes.list.filtered)) 
seu.topgenes <- subset(x = ss.seurat.SSC.CA3_unprocessed, features = top.genes)
cluster.and.plot <- function(seu.topgenes, top.genes, res = 0.2,
                             normalization=NULL, # 'RC', 'LogNormalize'
                             ){
  npcs = 5
  if (!(is.null(normalization))){
    seu.topgenes <- NormalizeData(object = seu.topgenes, 
                                  normalization.method = normalization)
  }
  seu.topgenes <- ScaleData(object = seu.topgenes)
  seu.topgenes <- Seurat:::RunPCA(object = seu.topgenes, 
                                  features = seu.topgenes@assays$RNA@counts@Dimnames[[1]])
  seu.topgenes <- FindNeighbors(object = seu.topgenes, dims = 1:npcs)
  seu.topgenes <- FindClusters(object = seu.topgenes, resolution = res)
  seu.topgenes <- RunTSNE(object = seu.topgenes, dims = 1:npcs)
  seu.topgenes <- Seurat:::RunUMAP(object = seu.topgenes, dims = 1:npcs)
  top.five.genes <- lapply(top.genes, function(x) x[1:5])
  
  return(list(seu.topgenes = seu.topgenes, top.five.genes = top.five.genes))
}
ln.scaled.list.clus = cluster.and.plot(seu.topgenes, top.genes.list, 
                                       normalization='LogNormalize') 

seu.topgenes = ln.scaled.list.clus$seu.topgenes
savestring = 'ssc_ca3_ln_scaled_ln'
top.five.genes = ln.scaled.list.clus$top.five.genes

pdf(paste0(results.dir, '/tsne_topgenes_', savestring, '.pdf'), height = 6, 
    width = 6)
DimPlot(object = seu.topgenes, reduction = "tsne",
        group.by = 'subclass_label')
dev.off()

pdf(paste0(results.dir, '/umap_topgenes_', savestring, '.pdf'), height = 6, 
    width = 6)
DimPlot(object = seu.topgenes, 
        reduction = "umap", group.by = 'subclass_label')
dev.off()

pdf(paste0(results.dir, '/umap_clus_topgenes_', savestring, '.pdf'), height = 6, 
    width = 6)
DimPlot(object = seu.topgenes, 
        reduction = "umap", group.by = 'seurat_clusters')
dev.off()
pdf(paste0(results.dir, '/featureplot_genes_', savestring, '.pdf'), height = 6, 
    width = 15)
FeaturePlot(seu.topgenes, features = ca3_de, ncol = 5)
dev.off()

i = 1
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 2
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 3
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 4
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 5
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 6
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 7
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 8
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 9
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()
i = 10
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', names(top.five.genes[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes, features = names(top.five.genes[[i]]), ncol = 1)
dev.off()

match.idx = match(ca3_de, names(unlist(lapply(top.five.genes, function(genes) genes[1]))))
top.five.genes.ord = top.five.genes[match.idx]

pdf(paste0(results.dir, '/featt_all3_', savestring, '.pdf'), height = length(ca3_de)*3, 
    width = 9)
FeaturePlot(seu.topgenes, features = names(unlist(lapply(top.five.genes.ord, 
                                                         function(genes) genes[1:3]))), ncol = 3)
# }
dev.off()

pdf(paste0(results.dir, '/featt_all4_', savestring, '.pdf'), height = length(ca3_de)*3, 
    width = 12)
FeaturePlot(seu.topgenes, features = names(unlist(lapply(top.five.genes.ord, 
                                                         function(genes) genes[1:4]))), ncol = 4)
dev.off()

# --------------------# --------------------# --------------------
# add violin plots and compute p values using wilcoxon rank sum test. 
# ident 3 is CA3, and idents 0:4 \ 3 are the SSC neurons

CA3.color = '#E45825'
SSC.color = '#664621'

top.three.genes = lapply(top.five.genes, function(x) x[1:3])
top.three.genes.unlisted = names(unlist(top.three.genes))
top.three.genes.unduplicated = top.three.genes.unlisted[!duplicated(top.three.genes.unlisted)]

top.three.genes.ssc = lapply(top.five.genes[1:5], function(x) x[1:3])
top.three.genes.ssc.unlisted = names(unlist(top.three.genes.ssc))
top.three.genes.ssc.unduplicated = top.three.genes.ssc.unlisted[!duplicated(top.three.genes.ssc.unlisted)]

top.three.genes.ca3 = lapply(top.five.genes[6:10], function(x) x[1:3])
top.three.genes.ca3.unlisted = names(unlist(top.three.genes.ca3))
top.three.genes.ca3.unduplicated = top.three.genes.ca3.unlisted[!duplicated(top.three.genes.ca3.unlisted)]

library(plyr)
# setdiff(rownames(bank_OD_mk_scaled@own.expr), OD.mk.de)
# setdiff(OD.mk.de, rownames(bank_OD_mk_scaled@own.expr))

gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[top.three.genes.ssc.unduplicated,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes)) # ident 3 is CA3, and idents 0:4 \ 3 are the SSC neurons
# We will change ident 3 CA3 to 0 and idents 0,1,2,4 ssc to 1
gene_data_od[, Cluster := as.numeric(Cluster %in% as.factor(c(0, 1, 2, 4)))]

gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('SSC','CA3'))
len <- length(top.three.genes.ssc.unduplicated)
wilcox.ssc <- vector(mode = "list", length = len)
for (i in 1:length(wilcox.ssc)){
  print(i)
  wilcox.ssc[[i]] = wilcox.test(unname(unlist(gene_data_od[gene_data_od$Cluster=='CA3',
                                                           top.three.genes.ssc.unduplicated[i], with = FALSE])),
                                unname(unlist(gene_data_od[gene_data_od$Cluster=='SSC',
                                                           top.three.genes.ssc.unduplicated[i], with = FALSE])),
                                paired = FALSE,alternative = 'l',
                                exact = FALSE)
  print(str(wilcox.ssc[[i]]))
}
names(wilcox.ssc)<-top.three.genes.ssc.unduplicated
p.values.ssc = unlist(lapply(wilcox.ssc, function(x) x$p.value))
print(p.values.ssc)

# Sparcl1          Pcp4         Sept7          Cd34        Igfbp6          Pdp1          Grm3         Ptprd          Fut9 
# 4.645920e-96 1.725468e-154 4.430863e-164  7.203888e-28  3.140053e-34  4.278194e-14 2.300482e-161 2.992883e-171 3.033771e-154 
# Egr1           Arc          Junb          Tbr1         Ttc9b 3110035E14Rik 
# 2.253887e-15  4.825409e-01  9.999999e-01 1.603674e-104 1.142874e-156 1.330724e-114 
sig.ssc = names(p.values.ssc[p.values.ssc<0.05])


gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[sig.ssc,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes)) # ident 3 is CA3, and idents 0:4 \ 3 are the SSC neurons
# We will change ident 3 CA3 to 0 and idents 0,1,2,4 ssc to 1
gene_data_od[, Cluster := as.numeric(Cluster %in% as.factor(c(0, 1, 2, 4)))]

gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('SSC','CA3'))
od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= CA3.color, '1' = SSC.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_sscca3_SSChigh_', 'scaled' ,'.pdf'), 
       height = 6, width = 6)

gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[top.three.genes.ca3.unduplicated,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes)) # ident 3 is CA3, and idents 0:4 \ 3 are the SSC neurons
# We will change ident 3 CA3 to 0 and idents 0,1,2,4 ssc to 1
gene_data_od[, Cluster := as.numeric(Cluster %in% as.factor(c(0, 1, 2, 4)))]

gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('SSC','CA3'))


len <- length(top.three.genes.ca3.unduplicated)
wilcox.ca3 <- vector(mode = "list", length = len)
for (i in 1:length(wilcox.ca3)){
  print(i)
  wilcox.ca3[[i]] = wilcox.test(unname(unlist(gene_data_od[gene_data_od$Cluster=='CA3',
                                                           top.three.genes.ca3.unduplicated[i], with = FALSE])),
                                unname(unlist(gene_data_od[gene_data_od$Cluster=='SSC',
                                                           top.three.genes.ca3.unduplicated[i], with = FALSE])),
                                paired = FALSE,alternative = 'g',
                                exact = FALSE)
  print(str(wilcox.ca3[[i]]))
}
names(wilcox.ca3)<-top.three.genes.ca3.unduplicated
p.values.ca3 = unlist(lapply(wilcox.ca3, function(x) x$p.value))
print(p.values.ca3)
print(p.values.ca3)
# Nefl          Ncdn        Pcdh20       Cacna1a          Meg3         Nlgn1         Cadm3         Epha6        Rnf182 
# 1.353047e-164 6.525021e-172  0.000000e+00  3.230947e-25  1.017147e-01  1.331388e-85  2.582151e-77 8.627087e-185 5.413125e-213 
# Clstn2          Tcf4         Cabp1          Gnaq 4930547E14Rik 
# 1.272372e-63 3.076018e-156  2.023537e-57  2.720620e-90 5.667812e-281 

sig.ca3 = names(p.values.ca3[p.values.ca3<0.05])

gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[sig.ca3,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes)) # ident 3 is CA3, and idents 0:4 \ 3 are the SSC neurons
# We will change ident 3 CA3 to 0 and idents 0,1,2,4 ssc to 1
gene_data_od[, Cluster := as.numeric(Cluster %in% as.factor(c(0, 1, 2, 4)))]

gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('SSC','CA3'))


od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= CA3.color, '1' = SSC.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_sscca3_CA3high_', 'scaled' ,'.pdf'), 
       height = 6, width = 6)




# --------------------# --------------------# --------------------

# --------------------# --------------------# --------------------

# --------------------# --------------------# --------------------

# --------------------# --------------------# --------------------
ss.seurat.OLIGO_unprocessed<-readRDS(file = paste0(data.dir2, '/ss_seurat_oligo_unprocessed.rds'))
ss.seurat.OLIGO.ln <- NormalizeData(object = ss.seurat.OLIGO_unprocessed, normalization.method = 'LogNormalize')
ss.seurat.OLIGO.ln.scaled <- ScaleData(object = ss.seurat.OLIGO.ln)
ss.seurat.OLIGO.unnorm.scaled <- ScaleData(object = ss.seurat.OLIGO_unprocessed)

fimbria_de <- readRDS('/home/ubuntu/banksy/VeraFISH/Oct25_Aug31_Jun8_dt0728_r1_5_0_L0_25_seed42/fimbria_de_Oct25.rds')



data.unnorm <- GetAssayData(ss.seurat.OLIGO_unprocessed, slot = 'counts')
data.counts = data.unnorm

seedgenes.OLIGO<-fimbria_de[fimbria_de %in% rownames(data.unnorm)]
A1 = as.matrix(data.unnorm[seedgenes.OLIGO,])
B1 = A1>0
c1 = rowSums(B1)
c1[1:10]

c1>cutoff*ncol(A1) # bottom 2 percent- 

genes.to.keep = names(which(c1>cutoff*ncol(A1))) # remove bottom 2% of genes
# > genes.to.keep
seedgenes.OLIGO<-genes.to.keep





# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
# --------------------# --------------------# --------------------
#
#
#
#
#
#

# --------------- Fimbria and Thalamic oligodendrocytes -------------------
results.dir <- paste0('fig4-hippocampus/out/scrna_oligos')

check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)

ss.seurat.OLIGO_unprocessed<-readRDS(file = paste0('fig4-hippocampus/data/ss_seurat_oligo_unprocessed.rds'))
ss.seurat.OLIGO.ln <- NormalizeData(object = ss.seurat.OLIGO_unprocessed, normalization.method = 'LogNormalize')
ss.seurat.OLIGO.ln.scaled <- ScaleData(object = ss.seurat.OLIGO.ln)
fimbria_de <- readRDS('fig4-hippocampus/data/fimbria_de.rds')
data.counts <- GetAssayData(ss.seurat.OLIGO_unprocessed, slot = 'counts')

fimbria_de<-fimbria_de[fimbria_de %in% rownames(data.counts)]
A1 = as.matrix(data.counts[fimbria_de,])
B1 = A1>0
c1 = rowSums(B1)
c1[1:10]
fimbria_de = names(which(c1>cutoff*ncol(A1))) # remove bottom 2% of genes

data.ln.scaled <- GetAssayData(ss.seurat.OLIGO.ln.scaled, slot = 'scale.data')
fimbria_de<-fimbria_de[fimbria_de %in% rownames(data.ln.scaled)]

cor.genes.OLIGO <- lapply(fimbria_de, find_corr_genes, data.mat = data.ln.scaled)
top.genes.list.OLIGO <- lapply(cor.genes.OLIGO, function(x) x[1:num_top_genes_to_use])

top.genes.list.OLIGO.filtered = lapply(top.genes.list.OLIGO, remove.lowly.expressed.genes)

top.genes <- unique((unlist(top.genes.list.OLIGO.filtered)))
seu.topgenes.OLIGO <- subset(x = ss.seurat.OLIGO_unprocessed, features = top.genes)

ln.scaled.list.clus.OLIGO = cluster.and.plot(seu.topgenes.OLIGO, top.genes.list.OLIGO.filtered, 
                                             normalization='LogNormalize') 

seu.topgenes.OLIGO = ln.scaled.list.clus.OLIGO$seu.topgenes
savestring = 'OLIGO_ln_scaled_ln'
top.five.genes.OLIGO = ln.scaled.list.clus.OLIGO$top.five.genes

pdf(paste0(results.dir, '/tsne_topgenes_', savestring, '.pdf'), height = 3, 
    width = 3)
DimPlot(object = seu.topgenes.OLIGO, reduction = "tsne",
        group.by = 'subclass_label')
dev.off()
pdf(paste0(results.dir, '/umap_topgenes_', savestring, '.pdf'), height = 3, 
    width = 3)
DimPlot(object = seu.topgenes.OLIGO, 
        reduction = "umap", group.by = 'subclass_label')
dev.off()
pdf(paste0(results.dir, '/umap_clus_topgenes_', savestring, '.pdf'), height = 3, 
    width = 3.5)
DimPlot(object = seu.topgenes.OLIGO, 
        reduction = "umap", group.by = 'seurat_clusters')
dev.off()
pdf(paste0(results.dir, '/featureplot_genes_', savestring, '.pdf'), height = 6, 
    width = 15)
FeaturePlot(seu.topgenes.OLIGO, features = fimbria_de, ncol = 5)
dev.off()
i = 1
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 2
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 3
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 4
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 5
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 6
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 7
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 8
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 9
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()
i = 10
pdf(paste0(results.dir, '/feat_seed_', savestring, '_', (top.five.genes.OLIGO[[i]])[1], '.pdf'), height = 15, 
    width = 4)
FeaturePlot(seu.topgenes.OLIGO, features = (top.five.genes.OLIGO[[i]]), ncol = 1)
dev.off()

match.idx = match(fimbria_de, (unlist(lapply(top.five.genes.OLIGO, function(genes) genes[1]))))
top.five.genes.OLIGO.ord = top.five.genes.OLIGO[match.idx]

pdf(paste0(results.dir, '/featt_all3_', savestring, '.pdf'), height = length(fimbria_de)*3, 
    width = 9)
FeaturePlot(seu.topgenes.OLIGO, features = (unlist(lapply(top.five.genes.OLIGO.ord, 
                                                          function(genes) genes[1:3]))), ncol = 3)
dev.off()

pdf(paste0(results.dir, '/featt_all4_', savestring, '.pdf'), height = length(fimbria_de)*3, 
    width = 12)
FeaturePlot(seu.topgenes.OLIGO, features = (unlist(lapply(top.five.genes.OLIGO.ord, 
                                                          function(genes) genes[1:4]))), ncol = 4)
dev.off()


# violin plots and wilcoxon test: change to oligos (code copied over from ca3/ssc)

# add violin plots and compute p values using wilcoxon rank sum test. 
# ident 0 is fimbra, and idents 1 are the thalamic oligos neurons

fimbria.color = '#008957'
thalamic.color = '#8E56A2'

top.three.genes = lapply(top.five.genes.OLIGO.ord, function(x) x[1:3])
top.three.genes.unlisted = (unlist(top.three.genes))
top.three.genes.unduplicated = top.three.genes.unlisted[!duplicated(top.three.genes.unlisted)]

top.three.genes.fimbria = lapply(top.five.genes.OLIGO.ord[1:3], function(x) x[1:3])
top.three.genes.fimbria.unlisted = (unlist(top.three.genes.fimbria))
top.three.genes.fimbria.unduplicated = top.three.genes.fimbria.unlisted[!duplicated(top.three.genes.fimbria.unlisted)]

top.three.genes.thalamic = lapply(top.five.genes.OLIGO.ord[4:9], function(x) x[1:3])
top.three.genes.thalamic.unlisted = (unlist(top.three.genes.thalamic))
top.three.genes.thalamic.unduplicated = top.three.genes.thalamic.unlisted[!duplicated(top.three.genes.thalamic.unlisted)]

library(plyr)

len <- length(top.three.genes.fimbria.unduplicated)
gene_data_od <- data.table(t(as.matrix(seu.topgenes.OLIGO@assays$RNA@scale.data[top.three.genes.fimbria.unduplicated,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes.OLIGO)) 
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('Thalamic OD','Fornix OD'))
wilcox.fimbria <- vector(mode = "list", length = len)
for (i in 1:length(wilcox.fimbria)){
  print(i)
  wilcox.fimbria[[i]] = wilcox.test(unname(unlist(gene_data_od[gene_data_od$Cluster=='Thalamic OD',
                                                               top.three.genes.fimbria.unduplicated[i], with = FALSE])),
                                    unname(unlist(gene_data_od[gene_data_od$Cluster=='Fornix OD',
                                                               top.three.genes.fimbria.unduplicated[i], with = FALSE])),
                                    paired = FALSE,alternative = 'l',
                                    exact = FALSE)
  print(str(wilcox.fimbria[[i]]))
}
names(wilcox.fimbria)<-top.three.genes.fimbria.unduplicated
p.values.fimbria = unlist(lapply(wilcox.fimbria, function(x) x$p.value))
print(p.values.fimbria)
# Mobp          Mog        Sept4         Plp1          Mbp          Cnp         Gfap       Etnppl      Gm33370 
# 1.314516e-38 1.288091e-37 1.348884e-37 2.316211e-35 2.854165e-33 2.310263e-29 5.935266e-02 5.787939e-01 3.025504e-01 
# setdiff(rownames(bank_OD_mk_scaled@own.expr), OD.mk.de)
# setdiff(OD.mk.de, rownames(bank_OD_mk_scaled@own.expr))
sig.fimbria = names(p.values.fimbria[p.values.fimbria<0.05])

gene_data_od <- data.table(t(as.matrix(seu.topgenes.OLIGO@assays$RNA@scale.data[sig.fimbria,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes.OLIGO)) 
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('Thalamic OD','Fornix OD'))
od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= fimbria.color, '1' = thalamic.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_fimbriathalamic_Fornixhigh_', 'scaled' ,'.pdf'), 
       height = 6, width = 6)



#-----


len <- length(top.three.genes.thalamic.unduplicated)
gene_data_od <- data.table(t(as.matrix(seu.topgenes.OLIGO@assays$RNA@scale.data[top.three.genes.thalamic.unduplicated,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes.OLIGO)) 
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('Thalamic OD','Fornix OD'))
wilcox.thalamic <- vector(mode = "list", length = len)
for (i in 1:length(wilcox.thalamic)){
  print(i)
  wilcox.thalamic[[i]] = wilcox.test(unname(unlist(gene_data_od[gene_data_od$Cluster=='Thalamic OD',
                                                                top.three.genes.thalamic.unduplicated[i], with = FALSE])),
                                     unname(unlist(gene_data_od[gene_data_od$Cluster=='Fornix OD',
                                                                top.three.genes.thalamic.unduplicated[i], with = FALSE])),
                                     paired = FALSE,alternative = 'g',
                                     exact = FALSE)
  print(str(wilcox.thalamic[[i]]))
}
names(wilcox.thalamic)<-top.three.genes.thalamic.unduplicated
p.values.thalamic = unlist(lapply(wilcox.thalamic, function(x) x$p.value))
print(p.values.thalamic)
# Bcas1          Bmp4         Mpzl1         Mapk9         Kdm7a C330011M18Rik       Sparcl1         C1ql1      Serpine2 
# 2.066564e-04  4.016780e-13  1.463489e-06  1.424220e-02  7.300949e-03  4.486893e-01  1.340790e-34  7.221627e-44  8.846435e-43 
# Atp1a2        Ptprz1        Atp1b2         Cspg5        Pdgfra          Nefl         Acta1          Dmkn 
# 4.790698e-36  8.451538e-39  5.253555e-23  1.688320e-34  2.112048e-40  6.901452e-01  7.158629e-01  1.000000e+00 

sig.thalamic = names(p.values.thalamic[p.values.thalamic<0.05])
# setdiff(rownames(bank_OD_mk_scaled@own.expr), OD.mk.de)
# setdiff(OD.mk.de, rownames(bank_OD_mk_scaled@own.expr))
gene_data_od <- data.table(t(as.matrix(seu.topgenes.OLIGO@assays$RNA@scale.data[sig.thalamic,])), #[OD.mk.de,]
                           Cluster = Idents(seu.topgenes.OLIGO)) 
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('Thalamic OD','Fornix OD'))
od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= fimbria.color, '1' = thalamic.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_fimbriathalamic_Thalamushigh_', 'scaled' ,'.pdf'), 
       height = 6, width = 6)

