# You need the (legacy) version of BANKSY (v 0.1.5) for recreating this analysis: 
# remotes::install_github("prabhakarlab/Banksy@main") # to install 0.1.5. 

library(Seurat)
library(scater)
library(scran)
library(data.table)
library(plyr) # mapvalues
library(gridExtra)
cutoff = 0.01
options(timeout=1200)

# these genes are from the MERFISH analysis (see src/process_hypothalamus_merfish.R)
DE_genes_raw =c("Mlc1","Dgkk","Cbln2","Syt4","Gad1","Plin3","Gnrh1",
                "Sln","Gjc3","Mbp","Lpar1","Trh","Ucn3","Cck") 

# https://satijalab.org/seurat/reference/readmtx
if (TRUE) { # set to TRUE to generate seurat object. This step can take a while (large file download)
expression_matrix <- ReadMtx(mtx = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_matrix.mtx.gz",
                             cells = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_barcodes.tsv.gz",
                             features = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_genes.tsv.gz")
ss.seurat <- CreateSeuratObject(counts = expression_matrix)
} # sometimes you get errors like 'Server returned nothing (no headers, no data)' when downloading. 
# We find just restarting R, or checking internet connection (etc) usually fixes this.. Not sure what causes the error. 
# in any case, you can manually download the files from the links above, and use ReadMtx on the local files to generate the 
# Seurat object. 

results.dir = 'fig3-hypothalamus/out/scrna/'
meta <- openxlsx::read.xlsx('fig3-hypothalamus/data/aau5324_Moffitt_Table-S1.xlsx', 
                            colNames = TRUE)

stopifnot(all(meta$Cell.name==names(Idents(ss.seurat))))
names(meta) <- c('Cell_name', 'Sex', 'Rep', 'Cell_class',
                 'Non_neuronal_cluster', 'Neuronal_cluster')
ss.seurat@meta.data <- cbind(ss.seurat@meta.data, 
                             meta[match(meta$Cell_name, 
                                        rownames(ss.seurat@meta.data)),])
Idents(ss.seurat) <- ss.seurat@meta.data$Cell_class

# -------- functions ---------- #
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

cluster.with.topgenes <- function(seu.topgenes, 
                                  normalization=NULL
                                  # 'RC', 'LogNormalize'
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
  seu.topgenes <- FindClusters(object = seu.topgenes, resolution = 0.05)
  seu.topgenes <- RunTSNE(object = seu.topgenes, dims = 1:npcs)
  seu.topgenes <- Seurat:::RunUMAP(object = seu.topgenes, dims = 1:npcs)
  return(seu.topgenes)
}
remove.genes <- function(x){
  Ax = as.matrix(data.counts[x,]) 
  Bx = Ax>0
  cx = rowSums(Bx)
  x.to.keep = names(which(cx>cutoff*ncol(Ax)))
  return(x.to.keep)
}
# -------- /functions ---------- #

ss.seurat[['pct.mito']] <- PercentageFeatureSet(ss.seurat, '^mt')
bqc <- FeatureScatter(ss.seurat, feature2 = 'pct.mito', feature1 = 'nFeature_RNA') 
ss.seurat <- subset(ss.seurat, subset = pct.mito < 20 & nFeature_RNA > 1000)
aqc <- FeatureScatter(ss.seurat, feature2 = 'pct.mito', feature1 = 'nFeature_RNA')
QC_ALL <- bqc + aqc

ss.seurat.oligo = subset(x = ss.seurat, subset = Cell_class=="Mature oligodendrocyte")
ss.seurat.oligo.ln <- NormalizeData(object = ss.seurat.oligo, 
                                     normalization.method = 'LogNormalize')
ss.seurat.oligo.ln <- ScaleData(object = ss.seurat.oligo.ln)


data.scaled <- GetAssayData(ss.seurat.oligo.ln, slot = 'scale.data')
data.counts <- GetAssayData(ss.seurat.oligo.ln, slot = 'counts')
genes.to.keep = remove.genes(DE_genes_raw) 
genes.to.keep <- genes.to.keep[genes.to.keep %in% rownames(data.scaled)] 
cor.genes <- lapply(genes.to.keep, find_corr_genes, data.mat = data.scaled)
top.genes.list <- lapply(cor.genes, function(x) names(x[1:25]))
top.genes.list <- top.genes.list[match(genes.to.keep, 
                                       (unlist(lapply(top.genes.list, 
                                                           function(x) x[1]))))]

top.genes.list = lapply(top.genes.list, remove.genes)
top.genes.undup = (unlist(top.genes.list))
top.genes.undup = top.genes.undup[!duplicated(top.genes.undup)]

# subset to create Seurat object
seu.topgenes.unprocessed <- subset(x = ss.seurat.oligo, features = top.genes.undup)
seu.topgenes = cluster.with.topgenes(seu.topgenes.unprocessed, 
                                       normalization='LogNormalize') 

dense.color = '#BE0032'
sparse.color = '#F38400'
color.map = c('0' = dense.color, '1' = sparse.color)

pdf(paste0(results.dir, '/umap_clusters', '.pdf'), height = 6, width = 6)
DimPlot(object = seu.topgenes, cols = color.map,
        reduction = "umap", group.by = 'seurat_clusters')
dev.off()

top.3.sparse = lapply(top.genes.list[c(1:8)], function(x){
  if (length(x)>=3){
    return(x[1:3])
  } else {
    return(x[1:length(x)])
  }
})
top.3.sparse.unlisted = (unlist(top.3.sparse))
top.3.sparse.undup = top.3.sparse.unlisted[!duplicated(top.3.sparse.unlisted)]

top.3.dense = lapply(top.genes.list[9:12], function(x) x[1:3])
top.3.dense.unlisted = (unlist(top.3.dense))
top.3.dense.undup = top.3.dense.unlisted[!duplicated(top.3.dense.unlisted)]

###------- violin
gene_data_od.sparse <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[top.3.sparse.undup,])), 
                                  Cluster = Idents(seu.topgenes))
gene_data_od.dense <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[top.3.dense.undup,])), 
                                 Cluster = Idents(seu.topgenes))

len <- length(top.3.sparse.undup)
wilcox.sparse <- vector(mode = "list", length = len)

len <- length(top.3.dense.undup)
wilcox.dense <- vector(mode = "list", length = len)

for (i in 1:length(wilcox.sparse)){
  wilcox.sparse[[i]] = wilcox.test(unname(unlist(gene_data_od.sparse[Cluster==0,
                                                                     top.3.sparse.undup[i], with = FALSE])),
                                   unname(unlist(gene_data_od.sparse[Cluster==1,
                                                                     top.3.sparse.undup[i], with = FALSE])),
                                   paired = FALSE,alternative = 'l',
                                   exact = TRUE)
}
names(wilcox.sparse)<-top.3.sparse.undup
p.values.sparse = unlist(lapply(wilcox.sparse, function(x) x$p.value))
print(p.values.sparse[p.values.sparse<0.05])

#          Mlc1        Atp1a2        Slc4a4          Dgkk        Elavl4         Tenm1         Cbln2         Cxxc4          Grm1          Syt4 
# 7.209978e-07  9.511541e-15  2.832355e-08  1.502088e-97 5.085965e-272  0.000000e+00  1.525389e-64 1.934872e-213 1.348877e-185 2.565480e-170 
# Meg3        Snhg11          Gad1         S100b 
# 1.933490e-173 2.175790e-183 3.422688e-128  2.933241e-07 

for (i in 1:length(wilcox.dense)){
  wilcox.dense[[i]] = wilcox.test(unname(unlist(gene_data_od.dense[Cluster==0,
                                                                   top.3.dense.undup[i], 
                                                                   with = FALSE])),
                                  unname(unlist(gene_data_od.dense[Cluster==1,
                                                                   top.3.dense.undup[i], 
                                                                   with = FALSE])),
                                  paired = FALSE,alternative = 'g',
                                  exact = TRUE)
}
names(wilcox.dense)<-top.3.dense.undup
p.values.dense = unlist(lapply(wilcox.dense, function(x) x$p.value))
print(p.values.dense[p.values.dense<0.05])
# Mbp         Mobp         Plp1        Lpar1       Cldn11 
# 1.314820e-17 1.246486e-06 3.361609e-22 3.782918e-03 8.168704e-12 

sig.sparse = names(p.values.sparse[p.values.sparse<0.05])
gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[sig.sparse,])), 
                           Cluster = Idents(seu.topgenes))
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('dense ODM', 'sparse ODM'))
od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= dense.color, '1' = sparse.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_OD_sparsehigh_', 'scaled' ,'.pdf'), 
       height = 6, width = 6)


sig.dense = names(p.values.dense[p.values.dense<0.05])
gene_data_od <- data.table(t(as.matrix(seu.topgenes@assays$RNA@scale.data[sig.dense,])), 
                           Cluster = Idents(seu.topgenes))
gene_data_od$Cluster <- mapvalues(gene_data_od[,Cluster], 
                                  from=unique(gene_data_od$Cluster), 
                                  to=c('dense ODM', 'sparse ODM'))
od_melted <- melt(gene_data_od, id.vars = "Cluster")
colnames(od_melted)[c(2,3)]<-c('Gene', 'Expression')
ggplot(od_melted, aes(x=Gene, y=Expression, fill=Cluster)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=unname(c('0'= dense.color, '1' = sparse.color)))+
  geom_violin(trim = TRUE) + facet_wrap(~Gene, scales = "free_x")+
  geom_boxplot(aes(x=Gene, y=Expression), width=0.1, 
               outlier.shape = NA, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2, 4))

ggsave(paste0(results.dir, '/vplot_OD_densehigh_', 'scaled' ,'.pdf'), 
       height = 4.5, width = 6)