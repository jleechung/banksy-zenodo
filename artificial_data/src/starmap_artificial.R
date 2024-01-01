source('artificial_data/src/process_artificial.R')
for (ngene in c(400, 600, 800, 1020)){
  library(scDesign2)
  library(mclust)
  library(Rtsne)
  library(scran)
  library(plyr)      # mapvalues
  library(reshape2)  # melt
  library(gridExtra) # arrangeGrob
  library(Seurat)
  library(data.table)
  library(cowplot)   # draw_plot_label
  library(ggplot2); theme_set(theme_bw());
  
  # You need the (legacy) version of BANKSY (v 0.1.5) for recreating this analysis: 
  # remotes::install_github("prabhakarlab/Banksy@main") # to install 0.1.5. 
  library(Banksy)
  
  
  data.dir = 'artificial_data/data/'
  check <- dir.exists('artificial_data/out/')
  if (!check) dir.create('artificial_data/out/')
  
  data.starmap.old.meta <- read.csv(file = paste0(data.dir, 'Starmap_BY3_1k_meta_annotated.csv'))
  data.starmap.meta <- read.csv(file = paste0(data.dir, 'Starmap_BY3_1k_meta_annotated_14oct22.csv'))
  data.starmap.gcm <- read.csv(file = paste0(data.dir, 'Starmap_BY3_1k_gcm.csv'))
  data.starmap.gcm.raw <- read.csv(file = paste0(data.dir, 'Starmap_BY3_1k_gcm_raw.csv'))
  data.starmap.clusters_p1 <- read.csv(file = paste0(data.dir, 'clusters_p_1_starmap.csv'))
  
  celltype.prop.perturb = 0
  scale.factor = 2
  n.cell = 0.1
  num.clus = 16
  
  # ------------------------------------------------------------------------------
  cell_ids = paste0('cell_', data.starmap.gcm.raw$X)
  data.starmap.gcm.raw.temp = data.starmap.gcm.raw[, 2:ncol(data.starmap.gcm.raw)]
  matrix.starmap.gcm.raw = t(data.starmap.gcm.raw.temp)
  colnames(matrix.starmap.gcm.raw) = cell_ids
  cell_types = names(table(data.starmap.old.meta$cluster_name))
  all_counts = matrix.starmap.gcm.raw
  colnames(all_counts)<- data.starmap.old.meta$cluster_name
  
  all.genes = rownames(matrix.starmap.gcm.raw)
  
  genes.subset.starmap.158 = c('Acss1', 'Adcyap1', 'Adgrl2', 'Aqp4', 'Arc', 'Arf5', 'Arhgap24', 
                               'Arl4d', 'Arx', 'Batf3', 'Bcl6', 'Bdnf', 'Bgn', 'Btg2', 'Calb2', 
                               'Car12', 'Car4', 'Cbln4', 'Cck', 'Cdh13', 'Cdk6', 'Chat', 'Chodl',
                               'Chrna6', 'Col6a1', 'Cplx3', 'Cpne5', 'Crh', 'Crispld2', 'Csrnp1',
                               'Ctgf', 'Ctxn3', 'Cux2', 'Cxcl14', 'Ddit4l', 'Deptor', 'Egr1', 
                               'Egr2', 'Egr4', 'Enpp2', 'eRNA1', 'eRNA2', 'eRNA3', 'eRNA4',
                               'eRNA5', 'Fam19a1', 'Flt1', 'Fos', 'Fosb', 'Fosl1', 'Fosl2', 
                               'Foxp2', 'Frmd7', 'Frmpd3', 'Gad1', 'Gad2', 'Gda', 'Gpc3', 
                               'Gpr3', 'Gpx3', 'Homer1', 'Hsd11b1', 'Htr3a', 'Igf1', 'Igtp', 
                               'Itgam', 'Junb', 'Kcna1', 'Kcnip2', 'Kcnk2', 'Lhx6', 'Lmo2',
                               'Mab21l1', 'Mdga1', 'Mgp', 'Mog', 'Mybpc1', 'Myh8', 'Myl4',
                               'Mylk', 'Myo1c', 'Ndnf', 'Nectin3', 'Nectin4', 'Ngb', 'Nos1',
                               'Nov', 'Npas4', 'Nptx2', 'Npy', 'Nr4a2', 'Nrn1', 'Obox3', 
                               'Osgin2', 'Otof', 'Parm1', 'Pax6', 'Pcdhac1', 'Pcdhac2', 
                               'Pcdhgc3', 'Pcdhgc4', 'Pcdhgc5', 'Pcp4', 'Pde1a', 'Pde1c',
                               'Pdgfra', 'Pdzrn3', 'Penk', 'Plcxd2', 'Pnoc_01', 'Ppm1h',
                               'Prok2', 'Prss22', 'Ptgs2', 'Pthlh', 'Pvalb', 'Qrfpr', 
                               'Rasgrf2', 'Reln', 'Rerg', 'Rgs10', 'Rgs12', 'Rgs2', 'Rorb', 
                               'Rprm', 'Rspo2', 'Scnn1a', 'Sema3c', 'Sema3e', 'Serpinb11', 
                               'Serpine1', 'Sla', 'Slc17a7', 'Slc25a36', 'Slc5a7', 'Smad3', 
                               'Sncg', 'Spp1', 'Sst', 'Stac', 'Sulf2', 'Synpr', 'Syt17', 'Syt6', 
                               'Tacr1', 'Tacr3', 'Tacstd2', 'Tbr1', 'Tcerg1l', 'Th', 'Thsd7a', 
                               'Tnfaip6', 'Tnfaip8l3', 'Tnmd', 'Tpbg', 'Tph2', 'Ucma', 'Vgf', 
                               'Vip', 'Wt1')
  
  
  
  remaining.genes = all.genes[which(!(all.genes %in% genes.subset.starmap.158))]
  
  set.seed(42)
  randgenes_48 = sample(remaining.genes, 48)
  set.seed(43)
  randgenes_248 = sample(remaining.genes, 248)
  set.seed(44)
  randgenes_448 = sample(remaining.genes, 448)
  set.seed(45)
  randgenes_648 = sample(remaining.genes, 648)
  
  switch(as.character(ngene), 
         '400'={
           genes.subset.starmap = union(genes.subset.starmap.158, randgenes_248)
           results.dir = 'artificial_data/out/ngene_400/'
           check <- dir.exists(results.dir)
           if (!check) dir.create(results.dir)
         },
         '600'={
           genes.subset.starmap = union(genes.subset.starmap.158, randgenes_448)
           results.dir = 'artificial_data/out/ngene_600/'
           check <- dir.exists(results.dir)
           if (!check) dir.create(results.dir)
         },
         '800'={
           genes.subset.starmap = union(genes.subset.starmap.158, randgenes_648)
           results.dir = 'artificial_data/out/ngene_800/'
           check <- dir.exists(results.dir)
           if (!check) dir.create(results.dir)
         },
         '1020'={
           # all genes
           genes.subset.starmap = union(genes.subset.starmap.158, remaining.genes)
           results.dir = 'artificial_data/out/ngene_1020/'
           check <- dir.exists(results.dir)
           if (!check) dir.create(results.dir)
         },
         {
           print('The ngene parameter can only be 400, 600, 800 or the full 1020.')
         }
  )
  
  dff = data.starmap.old.meta[, c("X", "domain", "smoothed_manual", "cluster_name", "x", "y")]
  names(dff)[colnames(dff)=='X'] = 'cell_ID'
  dff$cell_ID = paste0('cell_', dff$cell_ID)
  dtt = as.data.table(dff)
  
  
  dtt4 = dtt[, by = .(domain, cluster_name), .N
  ][, per := round(prop.table(N),3) , by = "domain"
  ][order(domain, -N)]
  
  
  dtt.minmax.e = dtt[, .(maxes = scale.factor*max(x), 
                         mins = scale.factor*min(x), 
                         num_cells_rectangles = .N*scale.factor^2), 
                     by = .(domain)][order(-mins)]
  
  dtt.manual = dtt[, .(num_cells_manual = .N*scale.factor^2), 
                   by = .(smoothed_manual)][
                     order(match(smoothed_manual, dtt.minmax.e$domain))]
  
  dtt.domains = cbind(dtt.minmax.e, dtt.manual)
  
  dtt4.ordered = dtt4[order(match(domain, dtt.domains$domain))]
  dtt4.wide = reshape(dtt4.ordered[, c('domain', 'cluster_name', 'per')], 
                      idvar = "domain", 
                      timevar = "cluster_name", 
                      direction = 'wide')
  dtt4.wide[is.na(dtt4.wide)]=0
  colnames(dtt4.wide) = gsub('^per.', '', names(dtt4.wide))
  
  dtt.domains[, num_cells_rand := round(runif(7, 1-n.cell, 1+n.cell)*num_cells_rectangles), ]
  dtt.domains2 = cbind(dtt.domains, 
                       round(dtt.domains$num_cells_rand*dtt4.wide[order(match(domain, 
                                                                              dtt.domains$domain)),
                                                                  -c("domain")]))
  
  
  barplot(colSums(dtt.domains2[,-c(1:6)])[-1]/sum(colSums(dtt.domains2[,-c(1:6)])[-1]))
  
  dtt.domains2$num_cells_rand = rowSums(dtt.domains2[,-c(1:7)])
  colnames(dtt.domains2) <- gsub('^per.', '', names(dtt.domains2)) # remove per.
  total_cells_per_category = colSums(dtt.domains2[,-c(1,2,3,5)])
  cell_type_prop2 = total_cells_per_category[-c(1,2,3)]
  cell_type_prop3 = (cell_type_prop2)[match( cell_types, names(cell_type_prop2))]
  copula_result.all <- fit_model_scDesign2(all_counts,
                                           cell_types,
                                           sim_method = 'copula',
                                           ncores = 4)
  saveRDS(copula_result.all, file = paste0(results.dir, 'copula_result.all.rds'))
  copula_result.all2 <- readRDS(file = paste0(results.dir, 'copula_result.all.rds'))
  
  
  
  rand.seed.idx = 1
  lam.celltype = c(0, 0.2)
  lam.domain = 0.8
  
  
  
  
  artificial_seed1 = process.artificial(scale.factor = 2,
                                        rand.seed.idx =1,
                                        lam.celltype = c(0, 0.2),
                                        lam.domain = 0.8,
                                        results.dir = results.dir)
  
  artificial_seed2 = process.artificial(scale.factor = 2,
                                        rand.seed.idx =2,
                                        lam.celltype = c(0, 0.2),
                                        lam.domain = 0.8,
                                        results.dir = results.dir)
  
  artificial_seed3 = process.artificial(scale.factor = 2,
                                        rand.seed.idx =3,
                                        lam.celltype = c(0, 0.2),
                                        lam.domain = 0.8,
                                        results.dir = results.dir)
  # 
  # 
  
  saveRDS(artificial_seed1 , file = paste0(results.dir, '/artificial_seed1.rds'))
  saveRDS(artificial_seed2 , file = paste0(results.dir, '/artificial_seed2.rds'))
  saveRDS(artificial_seed3 , file = paste0(results.dir, '/artificial_seed3.rds'))
  
  
  # -----------------
  
  
  artif.list = readRDS(file = paste0(results.dir, '/artificial_seed1.rds'))
  
  data.matrix = artif.list$bank.unnormalized@own.expr
  locs.df = artif.list$bank.unnormalized@cell.locs
  meta.celltype.df = artif.list$bank.celltype@meta.data
  meta.domain.df = artif.list$bank.domain@meta.data
  
  class(meta.domain.df)
  class(meta.celltype.df)
  class(locs.df)
  class(data.matrix)
  write.csv(data.matrix, file = paste0(results.dir, '/', 'datamatrix_seed1.csv' ))
  write.csv(locs.df, file = paste0(results.dir, '/', 'locsdf_seed1.csv' ))
  write.csv(meta.celltype.df, file = paste0(results.dir, '/', 'metacelltypedf_seed1.csv' ))
  write.csv(meta.domain.df, file = paste0(results.dir, '/', 'metadomaindf_seed1.csv' ))
  
  
  artif.list = readRDS(file = paste0(results.dir, '/artificial_seed2.rds'))
  
  data.matrix = artif.list$bank.unnormalized@own.expr
  locs.df = artif.list$bank.unnormalized@cell.locs
  meta.celltype.df = artif.list$bank.celltype@meta.data
  meta.domain.df = artif.list$bank.domain@meta.data
  
  class(meta.domain.df)
  class(meta.celltype.df)
  class(locs.df)
  class(data.matrix)
  write.csv(data.matrix, file = paste0(results.dir, '/', 'datamatrix_seed2.csv' ))
  write.csv(locs.df, file = paste0(results.dir, '/', 'locsdf_seed2.csv' ))
  write.csv(meta.celltype.df, file = paste0(results.dir, '/', 'metacelltypedf_seed2.csv' ))
  write.csv(meta.domain.df, file = paste0(results.dir, '/', 'metadomaindf_seed2.csv' ))
  
  
  artif.list = readRDS(file = paste0(results.dir, '/artificial_seed3.rds'))
  
  data.matrix = artif.list$bank.unnormalized@own.expr
  locs.df = artif.list$bank.unnormalized@cell.locs
  meta.celltype.df = artif.list$bank.celltype@meta.data
  meta.domain.df = artif.list$bank.domain@meta.data
  
  class(meta.domain.df)
  class(meta.celltype.df)
  class(locs.df)
  class(data.matrix)
  write.csv(data.matrix, file = paste0(results.dir, '/', 'datamatrix_seed3.csv' ))
  write.csv(locs.df, file = paste0(results.dir, '/', 'locsdf_seed3.csv' ))
  write.csv(meta.celltype.df, file = paste0(results.dir, '/', 'metacelltypedf_seed3.csv' ))
  write.csv(meta.domain.df, file = paste0(results.dir, '/', 'metadomaindf_seed3.csv' ))
  
  
}
