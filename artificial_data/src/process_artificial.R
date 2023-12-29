
process.artificial <- function(scale.factor = 2, 
                               rand.seed.idx =1,
                               lam.celltype = c(0, 0.2),
                               lam.domain = 0.8,
                               results.dir = 'test_dir_1'){
  
  set.seed(42*rand.seed.idx)
  sim_count_copula.all_1 <- simulate_count_scDesign2(copula_result.all2,
                                                     sum(cell_type_prop3),
                                                     sim_method = 'copula',
                                                     cell_type_prop = cell_type_prop3)
  
  # check this via a heatmap! To be sure that the matrix generated is correct. 
  rownames(sim_count_copula.all_1)<-rownames(all_counts)
  
  # SUBSETTING-->
  sim_count_copula.all_1=sim_count_copula.all_1[rownames(sim_count_copula.all_1)%in%genes.subset.starmap,]
  
  
  # generate positions for all the cells in the 7 domains
  
  # generate labels
  domains.vec2 = rep(dtt.domains2$domain, dtt.domains2$num_cells_manual)
  
  cell.types.domains.melt = melt(dtt.domains2[,-c(2:7)], 
                                 id.vars = 'domain', 
                                 variable.name = 'cell_type', 
                                 value.name = 'num_cells', 
                                 variable.factor = FALSE)[order(cell_type)] 
  # prevents cell types from being converted to factors. 
  
  domain_celltype_dt = data.table(domain = rep(cell.types.domains.melt$domain, 
                                               cell.types.domains.melt$num_cells), 
                                  cell_types = rep(cell.types.domains.melt$cell_type, 
                                                   cell.types.domains.melt$num_cells))
  domain_celltype_dt$cell_types[!duplicated(domain_celltype_dt$cell_types)]
  # [1] "Astro-1" "Astro-2" "Endo"    "HPC"     "Micro"   "Oligo"   "PVALB"   "Reln"    "SST"    
  # [10] "Smc"     "VIP"     "eL2/3"   "eL4"     "eL5"     "eL6-1"   "eL6-2" 
  colnames(sim_count_copula.all_1)[
    !duplicated(colnames(sim_count_copula.all_1))]
  # [1] "Astro-1" "Astro-2" "eL2/3"   "eL4"     "eL5"     "eL6-1"   "eL6-2"   "Endo"    "HPC"    
  # [10] "Micro"   "Oligo"   "PVALB"   "Reln"    "Smc"     "SST"     "VIP"  
  
  # Next, align the cell types and the expression vectors (ie, the above two arrays):
  domain_celltype_dt = domain_celltype_dt[
    order(match(cell_types,
                colnames(sim_count_copula.all_1)[
                  !duplicated(colnames(sim_count_copula.all_1))]))]
  
  all(domain_celltype_dt$cell_types==colnames(sim_count_copula.all_1))
  
  # Next, generate the domain positions and labels for each cell
  # xcoords = unlist(
  #   lapply(mapply(runif,dtt.domains2$num_cells_rand, 
  #                 dtt.domains2$mins, 
  #                 dtt.domains2$maxes),
  #          function(x) rnorm(length(x), mean = 0,sd = 100 ) + x))
  set.seed(43*rand.seed.idx)
  xcoords = unlist(mapply(runif,dtt.domains2$num_cells_rand, dtt.domains2$mins, dtt.domains2$maxes))
  ycoords = runif(sum(dtt.domains2$num_cells_rand), range(dtt$y)[1], 
                  range(dtt$y)[2]*scale.factor)
  
  # !todo later possibly:
  # generate an alternate set of x and y coords which enforce a minimum distance 
  # between cells (which is more biologically plausible). A simple approach is to 
  # sparsely randomly populate a lattice, where the internode distance enforces the 
  # minimum distance between cells. We can add a little bit of random noise to to the 
  # sampled positions
  
  domain.labels = rep(dtt.domains2$domain, dtt.domains2$num_cells_rand)
  # plot(xcoords, ycoords, cex = 1, pch=20, col = factor(domain.labels))
  domain.coords.dt = data.table(domain = domain.labels, 
                                sdimx = xcoords, 
                                sdimy = ycoords, stringsAsFactors = FALSE)
  domain.coords.dt = domain.coords.dt[order(domain.coords.dt$domain)] # this used radix
  
  all(domain.coords.dt[
    order(domain.coords.dt$domain)]$domain == domain.coords.dt$domain[
      order(domain.coords.dt$domain,method="radix")]) 
  # data.table uses radix. Need to force base R's order to also use radix to get the same results
  # https://stackoverflow.com/questions/72051762/order-in-data-frame-and-data-table
  
  all(domain_celltype_dt$domain[
    order(domain_celltype_dt$domain)] == domain.labels[
      order(domain.labels)]) # here both are not using radix sort, since both are vectors (not data tables)
  
  all(domain_celltype_dt$domain[
    order(domain_celltype_dt$domain)] == domain.coords.dt$domain) 
  # here one is using radix (the data.table order), and one is not (the base R order)
  # so we expect false. 
  
  # Thus, we force the radix sort on the Base R's order function above. 
  
  all(domain_celltype_dt$domain[
    order(domain_celltype_dt$domain, method = 'radix')] == domain.coords.dt$domain) 
  # [1] TRUE
  
  domain_ord_celltype_dt = domain_celltype_dt[order(domain_celltype_dt$domain)]
  all(domain_ord_celltype_dt$domain == domain.coords.dt$domain) 
  # [1] TRUE
  
  domain.coords.dt2 = cbind(domain.coords.dt, 
                            cell_types = domain_ord_celltype_dt$cell_types)
  
  domain.coords.dt3 = domain.coords.dt2[order(match(cell_types, 
                                                    colnames(sim_count_copula.all_1)[
                                                      !duplicated(colnames(sim_count_copula.all_1))]))]
  
  png(paste0(results.dir, '/spatial_ground_truths_raw_seed', 
             rand.seed.idx, '.png'), height = 3*scale.factor, 
      width = 8*scale.factor, res = 300, units = 'in')
  par(mfrow = c(1, 2))   
  pl1 = plot(domain.coords.dt3$sdimx, domain.coords.dt3$sdimy, 
             cex = .6*sqrt(scale.factor), pch=20,
             col = factor(domain.coords.dt3$cell_types), main = 'celltypes')
  
  pl2 = plot(domain.coords.dt3$sdimx, domain.coords.dt3$sdimy, 
             cex = .6*sqrt(scale.factor), pch=20,
             col = factor(domain.coords.dt3$domain), main = 'domain')
  
  # grid.arrange(list(pl1, pl2), ncol = 2)
  # do.call("grid.arrange", c(list(pl1, pl2), ncol = 2))
  dev.off()
  
  domain.celltypes.dt = domain.coords.dt3
  
  # we need the domain and cell type labels to be numeric too, for the purposes of 
  # cluster consensus
  domain.celltypes.dt[,clust_domain := as.numeric(factor(domain.celltypes.dt$domain))]
  domain.celltypes.dt[,clust_celltypes := as.numeric(factor(domain.celltypes.dt$cell_types))]
  
  # final check before constructing banksy object
  all(
    domain.celltypes.dt$cell_types[
      !duplicated(domain.celltypes.dt$cell_types)
    ]
    == 
      colnames(sim_count_copula.all_1)[
        !duplicated(colnames(sim_count_copula.all_1))
      ]
  )
  
  # 
  # (sim_count_copula.all_1[1:10,1:10])
  
  # banksy object, and make into function, and repeat and collect stats. 
  # compute ARIs. 
  locs1 = as.data.frame(domain.celltypes.dt[, c('sdimx', 'sdimy')])
  rownames(locs1)<-paste0('cell_', 1:nrow(locs1))
  # 
  # # add noise?
  # if (amount.noise>0){
  #   sim_count_copula.all_1.temp = round((1-amount.noise)*sim_count_copula.all_1 
  #                                  + amount.noise*sim_count_copula.all_1[,sample.int(ncol(sim_count_copula.all_1))])
  # }
  # 
  colnames(sim_count_copula.all_1) = rownames(locs1)
  bank.unnorm = BanksyObject(own.expr = sim_count_copula.all_1, 
                             cell.locs = locs1, 
                             meta.data = as.data.frame(domain.celltypes.dt))
  
  
  k_geom = c(15, 30)
  m.val = 1
  # lam = c(0, 0.2, 0.8)
  npcs = 20
  k_expr = 50
  res = seq(0.8, 1.2, 0.2)
  
  res.ct = seq(1, 1.5, 0.25)
  res.dom = seq(0.8, 1.2, 0.2)
  
  bank1 <- copy(bank.unnorm)
  bank1 <- NormalizeBanksy(bank1)
  bank1 <- ComputeBanksy(bank1, k_geom = k_geom, M = m.val)
  bank1 <- ScaleBanksy(bank1)
  
  bank.celltype = copy(bank1)
  bank.domain = copy(bank1)
  # do the cell typing and domain finding lambdas separately. 
  lam.celltype
  lam.domain
  
  # do cell types and domain separately
  set.seed(44*rand.seed.idx)
  bank.celltype <- Banksy:::RunBanksyPCA(bank.celltype,
                                         lambda = lam.celltype, 
                                         npcs = npcs, 
                                         M = m.val)
  bank.celltype@meta.data$clust_domain<-NULL
  bank.domain <- Banksy:::RunBanksyPCA(bank.domain,
                                       lambda = lam.domain, 
                                       npcs = npcs, 
                                       M = m.val)
  bank.domain@meta.data$clust_celltypes<-NULL
  
  
  set.seed(45*rand.seed.idx)
  bank.celltype <- Banksy:::RunBanksyUMAP(bank.celltype,
                                          lambda = lam.celltype, 
                                          npcs = npcs, 
                                          nneighbors = k_expr, 
                                          M = m.val)
  bank.domain <- Banksy:::RunBanksyUMAP(bank.domain,
                                        lambda = lam.domain, 
                                        npcs = npcs, 
                                        nneighbors = k_expr, 
                                        M = m.val)
  set.seed(46*rand.seed.idx)
  bank.celltype <- ClusterBanksy(bank.celltype,
                                 lambda = lam.celltype, 
                                 pca = TRUE, 
                                 npcs = npcs,
                                 M = m.val,
                                 method = 'leiden', 
                                 k.neighbors = k_expr, 
                                 resolution = res.ct)
  bank.domain <- ClusterBanksy(bank.domain,
                               lambda = lam.domain, 
                               pca = TRUE, 
                               npcs = npcs,
                               M = m.val,
                               method = 'leiden', 
                               k.neighbors = k_expr, 
                               resolution = res.dom)
  # saveRDS(bank1, file = paste0(results.dir, 'bank1.rds'))
  # bank1 = readRDS(file = paste0(results.dir, 'bank1.rds'))
  
  
  set.seed(47*rand.seed.idx)
  
  cnms.celltype <- clust.names(bank.celltype)
  cnms.domain <- clust.names(bank.domain)
  
  map.celltype <- 'clust_celltypes'
  map.domain <- 'clust_domain'
  
  set.seed(48*rand.seed.idx)
  bank.celltype <- ConnectClusters(bank.celltype, 
                                   map.to = map.celltype)
  bank.domain <- ConnectClusters(bank.domain, 
                                 map.to = map.domain)
  
  
  num_clusters.celltype<-max(bank.celltype@meta.data[,cnms.celltype])
  num_clusters.domain<-max(bank.domain@meta.data[,cnms.domain])
  
  # Compute ARIs
  ground.truth = bank.domain@meta.data$clust_domain
  query.clustering = bank.domain@meta.data$clust_M1_lam0.8_k50_res0.8
  ari.val.domains = unlist(lapply(clust.names(bank.domain), function(query.c){
    adjustedRandIndex(ground.truth, bank.domain@meta.data[[query.c]])
  }))
  
  ground.truth = bank.celltype@meta.data$clust_celltypes
  query.clustering = bank.celltype@meta.data$clust_M1_lam0.8_k50_res0.8 # whats this??
  ari.val.celltypes = unlist(lapply(clust.names(bank.celltype), function(query.c){
    adjustedRandIndex(ground.truth, bank.celltype@meta.data[[query.c]])
  }))
  
  # 
  # all.runs <- clust.names(bank.manual)
  # for (i in 1:length(all.runs)){
  #   subplot.nrow <- ceiling(length(unique(bank.manual@meta.data[[all.runs[i]]]))/4)
  #   print(subplot.nrow)
  #   png(paste0(results.dir, '/spatial_facet_',all.runs[i], '.png'),
  #       height = subplot.nrow*2, width = 15, res=200, units = 'in')
  #   pp<-plotSpatial(bank.manual,
  #                   by = all.runs[i],
  #                   pt.size = 1,
  #                   type = 'discrete') + facet_wrap(~feature, ncol = 4)
  #   print(pp)
  #   dev.off()
  # }
  
  
  # all.runs <- clust.names(bank.manual)
  # for (i in 1:length(all.runs)){
  #    png(paste0(results.dir, '/spatial_',all.runs[i], '.png'),
  #       height = 2, width = 4, res=200, units = 'in')
  #   pp<-plotSpatial(bank.manual,
  #                   by = all.runs[i],
  #                   pt.size = 1,
  #                   type = 'discrete')
  #   print(pp)
  #   dev.off()
  # }
  
  png(paste0(results.dir, '/spatial_celltypes_seed' , rand.seed.idx, '.png'), height = 1.2*6*scale.factor/1.5, 
      width = 1.2*12*scale.factor, res = 300, units = 'in')
  spatialplots <- mapply(FUN = function(by, ari.val) plotSpatial(bank.celltype,
                                                                 by = by, 
                                                                 main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                 type = 'discrete', 
                                                                 pt.size = 1.5*sqrt(scale.factor), 
                                                                 main.size = 18),
                         by = clust.names(bank.celltype)[2:length(clust.names(bank.celltype))], 
                         ari.val = ari.val.celltypes[2:length(ari.val.celltypes)], 
                         SIMPLIFY = FALSE)
  do.call("grid.arrange", c(spatialplots, ncol = 3))
  dev.off()
  
  all.runs = rep(clust.names(bank.celltype), 2)
  reductions <- rep(names(bank.celltype@reduction)[3:4], each = 7)
  png(paste0(results.dir, '/umap_all_celltype.png'), height = 8, 
      width = 28, res = 100, units = 'in')#, units = 'in',  res = 200)
  umapdims <- mapply(FUN = function(reduction, by, ari.val) plotReduction(bank.celltype,
                                                                          reduction = reduction,
                                                                          by = by,
                                                                          main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                          type = 'discrete',
                                                                          pt.size = 0.1, main.size = 10),
                     reduction = as.list(reductions),
                     by = as.list(all.runs),
                     ari.val = as.list(rep(round(ari.val.celltypes,3), 2)), 
                     SIMPLIFY = FALSE)
  do.call("grid.arrange", c(umapdims, ncol = 7))
  dev.off()
  
  png(paste0(results.dir, '/celltype_ground_truth_seed' , rand.seed.idx, '.png'), height = 4*scale.factor, 
      width = 2*4*scale.factor, res = 300, units = 'in')
  spatialplots <- mapply(FUN = function(by, ari.val) plotSpatial(bank.celltype,
                                                                 by = by, 
                                                                 main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                 type = 'discrete', 
                                                                 pt.size = 1.5*sqrt(scale.factor), 
                                                                 main.size = 18),
                         by = 'clust_celltypes', 
                         ari.val = ari.val.celltypes[1], 
                         SIMPLIFY = FALSE)
  do.call("grid.arrange", c(spatialplots, ncol = 1))
  dev.off()
  
  
  
  # domain
  png(paste0(results.dir, '/spatial_domain_seed' , rand.seed.idx, '.png'), height = 3*scale.factor/1.5, 
      width = 12*scale.factor, res = 300, units = 'in')
  spatialplots <- mapply(FUN = function(by, ari.val) plotSpatial(bank.domain,
                                                                 by = by, 
                                                                 main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                 type = 'discrete', 
                                                                 pt.size = 1.5*sqrt(scale.factor), 
                                                                 main.size = 18),
                         by = clust.names(bank.domain)[2:length(clust.names(bank.domain))], 
                         ari.val = ari.val.domains[2:length(ari.val.domains)], 
                         SIMPLIFY = FALSE)
  do.call("grid.arrange", c(spatialplots, ncol = 3))
  dev.off()
  
  png(paste0(results.dir, '/domain_ground_truth_seed' , rand.seed.idx, '.png'), height = 4*scale.factor, 
      width = 1.5*4*scale.factor, res = 300, units = 'in')
  spatialplots <- mapply(FUN = function(by, ari.val) plotSpatial(bank.domain,
                                                                 by = by, 
                                                                 main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                 type = 'discrete', 
                                                                 pt.size = 1.5*sqrt(scale.factor), 
                                                                 main.size = 18),
                         by = 'clust_domain', 
                         ari.val = ari.val.domains[1], 
                         SIMPLIFY = FALSE)
  do.call("grid.arrange", c(spatialplots, ncol = 1))
  dev.off()
  
  # saveRDS(bank.unnorm , file = paste0(results.dir, '/bank.unnormalized_seedidx_', rand.seed.idx, '.rds'))
  
  
  all.runs = rep(clust.names(bank.domain), 1)
  reductions <- rep(names(bank.domain@reduction)[2], each = 4)
  png(paste0(results.dir, '/umap_all_domain.png'), height = 4, 
      width = 16, res = 100, units = 'in')#, units = 'in',  res = 200)
  umapdims <- mapply(FUN = function(reduction, by, ari.val) plotReduction(bank.domain,
                                                                          reduction = reduction,
                                                                          by = by,
                                                                          main = paste0(by, ' ARI: ', round(ari.val, 3)),
                                                                          type = 'discrete',
                                                                          pt.size = 0.1, main.size = 10),
                     reduction = as.list(reductions),
                     by = as.list(all.runs),
                     ari.val = as.list(rep(round(ari.val.domains,3), 1)), 
                     SIMPLIFY = FALSE)
  do.call("grid.arrange", c(umapdims, ncol = 4))
  dev.off()
  
  return(list(bank.unnormalized = bank.unnorm, 
              bank.celltype = bank.celltype, 
              bank.domain = bank.domain, 
              ari.val.celltypes = ari.val.celltypes, 
              ari.val.domains = ari.val.domains))
}
