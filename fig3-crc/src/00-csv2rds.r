
#' Converts vz-ffpe-showcase csv files to rds

library(data.table)

msg <- message

# Get individual sample prefixes
all_samples <- list.files('../data', pattern = 'cell_by_gene')
all_samples <- gsub('_.*', '', all_samples)

# Process samples
process <- TRUE
if (process) {

invisible(lapply(all_samples, function (sam) 

{

msg('\n', sam, '\n')

path_gcm <- list.files('../data', pattern = paste0(sam, '_cell_by_gene'), full.names = TRUE)
path_mdata <- list.files('../data', pattern = paste0(sam, '_cell_metadata'), full.names = TRUE)

msg(path_gcm)
msg(path_mdata)

gcm <- fread(path_gcm)
mdata <- fread(path_mdata)
gcm <- gcm[order(gcm$cell),]
mdata <- mdata[order(mdata$V1),]

print(gcm[1:5,1:5])
print(mdata[1:5,1:5])

locs <- data.frame(mdata[,c('center_x', 'center_y')])
mdata <- data.frame(mdata[,c('fov', 'volume', 'V1')])
gcm <- as.matrix(t(gcm[,-1]))

colnames(locs) <- c('sdimx', 'sdimy')
colnames(gcm) <- rownames(locs) <- rownames(mdata) <- paste0('cell_', 1:nrow(locs))

print(gcm[1:5,1:5])
print(mdata[1:5,])
print(locs[1:5,])

path_save_gcm <- gsub( '_.*', '_gcm.rds', path_gcm)
path_save_locs <- gsub('_.*', '_locs.rds', path_gcm)
path_save_mdata <- gsub('_.*', '_metadata.rds', path_gcm)

message('Saving objects')
saveRDS(gcm, path_save_gcm)
saveRDS(locs, path_save_locs)
saveRDS(mdata, path_save_mdata)


}))

}

# Rm .csv files
all_samples <- list.files('../data', pattern = 'csv', full.names = TRUE)
file.remove(all_samples)
