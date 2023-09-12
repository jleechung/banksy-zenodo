
rm(list=ls())
graphics.off()

run_download <- FALSE
all_samples <- as.character(c(151507:151510, 151669:151676))
all_domains <- c(rep(7, 4), rep(5, 4), rep(7, 4))
sample_size <- length(all_samples)

# Fetch data -------------------------------------------------------------------

if (run_download) {
    
library(BiocFileCache)
library(assertthat)
library(RCurl)
library(SummarizedExperiment)

if (!(file.exists('../data'))) dir.create('../data', recursive = TRUE)

getRDS <- function(dataset, sample, cache=TRUE) {
    
    url <- "https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SpatialTranscriptomes/%s/%s.rds"
    url <- sprintf(url, dataset, sample)
    # assert_that(url.exists(url), message="Dataset/sample not available")
    
    if (cache) {
        bfc <- BiocFileCache()
        local.path <- bfcrpath(bfc, url)
    } else {
        local.path <- tempfile(fileext=".rds")
        download.file(url, local.path, quiet=TRUE, mode="wb")
    }
    
    readRDS(local.path)
}

for (curr_sample in all_samples) {
    
    message('Loading ', curr_sample)
    d <- getRDS('2020_maynard_prefrontal-cortex', curr_sample)
    
    # Manual annotation
    message('Fetch manual annotation')
    manual <- as.character(d$layer_guess_reordered_short)
    
    # Expression
    message('Get counts')
    gcm <- as.matrix(assays(d)$counts)
    
    # Locations
    message('Get locations')
    locs <- data.frame(sdimx = d$imagecol, sdimy = max(d$imagerow) - d$imagerow)
    rownames(locs) <- colnames(gcm)
    
    # Save
    message('Saving')
    saveRDS(gcm, paste0('../data/expression-', curr_sample, '.rds'))
    saveRDS(locs, paste0('../data/locations-', curr_sample, '.rds'))
    saveRDS(manual, paste0('../data/manual-', curr_sample, '.rds'))
    
}

}