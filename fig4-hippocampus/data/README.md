# Data download instructions. 

TLDR: If using the repo from GitHub, download the following files into 
fig4-hippocampus/data (we do not host data files on GitHub, obviously).
If you are using the repo from the Zenodo archive, then everything has already been already downloaded for you. 

-- bank.rds from
https://www.dropbox.com/scl/fi/7wwb3oge18f4dsd09t269/bank.rds?rlkey=kln41aetg3louna430r7l50ue&dl=0

-- bank_vfish_hip_0728_multi_domain_revision.rds
https://www.dropbox.com/scl/fi/zkjsznre08oywcf0aihm4/bank_vfish_hip_0728_multi_domain_revision.rds?rlkey=l8nludvrr04k2yu9ly9qytrtg&dl=0

-- ss_seurat_SSC_CA3_unprocessed.rds
https://www.dropbox.com/scl/fi/g8jtzh6smb33xxw966lcx/ss_seurat_SSC_CA3_unprocessed.rds?rlkey=mlyap7zcp8j9brx6mcaydp0s5&dl=0

-- ss_seurat_oligo_unprocessed.rds
https://www.dropbox.com/scl/fi/amu2eb07brpgqrmige942/ss_seurat_oligo_unprocessed.rds?rlkey=l6hq40ybktbgjozej0k9puphn&dl=0

-- Seurat.ss.rda
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/Seurat.ss.rda


------------------------------------------------------------------------------------------------------------

More detailed instructions: 

## VeraFISH data: 
This data is provided as part of the BANKSY package itself. Simply run the script 
fig4-hippocampus/src/process_hippocampus_verafish.R. 

If you do want the data, you can load the BANKSY package, and run the following commands: 

data(hippocampus)
expr <- hippocampus$expression
locs <- hippocampus$locations



## scRNA-seq data instructions 
(also found at the beginning of the file fig4-hippocampus/src/process_hippocampus_scrna.R)

Starting with the raw data from the literature involves an initial subsetting step which requires 32 GB RAM. If this is an issue, we recommend using our preprocessed Seurat objects, see below. 

For the full workflow, download the Seurat object from: 
https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq
and place it in the fig4-hippocampus/data directory. 

The relevant object can be found in the following line in the table: 

"Gene expression matrix (Seurat) |	5.4 GB | .ss.rda |	This Seurat object contains the cell-by-gene expression matrix, with introns and exons combined."

The specific link is: https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/Seurat.ss.rda

Uncomment the lines below between "START PREPROCESSING BLOCK" and "END PREPROCESSING BLOCK" to generate the individual seurat objects (ss.seurat.SSC.CA3_unprocessed for neurons, ss.seurat.oligo_unprocessed for oligodendrocytes.)

If you do not wish to subset this seurat object yourself, you can download the subsetted objects from the dropbox links below, and place them in fig4-hippocampus/data/ 
https://www.dropbox.com/scl/fi/g8jtzh6smb33xxw966lcx/ss_seurat_SSC_CA3_unprocessed.rds?rlkey=mlyap7zcp8j9brx6mcaydp0s5&dl=0
and 
https://www.dropbox.com/scl/fi/amu2eb07brpgqrmige942/ss_seurat_oligo_unprocessed.rds?rlkey=l6hq40ybktbgjozej0k9puphn&dl=0

