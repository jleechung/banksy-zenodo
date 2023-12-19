Download the following datasets into the data directory: 

1. Mouse Hypothalamus MERFISH data. Download the file Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv into the data directory (banksy-zenodo/fig3-hypothalamus/data).

Link to dataset: https://doi.org/10.5061/dryad.8t8s248

2. Mouse Hypothalamus scRNASeq data are available at Gene Expression Omnibus (GEO) (GSE113576).
These can be automatically downloaded using the script process_hypothalamus_scrna.R in the src directory. 

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_barcodes.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_genes.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_matrix.mtx.gz

------------------

After downloading the datasets, the scripts in the src directory can be run to regenerate the paper analyses, with the results stored in the /out directory. 

Note that the following intermediate .rds files can be found in the /data directory: 
banksyObj_naive.rds has the clustering results for the merfish analysis, and can be loaded 
for use with the process_hypothalamus_merfish.R script. 