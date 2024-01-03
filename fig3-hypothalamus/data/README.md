Download the following datasets into the fig3-hypothalamus/data directory: 

1. Mouse Hypothalamus MERFISH data. Download the file Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv into the data directory (banksy-zenodo/fig3-hypothalamus/data).

Link to dataset: https://doi.org/10.5061/dryad.8t8s248

2. Mouse Hypothalamus scRNASeq data are available at Gene Expression Omnibus (GEO) (GSE113576).
These can be automatically downloaded using the script process_hypothalamus_scrna.R in the src directory. 

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_barcodes.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_genes.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_matrix.mtx.gz

3. Semi processed data files to download (to help with exact reproducibility, 
see the respective scripts mentioned below)

-- banksyObj_provided.rds used in src/process_hypothalamus_merfish.R and others. 
https://www.dropbox.com/scl/fi/eq05c2gip0g61vc0i6n1e/banksyObj_provided.rds?rlkey=z7w2tywtn4h8jiapeymv0l54e&dl=0

-- bank_exc_rd2.rds, used in src/supp_hypothalamus_neurons.R
https://www.dropbox.com/scl/fi/bv451ijb5g7xoozy7iep0/bank_exc_rd2.rds?rlkey=2oa8inowc5omusnkbge9htmsg&dl=0

-- bank_inh_rd2.rds, used in src/supp_hypothalamus_neurons.R
https://www.dropbox.com/scl/fi/i6d0dieap1zrrfw5tl4ai/bank_inh_rd2.rds?rlkey=5plq543j1fqqszp0q6cmb058f&dl=0
 
-- bank_anim1_lamsweep.rds used in supp_hypothalamus_lambda.R
https://www.dropbox.com/scl/fi/72skj0ms4akw1s3szfghh/bank_anim1_lamsweep.rds?rlkey=9rd7cubivtaonaf7nwc1ipnsl&dl=0 
------------------

After downloading the datasets, the scripts in the src directory can be run to regenerate the paper analyses, with the results stored in the /out directory. 
