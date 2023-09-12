
import anndata as ad
import pandas as pd
from scipy import sparse, io

adata = ad.read_h5ad('../data/Full_obj_log_counts_soupx_v2.h5ad')

adata = adata[
    (adata.obs.Diagnosis == 'Healthy adult')
].copy()

io.mmwrite("../data/Colon-Healthy-adult.mtx", adata.X)
adata.obs.to_csv('../data/Colon-Healthy-adult-metadata.csv.gz', compression='gzip', index=False)
pd.DataFrame(adata.var_names, columns=['genes']).to_csv('../data/Colon-Healthy-adult-genes.csv.gz', compression='gzip', index=False)

