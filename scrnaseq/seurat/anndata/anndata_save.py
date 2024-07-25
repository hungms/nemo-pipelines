# change directory and final file name

import sys
import pandas as pd
import ipywidgets as widgets
import os
import hdf5plugin
import anndata as ad
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
os.chdir('/nemo/lab/caladod/working/Matthew/project/oscar/OA_ZND6629388/')

meta = pd.read_csv('anndata/metadata.csv',index_col=0)
meta = meta.iloc[:, :-4]
counts = pd.read_csv("anndata/rna_counts.csv",index_col=0)
counts.columns = [col.replace('.', '-') for col in counts.columns]
adata = sc.AnnData(X=counts.transpose(), obs=meta)

print(adata)
print(adata.X)
print(adata.obs.head())
print(adata.obs_names)

#lognormalized = pd.read_csv('anndata/rna_data.csv',index_col=0)
#lognormalized = lognormalized.transpose()
#lognormalized = csr_matrix(lognormalized.values, dtype=np.float32)
#adata.layers["lognormalized"] = lognormalized
#adata.X = adata.to_df(layer="lognormalized")
#adata.raw = adata

X_pca = pd.read_csv("anndata/pca.csv",index_col=0)
adata.obsm["X_pca"] = X_pca.values

X_umap = pd.read_csv("anndata/reduction.csv",index_col=0)
adata.obsm["X_umap"] = X_umap.values

#snn_df = pd.read_csv('anndata/snn.csv',index_col=0)
#nn_df = pd.read_csv('anndata/nn.csv',index_col=0)
#snn_df = snn_df.T
#nn_df = nn_df.T

#snn_df.index = snn_df.index.str.replace('.', '-')
#nn_df.index = nn_df.index.str.replace('.', '-')
#snn_df.index = snn_df.index.str.replace('X', '')
#nn_df.index = nn_df.index.str.replace('X', '')

#adata.obsp['connectivities'] = snn_df.to_numpy()
#adata.obsp['distances'] = nn_df.to_numpy()

#adata.uns['neighbors'] = {}
#adata.uns['neighbors']['distances'] = adata.obsp['distances']
#adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

adata.write_h5ad("anndata/ZND6629388_B.h5ad") #, compression=hdf5plugin.FILTERS["zstd"]

