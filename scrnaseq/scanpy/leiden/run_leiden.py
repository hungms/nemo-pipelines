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
os.chdir('/nemo/lab/caladod/working/Matthew/project/anqi/AX_GSE230705')

adata = sc.read_h5ad('anndata/GSE230705.h5ad')
for x in range(1,21,1):
  no = x/10
  key = 'leiden_' + str(no)
  sc.tl.leiden(adata, resolution=no, key_added=key)
  
adata.obs.to_csv("anndata/metadata.csv")
adata.write_h5ad(filename = "anndata/GSE230705.h5ad", compression=hdf5plugin.FILTERS["zstd"])
