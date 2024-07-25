import os
os.chdir('/nemo/lab/caladod/working/Matthew/project/anqi/20230719_SC22272')

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as scanpy
import scipy
from scipy.sparse import csr_matrix
import hdf5plugin
try:
    import scvelo as scv
except Exception:
    pass
try:
    import scvelo as scv
except Exception:
    pass
try:
    import scvelo as scv
except Exception:
    pass

adata = scanpy.read_h5ad("anndata/SC22272_PC_adata.h5ad")
ldata = scanpy.read_h5ad("anndata/SC22272_PC_ldata.h5ad")
adata = scv.utils.merge(adata, ldata)
adata.obs["SCT_rpca2_0.5"] = pd.Categorical(adata.obs["SCT_rpca2_0.5"])
scv.tl.recover_dynamics(adata, n_jobs = 30)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write_h5ad(filename = "/camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/anndata/SC22272_PC_scvelo.h5ad", compression=hdf5plugin.FILTERS["zstd"])
adata = scv.read('/camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/anndata/SC22272_PC_scvelo.h5ad')

