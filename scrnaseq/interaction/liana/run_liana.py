import sys
import pandas as pd
import os
import hdf5plugin
import anndata
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
exec(open('/camp/home/hungm/working/Matthew/library/pythonlib/liana/liana_functions.py').read())

### modify here ###
os.chdir('/nemo/lab/caladod/working/Matthew/project/lingling/LZ_SC23355')
SC23355_Tcell = sc.read_h5ad("/camp/home/hungm/working/Matthew/project/lingling/LZ_SC23355/anndata/Tcell/LZ_SC23355.Tcell.h5ad")
SC23355_Bcell = sc.read_h5ad("/camp/home/hungm/working/Matthew/project/lingling/LZ_SC23355/anndata/Bcell/LZ_SC23355.Bcell.h5ad")
adatas = [SC23355_Tcell, SC23355_Bcell]
adata = anndata.concat(adatas)
del(SC23355_Tcell, SC23355_Bcell, adatas)

sample_key = 'batch'
condition_key = 'genot'
groupby = 'celltype_v2'
resource_name = 'mouseconsensus' ## selection = ['baccin2019','cellcall', 'cellchatdb', 'cellinker', 'cellphonedb', 'celltalkdb', 'connectomedb2020', 'consensus', 'embrace', 'guide2pharma', 'hpmr', 'icellnet', 'italk', 'kirouac2010', 'lrdb', 'mouseconsensus', 'ramilowski2015']
save = '/data/interaction/liana/SC23355_tensor.pkl'
###################

# run INTRACELLULAR CONTEXT FACTORIZATION (TENSOR-C2C) 
liana_tensor_c2c(adata, sample_key = sample_key, condition_key = condition_key, groupby = groupby, resource_name = resource_name, save = save, n_job = 20)

