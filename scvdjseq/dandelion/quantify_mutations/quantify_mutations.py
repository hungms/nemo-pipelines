import sys
import pandas as pd
import ipywidgets as widgets
import os
import warnings
import hdf5plugin
import dandelion as ddl
import scanpy as sc
import anndata
os.chdir('/nemo/lab/caladod/working/Matthew/project/anqi/AX_GSE230705')

SC23355_Tcell = sc.read_h5ad("/camp/home/hungm/working/Matthew/project/anqi/LZ_SC23355/anndata/Tcell/SC23355.Tcell.h5ad")
adata = SC23355_Tcell
obs_levels = adata.obs['sample.id'].unique()
samples = []

for dirname in obs_levels:
    filename = ''.join(['/camp/home/hungm/scratch/hungm/lingling/LZ_SC23355/dandelion/tcr/', dirname, '/dandelion/all_contig_dandelion.tsv'])
    samples.append(filename)
    
adata_list = []
for level in obs_levels:
    mask = adata.obs['sample.id'] == level
    subset_adata = adata[mask].copy()  # make a copy to avoid modifying the original object
    adata_list.append(subset_adata)

vdj_list = []
for i in range(len(samples)):
    vdj = ddl.read_10x_airr(samples[i])
    ddl.tl.transfer(adata_list[i], vdj)
    vdj_list.append(vdj)
adata_vdj = anndata.concat(adata_list)
vdj = ddl.concat(vdj_list)

vdj, adata_vdj = ddl.pp.filter_contigs(vdj, adata_vdj, library_type="tr-ab")
ddl.tl.find_clones(vdj, identity = 1, key = "junction") #, identity = 1, key = "junction"
#ddl.tl.generate_network(vdj, layout_method = 'sfdp')
ddl.tl.clone_size(vdj)
ddl.tl.clone_size(vdj, max_size=4)
ddl.tl.transfer(adata_vdj, vdj)

vdj.write_h5ddl("output/vdj_analysis/dandelion/SC22272_Tcell_vdj.h5ddl")
vdj.write_airr(filename='output/vdj_analysis/dandelion/SC22272_Tcell_airr.tsv')
adata_vdj.write_h5ad(filename = "anndata/Tcell/SC23355.Tcell.h5ad")
adata_vdj.obs.to_csv("anndata/Tcell/metadata.csv")
