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

adata = sc.read_h5ad("anndata/GSE230705.h5ad")
obs_levels = adata.obs['sample.id'].unique()

samples = []
for dirname in obs_levels:
    filename = ''.join(['/camp/home/hungm/scratch/hungm/anqi/AX_GSE230705/dandelion/', dirname, '/dandelion/filtered_contig_dandelion.tsv'])
    samples.append(filename)
    
samples = [x for x in samples if "pop" in x]

################## if tigger is ran ######################
#fastas = []
#for dirname in obs_levels:
#    filename = ''.join(['/camp/home/hungm/scratch/hungm/anqi/AX_GSE230705/dandelion/', dirname, "/", dirname, "_heavy_igblast_db-pass_genotype.fasta"])
#    fastas.append(filename)
#fastas =[y for y in fastas if "_heavy_igblast_db-pass_genotype.fasta" in y]
fastas = "/camp/home/hungm/scratch/hungm/anqi/AX_GSE230705/dandelion/tigger/tigger_heavy_igblast_db-pass_genotype.fasta"
################## if tigger is ran ######################

adata_list = []
for level in obs_levels:
    mask = adata.obs['sample.id'] == level
    subset_adata = adata[mask].copy()  # make a copy to avoid modifying the original object
    adata_list.append(subset_adata)

vdj_list = []
for i in range(len(samples)):
    vdj = ddl.read_10x_airr(samples[i])
    vdj.store_germline_reference(
      corrected=fastas,
      germline="/camp/home/hungm/scratch/hungm/reference/alignment/dandelion_ref/container/database/germlines/imgt/human/vdj/",
      org="human",)
    ddl.pp.create_germlines(vdj)
    ddl.pp.quantify_mutations(vdj)
    ddl.pp.quantify_mutations(vdj, frequency = True)
    ddl.pp.quantify_mutations(vdj, split_locus=True)
    ddl.pp.quantify_mutations(vdj, split_locus=True, frequency = True)
    ddl.tl.transfer(adata_list[i], vdj)
    vdj_list.append(vdj)
adata_vdj = anndata.concat(adata_list)
vdj = ddl.concat(vdj_list)

vdj, adata_vdj = ddl.pp.filter_contigs(vdj, adata_vdj, library_type="ig")
ddl.tl.find_clones(vdj)
#ddl.tl.generate_network(vdj, layout_method = 'sfdp')
ddl.tl.clone_size(vdj)
ddl.tl.clone_size(vdj, max_size=4)
ddl.tl.transfer(adata_vdj, vdj)

vdj.write_h5ddl("data/bcr_vdj/dandelion/GSE230705_vdj.h5ddl")
vdj.write_airr(filename='data/bcr_vdj/dandelion/GSE230705_airr.tsv')
adata_vdj.write_h5ad(filename = "anndata/GSE230705.h5ad")
adata_vdj.obs.to_csv("anndata/vdj_metadata.csv")

