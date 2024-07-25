iteration_path = "/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/data/functional/pyscenic/iterations/"
reference_path = "/camp/home/hungm/scratch/hungm/reference/gene_annotations/pyscenic_tfs/mus_musculus/"

################################################
import os
import numpy as np
import pandas as pd
import scanpy as sc
import hdf5plugin
import loompy as lp
import ipywidgets as widgets
import matplotlib.pyplot as plt
import seaborn as sns
import glob

f_loom = str(iteration_path) + "input/*.loom"
f_pyscenic = [str(iteration_path) + "output/pyscenic_output_" + str(s) + '.loom' for s in range(1,21)]
adj = [str(iteration_path) + "adj/adj_" + str(s) + ".csv" for s in range(1,21)]
reg = [str(iteration_path) + "reg/reg_" + str(s) + ".csv" for s in range(1,21)]

f_tfs = str(reference_path) + "allTFs_mm.txt"
f_db_glob = str(reference_path) + "*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )
f_motif = str(reference_path) + "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

import subprocess
command1 = 'pyscenic grn {f_loom} {f_tfs} -o {f_adj} --num_workers 32'
command2 = 'pyscenic ctx {f_adj} {f_db_names} --annotations_fname {f_motif} --expression_mtx_fname {f_loom} --output {f_reg} --mask_dropouts --num_workers 32'
command3 = 'pyscenic aucell {f_loom} {f_reg} --output {f_pyscenic} --num_workers 32'

for i in range(0, 20):
        subprocess.run(command1.format(f_loom = f_loom, f_tfs = f_tfs, f_adj = adj[i]), shell=True, stdout=subprocess.PIPE)
        subprocess.run(command2.format(f_adj = adj[i], f_db_names = f_db_names, f_motif = f_motif, f_loom = f_loom, f_reg = reg[i]), shell=True, stdout=subprocess.PIPE)
        subprocess.run(command3.format(f_loom = f_loom, f_reg = reg[i], f_pyscenic = f_pyscenic[i]), shell=True, stdout=subprocess.PIPE)
