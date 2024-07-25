import os
import numpy as np
import pandas as pd
import scanpy as sc
import hdf5plugin
import loompy as lp
os.chdir("/camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/anndata")

f_tfs = "/camp/home/hungm/working/Matthew/sctools/cisTarget_databases/allTFs_mm.txt"
f_loom = "/camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/anndata/SC22272_PC.loom"

import subprocess
command = 'pyscenic grn {f_loom} {f_tfs} -o adj.csv --num_workers 20'
subprocess.run(command.format(f_loom = f_loom, f_tfs = f_tfs), shell=True, stdout=subprocess.PIPE)
