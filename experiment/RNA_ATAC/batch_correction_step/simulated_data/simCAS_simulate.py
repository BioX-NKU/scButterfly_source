import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import umap
from sklearn.decomposition import PCA
import episcanpy.api as epi
warnings.filterwarnings("ignore")
from simCAS.code.simCAS import *
import scanpy as sc
import numpy as np
import scipy
import anndata as ad
import pandas as pd

adata_dir = '/home/atac2rna/data/atac2rna/data/scCAT/cellline/ATAC_data.h5ad'
adata = sc.read(adata_dir)

batches = adata.obs.batch.unique()

adata_batch = adata[adata.obs.batch == batches[0]]
data_final=simCAS_generate(adata=adata_batch,use_noise=True,seed_noise=0,n_peak=adata.shape[1],rand_seed=2022,simu_type='cell_type',activation='exp_linear',stat_estimation='one_logser',cell_scale=3000/adata.shape[0])
data_final.obs['batch'] = batches[0]
data_final.var = adata.var
adata_simulate = data_final.copy()

for i,batch in enumerate(batches[1:]):
    print(batch)
    adata_batch = adata[adata.obs.batch == batch]
    data_final=simCAS_generate(adata=adata_batch,use_noise=True,seed_noise=i,n_peak=adata.shape[1],
                            n_cell_total=adata.shape[0],rand_seed=2022,simu_type='cell_type',
                    activation='exp_linear',stat_estimation='one_logser',cell_scale=3000/adata.shape[0])
    data_final.obs['batch'] = batch
    data_final.var = adata.var
    adata_simulate = ad.concat([adata_simulate,data_final])
    
adata_simulate.write('CL_simulated_atac_3000.h5ad')