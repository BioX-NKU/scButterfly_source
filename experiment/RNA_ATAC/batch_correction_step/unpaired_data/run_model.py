import getopt
import sys
import gc
import os
sys.path.append('/home/atac2rna/program/atac2rna/Model/butterfly/Butterfly/')
from data_processing import RNA_data_preprocessing, ATAC_data_preprocessing
import scanpy as sc
import anndata as ad
import pandas as pd

opts, args = getopt.gnu_getopt(sys.argv[1:], 't:d:f:n', ['model_type=', 'data=', 'file_path=', 'number='])
model_type = opts[0][1]
data = opts[1][1]
file_path = opts[2][1]
number = opts[3][1]

if data == 'UP_HK':
    
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Muto-2021-ATAC.h5ad')
    RNA_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Muto-2021-RNA.h5ad')
    RNA_data.obs.index = pd.Series([str(i) for i in range(len(RNA_data.obs.index))])
    ATAC_data.obs.index = pd.Series([str(i) for i in range(len(ATAC_data.obs.index))])

elif data == 'UP_MPMC':

    RNA_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Yao-2021-RNA.h5ad')
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Yao-2021-ATAC.h5ad')
    RNA_data.obs.index = pd.Series([str(i) for i in range(len(RNA_data.obs.index))])
    ATAC_data.obs.index = pd.Series([str(i) for i in range(len(ATAC_data.obs.index))])

elif data == 'simulated data':

    RNA_data = sc.read_h5ad('/RNA_ATAC_output/batch_correction_step/simulated_data/CL_simulated_rna_3000.h5ad')
    ATAC_data = sc.read_h5ad('/RNA_ATAC_output/batch_correction_step/simulated_data/CL_simulated_atac_3000.h5ad')
    ATAC_data_temp = sc.read_h5ad('/home/atac2rna/data/atac2rna/data/scCAT/cellline/ATAC_data.h5ad')
    RNA_data.obs.index = pd.Series([str(i) for i in range(len(RNA_data.obs.index))])
    ATAC_data.obs.index = pd.Series([str(i) for i in range(len(ATAC_data.obs.index))])
    ATAC_data.var = ATAC_data_temp.var

elif data in ['UP_eye', 'UP_pancreas', 'UP_muscle', 'UP_spleen', 'UP_stomach', 'UP_thymus']:
    
    organ = data.split('_')[1]
    RNA_data = sc.read_h5ad("/data/cabins/chenshengquan/hcaRNA/hca_RNA.h5ad")
    RNA_data = RNA_data[RNA_data.obs.Organ == organ]
    gc.collect()
    RNA_data.obs['cell_type'] = list(RNA_data.obs.Main_cluster_name)
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scCAS220601/FCA_fetal_33184180_'+organ+'.h5ad')
    RNA_data.obs.index = pd.Series([str(i) for i in range(len(RNA_data.obs.index))])
    ATAC_data.obs.index = pd.Series([str(i) for i in range(len(ATAC_data.obs.index))])



############################################################
# Part 1 data processing
RNA_data = RNA_data_preprocessing(
    RNA_data,
    normalize_total=True,
    log1p=True,
    use_hvg=True,
    n_top_genes=3000,
    save_data=False,
    file_path=file_path,
    logging_path=file_path
    )
ATAC_data = ATAC_data_preprocessing(
    ATAC_data,
    binary_data=True,
    filter_features=True,
    fpeaks=0.005,
    tfidf=True,
    normalize=True,
    save_data=False,
    file_path=file_path,
    logging_path=file_path
)[0]

############################################################
# batch correction step before training scButterfly
import scanpy as sc
import scgen
import inmoose
import numpy as np
if 'batch' in RNA_data.obs.keys():
    if len(RNA_data.obs.batch.cat.categories) > 1:

        scgen.SCGEN.setup_anndata(RNA_data, batch_key="batch", labels_key="cell_type")
        model = scgen.SCGEN(RNA_data)
        model.train(
            max_epochs=100,
            batch_size=32,
            early_stopping=True,
            early_stopping_patience=25,
        )
        RNA_data.X = model.batch_removal().X
        RNA_data.X = RNA_data.X - np.min(RNA_data.X)
if 'batch' in ATAC_data.obs.keys():
    if len(ATAC_data.obs.batch.cat.categories) > 1:
        corrected_data = inmoose.pycombat.pycombat_norm(ATAC_data.X.toarray().T, list(ATAC_data.obs.batch))
        ATAC_data.X = corrected_data.T
        ATAC_data.X = (ATAC_data.X - np.min(ATAC_data.X))/(np.max(ATAC_data.X) - np.min(ATAC_data.X))
############################################################
# Part 2 split datasets
from split_datasets import *

id_list = unpaired_split_dataset(RNA_data, ATAC_data)


train_id_r, train_id_a, validation_id_r, validation_id_a, test_id_r, test_id_a = id_list[int(number) - 1]

############################################################
# Part 3 calculate chrom list
if data == 'UP_HK':

    chrom_list = []
    last_one = ''
    for i in range(len(ATAC_data.var.chrom)):
        temp = ATAC_data.var.chrom[i]
        if temp[0 : 3] == 'chr':
            if not temp == last_one:
                chrom_list.append(1)
                last_one = temp
            else:
                chrom_list[-1] += 1
        else:
            chrom_list[-1] += 1

elif data == 'UP_MPMC':
    chrom_list = []
    last_one = ''
    for i in range(len(ATAC_data.var.chrom)):
        temp = ATAC_data.var.chrom[i]
        if temp[0 : 3] == 'chr':
            if not temp == last_one:
                chrom_list.append(1)
                last_one = temp
            else:
                chrom_list[-1] += 1
        else:
            chrom_list[-1] += 1

elif data in ['UP_eye', 'UP_pancreas', 'UP_muscle', 'UP_spleen', 'UP_stomach', 'UP_thymus']:
    chrom_list = []
    last_one = ''
    k = -1
    for i in range(len(ATAC_data.var_names)):
        temp = ATAC_data.var_names[i].split('-')[0][3]
        if not temp == last_one:
            chrom_list.append(1)
            k += 1
            last_one = temp
        else:
            chrom_list[k] += 1


from train_model import Model
import torch
import torch.nn as nn

RNA_input_dim = len([i for i in RNA_data.var['highly_variable'] if i])
ATAC_input_dim = ATAC_data.X.shape[1]

R_kl_div = 1 / RNA_input_dim * 20
A_kl_div = 1 / ATAC_input_dim * 20
kl_div = R_kl_div + A_kl_div

############################################################
# Part 4 construct model
model = Model(
    R_encoder_nlayer = 2, 
    A_encoder_nlayer = 2,
    R_decoder_nlayer = 2, 
    A_decoder_nlayer = 2,
    R_encoder_dim_list = [RNA_input_dim, 256, 128],
    A_encoder_dim_list = [ATAC_input_dim, 32 * len(chrom_list), 128],
    R_decoder_dim_list = [128, 256, RNA_input_dim],
    A_decoder_dim_list = [128, 32 * len(chrom_list), ATAC_input_dim],
    R_encoder_act_list = [nn.LeakyReLU(), nn.LeakyReLU()],
    A_encoder_act_list = [nn.LeakyReLU(), nn.LeakyReLU()],
    R_decoder_act_list = [nn.LeakyReLU(), nn.LeakyReLU()],
    A_decoder_act_list = [nn.LeakyReLU(), nn.Sigmoid()],
    translator_embed_dim = 128, 
    translator_input_dim_r = 128,
    translator_input_dim_a = 128,
    translator_embed_act_list = [nn.LeakyReLU(), nn.LeakyReLU(), nn.LeakyReLU()],
    discriminator_nlayer = 1,
    discriminator_dim_list_R = [128],
    discriminator_dim_list_A = [128],
    discriminator_act_list = [nn.Sigmoid()],
    dropout_rate = 0.1,
    R_noise_rate = 0.5,
    A_noise_rate = 0.3,
    chrom_list = chrom_list,
    logging_path = file_path,
    RNA_data = RNA_data,
    ATAC_data = ATAC_data
)
############################################################
# Part 5 train model    
model.train(
    R_encoder_lr = 0.001,
    A_encoder_lr = 0.001,
    R_decoder_lr = 0.001,
    A_decoder_lr = 0.001,
    R_translator_lr = 0.001,
    A_translator_lr = 0.001,
    translator_lr = 0.001,
    discriminator_lr = 0.005,
    R2R_pretrain_epoch = 100,
    A2A_pretrain_epoch = 100,
    lock_encoder_and_decoder = False,
    translator_epoch = 200,
    patience = 50,
    batch_size = 64,
    r_loss = nn.MSELoss(size_average=True),
    a_loss = nn.BCELoss(size_average=True),
    d_loss = nn.BCELoss(size_average=True),
    loss_weight = [1, 2, 1, R_kl_div, A_kl_div, kl_div],
    train_id_r = train_id_r,
    train_id_a = train_id_a,
    validation_id_r = validation_id_r, 
    validation_id_a = validation_id_a, 
    output_path = file_path,
    seed = 19193,
    kl_mean = True,
    R_pretrain_kl_warmup = 50,
    A_pretrain_kl_warmup = 50,
    translation_kl_warmup = 50,
    load_model = None,
    logging_path = file_path
)
############################################################
# Part 6 test model 
model.test(
    test_id_r = test_id_r,
    test_id_a = test_id_a, 
    model_path = file_path,
    load_model = True,
    output_path = file_path,
    test_evaluate = True,
    test_cluster = True,
    test_figure = False,
    output_data = False,
    return_predict = False
)