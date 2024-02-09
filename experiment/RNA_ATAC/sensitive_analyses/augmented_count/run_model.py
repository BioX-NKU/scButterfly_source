import getopt
import sys
import gc
import os
sys.path.append('/home/atac2rna/program/atac2rna/Model/butterfly/Butterfly/')
from data_processing import RNA_data_preprocessing, ATAC_data_preprocessing
import scanpy as sc
import anndata as ad

opts, args = getopt.gnu_getopt(sys.argv[1:], 't:d:f:n:a', ['model_type=', 'data=', 'file_path=', 'number=', 'aug_count='])
model_type = opts[0][1]
data = opts[1][1]
file_path = opts[2][1]
number = opts[3][1]
aug_count = float(opts[4][1])
float_count, int_count = aug_count - math.floor(aug_count), math.floor(aug_count)

if data == 'MCC':
    
    RNA_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-RNA.h5ad')
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-ATAC.h5ad')

    
if model_type == 'multiVI_augmentation':
    from scvi_colab import install
    import pandas as pd
    install()
    import scvi
    import sys
    from calculate_cluster import *
    import scipy.sparse as sp
    from split_datasets import *

    adata = ad.AnnData(sp.hstack((RNA_data.X, ATAC_data.X)))
    adata.X = adata.X.tocsr()
    adata.obs = RNA_data.obs

    m = len(RNA_data.var.index)
    n = len(ATAC_data.var.index)
    adata.var.index = pd.Series([RNA_data.var.index[i] if i<m else ATAC_data.var.index[i-m] for i in range(m+n)], dtype='object')
    adata.var['modality'] = pd.Series(['Gene Expression' if i<m else 'Peaks' for i in range(m+n)], dtype='object').values

    adata.var_names_make_unique()

    adata_mvi = scvi.data.organize_multiome_anndatas(adata)

    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()

    sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))

    id_list = five_fold_split_dataset(RNA_data, ATAC_data, seed = 19193)

    train_id, validation_id, test_id = id_list[int(number) - 1]

    train_adata = idx2adata_multiVI(adata_mvi, train_id, validation_id, test_id)[0]

    scvi.model.MULTIVI.setup_anndata(train_adata, batch_key='modality')

    mvi = scvi.model.MULTIVI(
        train_adata,
        n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
        n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
    )

    mvi.train()

    train_adata.obsm["MultiVI_latent"] = mvi.get_latent_representation()
    leiden_adata = ad.AnnData(train_adata.obsm["MultiVI_latent"])
    sc.pp.neighbors(leiden_adata)
    sc.tl.leiden(leiden_adata, resolution=3)

    del train_adata, adata_mvi
    gc.collect()

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
# Part 2 split datasets
from split_datasets import *
id_list = five_fold_split_dataset(RNA_data, ATAC_data)
train_id, validation_id, test_id = id_list[int(number) - 1]
train_id_r = train_id.copy()
train_id_a = train_id.copy()
validation_id_r = validation_id.copy()
validation_id_a = validation_id.copy()
test_id_r = test_id.copy()
test_id_a = test_id.copy()

############################################################
# Part 2.5 amplification
import random
if model_type == 'multiVI_augmentation':
    copy_count = 3
    random.seed(19193)
    cell_type = leiden_adata.obs.leiden
    for i in range(len(cell_type.cat.categories)):
        cell_type_name = cell_type.cat.categories[i]
        idx_temp = list(cell_type[cell_type == cell_type_name].index.astype(int))
        for j in range(int_count - 1):
            random.shuffle(idx_temp)
            for each in idx_temp:
                train_id_r.append(train_id[each])
            random.shuffle(idx_temp)
            for each in idx_temp:
                train_id_a.append(train_id[each])
        random.shuffle(idx_temp)
        for each in idx_temp[0:int(float_count * len(idx_temp))]:
            train_id_r.append(train_id[each])
        random.shuffle(idx_temp)
        for each in idx_temp[0:int(float_count * len(idx_temp))]:
            train_id_a.append(train_id[each])

if model_type == 'celltype_augmentation':
    copy_count = 3
    random.seed(19193)
    ATAC_data.obs.index = [str(i) for i in range(len(ATAC_data.obs.index))]
    cell_type = ATAC_data.obs.cell_type.iloc[train_id]
    for i in range(len(cell_type.cat.categories)):
        cell_type_name = cell_type.cat.categories[i]
        idx_temp = list(cell_type[cell_type == cell_type_name].index.astype(int))
        for j in range(int_count - 1):
            random.shuffle(idx_temp)
            train_id_r.extend(idx_temp)
            random.shuffle(idx_temp)
            train_id_a.extend(idx_temp)
        random.shuffle(idx_temp)
        train_id_r.extend(idx_temp[0:int(float_count * len(idx_temp))])
        random.shuffle(idx_temp)
        train_id_a.extend(idx_temp[0:int(float_count * len(idx_temp))])

############################################################
# Part 3 calculate chrom list
if data == 'MCC':
    
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
    test_figure = True,
    output_data = False,
    return_predict = False
)