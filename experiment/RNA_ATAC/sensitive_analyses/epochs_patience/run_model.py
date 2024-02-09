import getopt
import sys
import gc
import os
sys.path.append('/home/atac2rna/program/atac2rna/Model/butterfly/Butterfly/')
from data_processing import RNA_data_preprocessing, ATAC_data_preprocessing
import scanpy as sc
import anndata as ad

opts, args = getopt.gnu_getopt(sys.argv[1:], 't:d:f:n:r:a:i:p', ['model_type=', 'data=', 'file_path=', 'number=', 'rna_pre_param=', 'atac_pre_param=', 'integ_param=', 'patience_param='])
model_type = opts[0][1]
data = opts[1][1]
file_path = opts[2][1]
number = opts[3][1]
rna_pre_param = int(opts[4][1])
atac_pre_param = int(opts[5][1])
integ_param = int(opts[6][1])
patience_param = int(opts[7][1])

if data == 'MCC':
    
    RNA_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-RNA.h5ad')
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-ATAC.h5ad')


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
# Part 3 calculate chrom list
    
elif data == 'MCC':
    
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
    R2R_pretrain_epoch = rna_pre_param,
    A2A_pretrain_epoch = atac_pre_param,
    lock_encoder_and_decoder = False,
    translator_epoch = integ_param,
    patience = patience_param,
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