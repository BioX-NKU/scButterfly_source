import getopt
import sys
import gc
import os
sys.path.append('/home/atac2rna/program/atac2rna/Model/butterfly/Butterfly/')
from data_processing import RNA_data_preprocessing, ATAC_data_preprocessing
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import random

opts, args = getopt.gnu_getopt(sys.argv[1:], 't:d:f:n', ['model_type=', 'data=', 'file_path=', 'number='])
model_type = opts[0][1]
data = opts[1][1]
file_path = opts[2][1]
number = opts[3][1]

    
if data == 'MCC':
    
    RNA_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-RNA.h5ad')
    ATAC_data = sc.read_h5ad('/data/cabins/chenshengquan/scglue/Chen-2019-ATAC.h5ad')


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
# Anchor
if model_type == 'scButterfly_A':
    anchors = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/anchors_fold'+number+'.scv', index_col=0)
    embeddings = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/embeddings_fold'+number+'.scv', index_col=0)
    query = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/query_fold'+number+'.scv', index_col=0)
    reference = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/reference_fold'+number+'.scv', index_col=0)

    RNA_data_train = RNA_data[train_id]
    ATAC_data_train = ATAC_data[train_id]

    RNA_index_temp = list(RNA_data_train.obs_names)
    RNA_index_temp = np.array([RNA_index_temp.index(each[:-10]) for each in reference.x])
    ATAC_index_temp = list(ATAC_data_train.obs_names)
    ATAC_index_temp = np.array([ATAC_index_temp.index(each[:-6]) for each in query.x])
    RNA_index_temp_inverse = list(RNA_index_temp)
    ATAC_index_temp_inverse = list(ATAC_index_temp)
    RNA_index_temp_inverse = np.array([RNA_index_temp_inverse.index(each) for each in range(len(train_id))])
    ATAC_index_temp_inverse = np.array([ATAC_index_temp_inverse.index(each) for each in range(len(train_id))])

    RNA_embeddings = np.array(embeddings.iloc[0:len(train_id), :])
    ATAC_embeddings = np.array(embeddings.iloc[len(train_id):, :])

    RNA_anchors = np.array(anchors.cell1)
    ATAC_anchors = np.array(anchors.cell2)

    RNA_dist = np.zeros((len(train_id), len(RNA_anchors)))
    ATAC_dist = np.zeros((len(train_id), len(ATAC_anchors)))

    for i in range(len(train_id)):
        if i % 100 == 0:
            print(i / 5882)
        for j in range(len(RNA_anchors)):
            target = RNA_embeddings[i, :]
            anchor = RNA_embeddings[RNA_anchors[j]-1, :]
            RNA_dist[i, j] = np.linalg.norm(target-anchor)
    for i in range(len(train_id)):
        if i % 100 == 0:
            print(i / 5882)
        for j in range(len(ATAC_anchors)):
            target = ATAC_embeddings[i, :]
            anchor = ATAC_embeddings[ATAC_anchors[j]-1, :]
            ATAC_dist[i, j] = np.linalg.norm(target-anchor)

    train_id_temp = np.array(train_id)        
    ATAC_list = list(train_id_temp[ATAC_index_temp[ATAC_anchors[pd.DataFrame(RNA_dist[RNA_index_temp_inverse]).T.idxmin()]-1]])
    RNA_list = list(train_id_temp[RNA_index_temp[RNA_anchors[pd.DataFrame(ATAC_dist[ATAC_index_temp_inverse]).T.idxmin()]-1]])

    train_id_r.extend(train_id)
    train_id_a.extend(ATAC_list)
    train_id_r.extend(RNA_list)
    train_id_a.extend(train_id)

############################################################
# Type Anchor
if model_type == 'scButterfly_TA':
    RNA_data_obs_names = list(RNA_data.obs_names)
    cell_type_list = RNA_data.obs.cell_type.cat.categories
    random.seed(19193)
    for ctype in cell_type_list:
        # find the cell ids for each cell type
        train_id_ctype = RNA_data[train_id][RNA_data[train_id].obs.cell_type == ctype].obs_names
        train_id_ctype = [RNA_data_obs_names.index(each) for each in train_id_ctype]
        try:
            print(ctype)
            anchors = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/'+ctype+'_anchors_fold'+number+'.csv', index_col=0)
            embeddings = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/'+ctype+'_embeddings_fold'+number+'.csv', index_col=0)
            query = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/'+ctype+'_query_fold'+number+'.csv', index_col=0)
            reference = pd.read_csv('/home/atac2rna/program/atac2rna/jqfile/data/'+ctype+'_reference_fold'+number+'.csv', index_col=0)

            # filter anchors
            anchors_filtered = []
            for i in range(len(train_id_ctype)):
                anchors_set = anchors[anchors.cell1 == i+1]
                score_for_each_anchor = anchors_set.score
                if len(score_for_each_anchor) > 0:
                    max_score_anchor = np.argmax(score_for_each_anchor)
                    anchors_filtered.append(list(anchors_set.iloc[max_score_anchor, :]))
            anchors = pd.DataFrame(anchors_filtered, columns=anchors.columns)
            anchors.cell1 = anchors.cell1.astype('int64')
            anchors.cell2 = anchors.cell2.astype('int64')

            RNA_data_train = RNA_data[train_id_ctype]
            ATAC_data_train = ATAC_data[train_id_ctype]

            RNA_index_temp = list(RNA_data_train.obs_names)
            RNA_index_temp = np.array([RNA_index_temp.index(each[:-10]) for each in reference.x])
            ATAC_index_temp = list(ATAC_data_train.obs_names)
            ATAC_index_temp = np.array([ATAC_index_temp.index(each[:-6]) for each in query.x])
            RNA_index_temp_inverse = list(RNA_index_temp)
            ATAC_index_temp_inverse = list(ATAC_index_temp)
            RNA_index_temp_inverse = np.array([RNA_index_temp_inverse.index(each) for each in range(len(train_id_ctype))])
            ATAC_index_temp_inverse = np.array([ATAC_index_temp_inverse.index(each) for each in range(len(train_id_ctype))])

            RNA_embeddings = np.array(embeddings.iloc[0:len(train_id_ctype), :])
            ATAC_embeddings = np.array(embeddings.iloc[len(train_id_ctype):, :])

            RNA_anchors = np.array(anchors.cell1)
            ATAC_anchors = np.array(anchors.cell2)

            RNA_dist = np.zeros((len(train_id_ctype), len(RNA_anchors)))
            ATAC_dist = np.zeros((len(train_id_ctype), len(ATAC_anchors)))

            train_id_temp = np.array(train_id_ctype)

            RNA_anchors_list = list(RNA_anchors)
            ATAC_anchors_list = list(ATAC_anchors)

            best_anchors_RNA = []
            for i in range(len(train_id_ctype)):
                if i+1 in RNA_anchors:
                    best_anchors_RNA.append(RNA_anchors_list.index(i+1))
                else:
                    temp = []
                    target = RNA_embeddings[i, :]
                    for j in range(len(RNA_anchors)):
                        anchor = RNA_embeddings[RNA_anchors[j]-1, :]
                        temp.append(np.linalg.norm(target-anchor))
                    best_anchors_RNA.append(np.argmax(temp))

            ATAC_list = list(train_id_temp[ATAC_index_temp[ATAC_anchors[best_anchors_RNA]-1]])

            best_anchors_ATAC = []
            for i in range(len(train_id_ctype)):
                if i+1 in ATAC_anchors:
                    best_anchors_ATAC.append(ATAC_anchors_list.index(i+1))
                else:
                    temp = []
                    target = ATAC_embeddings[i, :]
                    for j in range(len(ATAC_anchors)):
                        anchor = ATAC_embeddings[ATAC_anchors[j]-1, :]
                        temp.append(np.linalg.norm(target-anchor))
                    best_anchors_ATAC.append(np.argmax(temp))

            RNA_list = list(train_id_temp[RNA_index_temp[RNA_anchors[best_anchors_ATAC]-1]])

            train_id_r.extend(train_id_ctype)
            train_id_a.extend(ATAC_list)
            train_id_r.extend(RNA_list)
            train_id_a.extend(train_id_ctype)

        except:
            for j in range(2):
                idx_temp = train_id_ctype.copy()
                random.shuffle(idx_temp)
                train_id_r.extend(idx_temp)
                random.shuffle(idx_temp)
                train_id_a.extend(idx_temp)


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