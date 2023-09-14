#!/bin/bash

FILE_PATH="/RNA_ATAC_output/cross_validation_by_celltype/basic"
mkdir "$FILE_PATH"
for dataset in "MB" "MCC" "MK" "PBMC"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/cross_validation_by_celltype/run_model.py --model basic --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done
FILE_PATH="/RNA_ATAC_output/cross_validation_by_celltype/celltype_aug"
mkdir "$FILE_PATH"
for dataset in "MB" "MCC" "MK" "PBMC"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/cross_validation_by_celltype/run_model.py --model celltype_augmentation --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done
FILE_PATH="/RNA_ATAC_output/cross_validation_by_celltype/multiVI_aug"
mkdir "$FILE_PATH"
for dataset in "MB" "MCC" "MK" "PBMC"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/cross_validation_by_celltype/run_model.py --model multiVI_augmentation --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done