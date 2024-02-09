#!/bin/bash

FILE_PATH="/RNA_ATAC_output/down_sampling/features/B"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for down_sample_rate in "0.2" "0.4" "0.6" "0.8" "1"
    do
        for number in "1" "2" "3" "4" "5"
        do
        {
            mkdir "$FILE_PATH/${down_sample_rate}_${number}"
            python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/down_sampling/features/run_model.py --model basic --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number --down_sample_rate $down_sample_rate
        }
        done
    done
done
FILE_PATH="/RNA_ATAC_output/down_sampling/features/C"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for down_sample_rate in "0.2" "0.4" "0.6" "0.8" "1"
    do
        for number in "1" "2" "3" "4" "5"
        do
        {
            mkdir "$FILE_PATH/${down_sample_rate}_${number}"
            python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/down_sampling/features/run_model.py --model multiVI_augmentation --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number --down_sample_rate $down_sample_rate
        }
        done
    done
done
FILE_PATH="/RNA_ATAC_output/down_sampling/features/T"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for down_sample_rate in "0.2" "0.4" "0.6" "0.8" "1"
    do
        for number in "1" "2" "3" "4" "5"
        do
        {
            mkdir "$FILE_PATH/${down_sample_rate}_${number}"
            python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/down_sampling/features/run_model.py --model celltype_augmentation --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number --down_sample_rate $down_sample_rate
        }
        done
    done
done