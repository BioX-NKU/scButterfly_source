#!/bin/bash

FILE_PATH="/RNA_ATAC_output/unpaired_data/basic"
mkdir "$FILE_PATH"
for dataset in "UP_HK" "UP_MPMC" "UP_eye" "UP_pancreas" "UP_muscle" "UP_spleen" "UP_stomach" "UP_thymus"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3" "4" "5"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/unpaired_data/run_model.py --model basic --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done

