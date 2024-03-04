#!/bin/bash

FILE_PATH="/RNA_ATAC_output/batch_correction_step/unpaired_data"
mkdir "$FILE_PATH"
for dataset in "UP_HK" "UP_MPMC" "simulated data" "UP_eye" "UP_pancreas" "UP_muscle" "UP_spleen" "UP_stomach" "UP_thymus"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3" "4" "5"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/batch_correction_step/unpaired_data/run_model.py --model basic --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done

