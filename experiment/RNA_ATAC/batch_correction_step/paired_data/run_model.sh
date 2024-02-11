#!/bin/bash

FILE_PATH="/RNA_ATAC_output/batch_correction_step/paired_data"
mkdir "$FILE_PATH"
for dataset in "BMMC" "CL" "MDS"
do
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3" "4"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/batch_correction_step/paired_data/run_model.py --model basic --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done
