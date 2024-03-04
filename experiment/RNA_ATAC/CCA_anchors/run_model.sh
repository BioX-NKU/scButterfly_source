#!/bin/bash

FILE_PATH="/RNA_ATAC_output/CCA_anchors/TA"
mkdir "$FILE_PATH"
for dataset in "MCC"
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3" "4" "5"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/CCA_anchors/run_model.py --model scButterfly_TA --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done
FILE_PATH="/RNA_ATAC_output/CCA_anchors/A"
mkdir "$FILE_PATH"
for dataset in "MCC"
    mkdir "$FILE_PATH/$dataset"
    for number in "1" "2" "3" "4" "5"
    do
    {
        mkdir "$FILE_PATH/$dataset/$number"
        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/CCA_anchors/run_model.py --model scButterfly_A --data $dataset --file "$FILE_PATH/$dataset/$number" --number $number
    }
    done
done