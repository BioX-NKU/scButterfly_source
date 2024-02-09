#!/bin/bash

FILE_PATH="/RNA_ATAC_output/sensitive_analyses/augmented_count/C/"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for aug_count in "1.5" "2" "2.5" "3" "3.5" "4" "4.5"
    do
        for number in "1" "2" "3" "4" "5"
        do
        {
            mkdir "$FILE_PATH/${aug_count}_${number}"
            python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/augmented_count/run_model.py --model multiVI_augmentation --data $dataset --file "$FILE_PATH/${aug_count}_${number}" --number $number --aug_count $aug_count
        }
        done
    done
done
FILE_PATH="/RNA_ATAC_output/sensitive_analyses/augmented_count/T/"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for aug_count in "1.5" "2" "2.5" "3" "3.5" "4" "4.5"
    do
        for number in "1" "2" "3" "4" "5"
        do
        {
            mkdir "$FILE_PATH/${aug_count}_${number}"
            python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/augmented_count/run_model.py --model celltype_augmentation --data $dataset --file "$FILE_PATH/${aug_count}_${number}" --number $number --aug_count $aug_count
        }
        done
    done
done