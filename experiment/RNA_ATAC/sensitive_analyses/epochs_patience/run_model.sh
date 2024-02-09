#!/bin/bash

FILE_PATH="/RNA_ATAC_output/sensitive_analyses/epochs_patience"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for rna_pre_param in "100"
    do
        for atac_pre_param in "100"
        do
            for integ_param in "200"
            do
                for patience_param in "0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_pre_param}_${atac_pre_param}_${integ_param}_${patience_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/epochs_patience/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_pre_param $rna_pre_param --atac_pre_param $atac_pre_param --integ_param $integ_param --patience_param $patience_param
                    }&
                    done
                    wait
                done
            done
        done
    done
done
for dataset in "MCC"
do
    for rna_pre_param in "100"
    do
        for atac_pre_param in "100"
        do
            for integ_param in "100" "120" "140" "160" "180" "200" "220" "240" "260" "280" "300"
            do
                for patience_param in "50"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_pre_param}_${atac_pre_param}_${integ_param}_${patience_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/epochs_patience/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_pre_param $rna_pre_param --atac_pre_param $atac_pre_param --integ_param $integ_param --patience_param $patience_param
                    }&
                    done
                    wait
                done
            done
        done
    done
done
for dataset in "MCC"
do
    for rna_pre_param in "100"
    do
        for atac_pre_param in "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150"
        do
            for integ_param in "200"
            do
                for patience_param in "50"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_pre_param}_${atac_pre_param}_${integ_param}_${patience_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/epochs_patience/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_pre_param $rna_pre_param --atac_pre_param $atac_pre_param --integ_param $integ_param --patience_param $patience_param
                    }&
                    done
                    wait
                done
            done
        done
    done
done
for dataset in "MCC"
do
    for rna_pre_param in "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150"
    do
        for atac_pre_param in "100"
        do
            for integ_param in "200"
            do
                for patience_param in "50"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_pre_param}_${atac_pre_param}_${integ_param}_${patience_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/epochs_patience/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_pre_param $rna_pre_param --atac_pre_param $atac_pre_param --integ_param $integ_param --patience_param $patience_param
                    }&
                    done
                    wait
                done
            done
        done
    done
done