#!/bin/bash

FILE_PATH="/RNA_ATAC_output/sensitive_analyses/loss_weights"
mkdir "$FILE_PATH"
for dataset in "MCC"
do
    for rna_param in "1"
    do
        for atac_param in "2"
        do
            for dis_param in "1"
            do
                for kl_param in "15" "16" "17" "18" "19" "21" "22" "23" "24" "25"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/loss_weights/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_param $rna_param --atac_param $atac_param --kl_param $kl_param --dis_param $dis_param
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
    for rna_param in "1"
    do
        for atac_param in "2"
        do
            for dis_param in "0.5" "0.6" "0.7" "0.8" "0.9" "1.1" "1.2" "1.3" "1.4" "1.5"
            do
                for kl_param in "20"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/loss_weights/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_param $rna_param --atac_param $atac_param --kl_param $kl_param --dis_param $dis_param
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
    for rna_param in "1"
    do
        for atac_param in "1.5" "1.6" "1.7" "1.8" "1.9" "2.1" "2.2" "2.3" "2.4" "2.5"
        do
            for dis_param in "1"
            do
                for kl_param in "20"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/loss_weights/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_param $rna_param --atac_param $atac_param --kl_param $kl_param --dis_param $dis_param
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
    for rna_param in "0.5" "0.6" "0.7" "0.8" "0.9" "1.1" "1.2" "1.3" "1.4" "1.5"
    do
        for atac_param in "2"
        do
            for dis_param in "1"
            do
                for kl_param in "20"
                do
                    for number in "1" "2" "3" "4" "5"
                    do
                    {
                        mkdir "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}"
                        python /home/atac2rna/program/atac2rna/Model/butterfly/experiment/RNA_ATAC/sensitive_analyses/loss_weights/run_model.py --model basic --data $dataset --file "$FILE_PATH/${rna_param}_${atac_param}_${kl_param}_${dis_param}_${number}" --number $number --rna_param $rna_param --atac_param $atac_param --kl_param $kl_param --dis_param $dis_param
                    }&
                    done
                    wait
                done
            done
        done
    done
done