#!/bin/bash

set -e


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"


echo "Method           Language     Runtime"
echo "--------------------------------------"


methods=("global" "local" "fitting" "affine")
languages=("python" "codon")


for method in "${methods[@]}"; do
    for lang in "${languages[@]}"; do
        # test q1-q5 vs t1-t5
        for i in {1..5}; do
            if [ "$lang" == "python" ]; then
                python3 code/align.py "${method}-q${i}" "$lang" "data/q${i}.fa" "data/t${i}.fa"
            else
                codon run code/align.py "${method}-q${i}" "$lang" "data/q${i}.fa" "data/t${i}.fa"
            fi
        done
        
        # test MT_human vs MT_orang
        if [ "$lang" == "python" ]; then
            :
        else
            codon run code/align.py "${method}-mt_human" "$lang" "data/MT-human.fa" "data/MT-orang.fa"
        fi
    done
done