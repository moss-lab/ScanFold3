#!/bin/bash

conda activate ScanFold3

sf3_root_dir=$(dirname "$(dirname "$PWD")")
sf3_py="$sf3_root_dir/ScanFold.py"
start_dir="$PWD"

#run across a range of window sizes
winsizes=("60" "80" "100" "120" "140" "150" "160" "180" "200" "250" "300")
for win in "${winsizes[@]}"; do
    echo "$win"
    for fasta in "$PWD/S2"/*; do
        echo "$fasta"
        fasta_nopath="${fasta##*/}"
        echo "$fasta_nopath"
        odirname="${fasta_nopath%.*}"
        echo "$odirname"
        ocwd="$PWD/output/winsizes/$win"
        echo "$ocwd"
        cd "$ocwd"
        python "$sf3_py" "$fasta" --folder_name "$odirname" --window "$win"
        cd "$start_dir"
    done
done
#run with global refolding
for fasta in "$PWD/S2"/*; do
    echo "$fasta"
    fasta_nopath="${fasta##*/}"
    echo "$fasta_nopath"
    odirname="${fasta_nopath%.*}"
    echo "$odirname"
    ocwd="$PWD/output/refold"
    echo "$ocwd"
    cd "$ocwd"
    python "$sf3_py" "$fasta" --folder_name "$odir" --benchmark --global_refold
    cd "$start_dir"
done

