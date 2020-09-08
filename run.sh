#!/usr/bin/env bash

i=0

# don't forget to run from ubuntu (type bash in powershell) or WSL
for file in ./data/NO3_CO_vdw5/*.cif; do
    
    # 2>&1 makes sure it prints both stderr and stdout to file
    python3 main.py "$file" "results/NO3_CO_vdw5_9-8-20.csv" >> log.txt 2>&1
    echo i: $i "$file"
    ((i=i+1))
done