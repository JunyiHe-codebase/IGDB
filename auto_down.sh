#!/bin/bash

if [ $# -eq 0 ]; then
    echo "用法: $0 <accession_list_file>"
    exit 1
fi

input_file="$1"

while IFS= read -r i; do
    [ -z "$i" ] && continue
    prefetch "$i" -O . &&
    fasterq-dump "$i/$i.sra" --split-3 -e 12 -O . &&
    pigz -p 12 "${i}_1.fastq" &&
    pigz -p 12 "${i}_2.fastq"
done < "$input_file"
