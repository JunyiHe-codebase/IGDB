# /usr/bin python3

# -*- coding: utf-8 -*-
# @time    : 2026/01/03 13:43
# @author  : HeJunyi
# @file    : split_ref_genome.py

import argparse

def split_reference(ref_genome: str, split_num: int):
    genome_dict = {}
    ref_name = ref_genome.replace("fasta", "").replace("fa", "")
    with open(ref_genome, "r") as f:
        tmp_genome_dict = {}
        for line in f:
            line = line.strip("\n")
            if line.startswith(">"):
                chr = line[1:]
                tmp_genome_dict[chr] = []
            else:
                seq = line
                tmp_genome_dict[chr].append(seq)

    for k, v in tmp_genome_dict.items():
        genome_dict[k] = "".join(v)

    if len(genome_dict) <= split_num:
        for k, v in genome_dict.items():
            with open(ref_name + f"{k}.fasta", "w") as f:
                f.write(f">{k}\n")
                seq_lines = [v[i:i + 80] for i in range(0, len(v), 80)]
                f.write("\n".join(seq_lines))
    else:
        count = 1
        others_genome_dict = {}
        for k, v in genome_dict.items():
            if count <= split_num - 1:
                with open(ref_name + f"{k}.fasta", "w") as f:
                    f.write(f">{k}\n")
                    seq_lines = [v[i:i + 80] for i in range(0, len(v), 80)]
                    f.write("\n".join(seq_lines))
                count += 1
            else:
                others_genome_dict[k] = v
        with open(ref_name + f"others.fasta", "w") as f:
            for k, v in others_genome_dict.items():
                f.write(f">{k}\n")
                seq_lines = [v[i:i + 80] for i in range(0, len(v), 80)]
                f.write("\n".join(seq_lines))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref_genome', type=str, required=True,
                        help='Input reference genome')
    parser.add_argument('-s', '--split_num', type=int, required=True,
                        help='Input the number of split reference genome')

    args = parser.parse_args()

    split_reference(args.ref_genome, args.split_num)
    