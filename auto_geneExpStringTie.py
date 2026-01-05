# /usr/bin python3

# -*- coding: utf-8 -*-
# @time    : 2025/12/25 09:07
# @author  : HeJunyi
# @file    : auto_geneExpStringTie.py

import subprocess
import datetime
import time
import concurrent.futures
import os
import argparse


    # 1) hisat2 match
def run_hisat(fastq_1: str, fastq_2: str, sample_name: str, fastq_pathway: str, index: str, threads_cmd: int):
    in_fastq_1 = os.path.join(fastq_pathway, fastq_1)
    in_fastq_2 = os.path.join(fastq_pathway, fastq_2)
    out_dir_sample = os.path.join(fastq_pathway, sample_name)

    if not os.path.exists(in_fastq_1) or not os.path.exists(in_fastq_2):
        print(f"错误：{sample_name} 的 fastq 文件未找到，路径：{in_fastq_1} 或 {in_fastq_2}")
        return 1, sample_name

    
    if os.path.exists(f"{out_dir_sample}.sam"):
        now = datetime.datetime.now()
        print(f"{sample_name} - {now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 已存在 hisat2 结果，跳过此步骤...")
        return 0, sample_name
    else:
        start = time.time()
        now = datetime.datetime.now()
        print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 开始运行 hisat2...")
        cmd = (f'hisat2 -p {threads_cmd} -q '
            f'-x {index} '
            f'-1 {in_fastq_1} '
            f'-2 {in_fastq_2} '
            f'-S {out_dir_sample}.sam 2> {out_dir_sample}.hisat2.mappingstat'
            )
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(cmd)
        now = datetime.datetime.now()
        if result.returncode != 0:
            print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - {sample_name} 运行 hisat2 失败，详见：{out_dir_sample}.hisat2.mappingstat")
            return 1, sample_name
        else:
            end = time.time()
            duration = end - start
            print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - hisat2 完成，用时{duration:.1f}s...")
            return 0, sample_name


    # 2) sort and MQ30
def run_sort_mq30(sample_name: str, fastq_pathway: str, threads_cmd: int):
    out_dir_sample = os.path.join(fastq_pathway, sample_name)
    if os.path.exists(f"{out_dir_sample}_sortMq30.bam"):
        now = datetime.datetime.now()
        print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 已存在 sort and MQ30 结果，跳过此步骤...")
        return 0, sample_name
    else:
        now = datetime.datetime.now()
        print(
            f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 开始运行 sort and MQ30...")
        start = time.time()
        cmd_1 = f'samtools view -@ {threads_cmd} -bS {out_dir_sample}.sam -o {out_dir_sample}.bam'
        cmd_2 = f'samtools flagstat {out_dir_sample}.bam > {out_dir_sample}.flagstat'
        cmd_3 = f'samtools sort -@ {threads_cmd} -o {out_dir_sample}_sort.bam {out_dir_sample}.bam'
        cmd_4 = f'samtools view -@ {threads_cmd} -b -q 30 {out_dir_sample}_sort.bam -o {out_dir_sample}_sortMq30.bam'
        cmd_5 = f'samtools flagstat {out_dir_sample}_sortMq30.bam > {out_dir_sample}_sortMq30.flagstat'
        cmds = [cmd_1, cmd_2, cmd_3, cmd_4, cmd_5]
        for idx, cmd in enumerate(cmds, start=1):
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            print(cmd)
            now = datetime.datetime.now()
            if proc.returncode != 0:
                end = time.time()
                duration = end - start
                print(
                    f"\033[1;31mStage {idx} of sort and MQ30 failed for {sample_name}.\ncmd: {cmd}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\nTime cost:{duration:.1f}s.\033[0m")
                print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - {sample_name} 运行 sort and MQ30 失败，详见：{out_dir_sample}.sortMq30.flagstat")
                return 1, sample_name
        end = time.time()
        duration = end - start
        now = datetime.datetime.now()
        print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - sort and MQ30 完成，用时{duration:.1f}s...")
        return 0, sample_name


    # 3) calculate reads count using stringTie
def run_stringtie(sample_name: str, fastq_pathway: str, ref_gff: str, threads_cmd: int):
    out_dir_sample = os.path.join(fastq_pathway, sample_name)
    if os.path.exists(f"{out_dir_sample}.gtf") and os.path.exists(f"{out_dir_sample}.gene_abund.tab"):
        now = datetime.datetime.now()
        print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 已存在 stringtie 结果，跳过此步骤...")
        return 0, sample_name
    else:
        now = datetime.datetime.now()
        print(
            f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - 开始运行 stringtie...")
        start = time.time()
        cmd = f'stringtie -p {threads_cmd} {out_dir_sample}_sortMq30.bam -G {ref_gff} -o {out_dir_sample}.gtf -A {out_dir_sample}.gene_abund.tab -e >{out_dir_sample}.stringtie.flagstat'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(cmd)
        now = datetime.datetime.now()
        if result.returncode != 0:
            print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - {sample_name} 运行 stringtie 失败，详见：{out_dir_sample}.stringtie.flagstat")
            return 1, sample_name
        else:
            end = time.time()
            duration = end - start
            print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - stringtie 完成，用时{duration:.1f}s...")
            return 0, sample_name

def run(sample_list_file: str, fastq_pathway: str, ref_gff: str,index: str, threads: int, threads_cmd: int):

    sample_list = []
    with open(sample_list_file, "r") as in_f:
        for line in in_f:
            line = line.strip("\n").split("\t")
            sample_list.append([line[0], line[1], line[2]])

    access_hisat_sample = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for sample in sample_list:
            fastq1, fastq2, sample_name = sample[0], sample[1], sample[2]
            future = executor.submit(
                run_hisat,
                fastq1,
                fastq2,
                sample_name,
                fastq_pathway,
                index,
                threads_cmd,
            )
            futures[future] = sample
        for future in concurrent.futures.as_completed(futures):
            tmp_label = futures[future]
            try:
                result_code, sample_name = future.result()
                if result_code == 0:
                    print(f"{sample_name} 转录组分析结束。")
                    access_hisat_sample.append(sample_name)
                else:
                    print(f"{sample_name} 转录组分析失败。")
            except Exception as e:
                print(f"{tmp_label} 转录组分析失败：{e}")
    
    access_sort_mq30_sample = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for sample_name in access_hisat_sample:
            future = executor.submit(
                run_sort_mq30,
                sample_name,
                fastq_pathway,
                threads_cmd,
            )
            futures[future] = sample_name
        for future in concurrent.futures.as_completed(futures):
            tmp_label = futures[future]
            try:
                result_code, sample_name = future.result()
                if result_code == 0:
                    print(f"{sample_name} 排序及MQ30过滤结束。")
                    access_sort_mq30_sample.append(sample_name)
                else:
                    print(f"{sample_name} 排序及MQ30过滤失败。")
            except Exception as e:
                print(f"{tmp_label} 排序及MQ30过滤失败：{e}")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for sample_name in access_sort_mq30_sample:
            future = executor.submit(
                run_stringtie,
                sample_name,
                fastq_pathway,
                ref_gff,
                threads_cmd,
            )
            futures[future] = sample_name
        for future in concurrent.futures.as_completed(futures):
            tmp_label = futures[future]
            try:
                result_code, sample_name = future.result()
                if result_code == 0:
                    print(f"{sample_name} 基因表达量计算结束。")
                    path = os.path.join(fastq_pathway, f"{sample_name}.sam")
                    if os.path.exists(path):
                        os.remove(path)
                    if os.path.exists(os.path.join(fastq_pathway, f"{sample_name}.bam")):
                        os.remove(os.path.join(fastq_pathway, f"{sample_name}.bam"))
                    if os.path.exists(os.path.join(fastq_pathway, f"{sample_name}_sort.bam")):
                        os.remove(os.path.join(fastq_pathway, f"{sample_name}_sort.bam"))
                else:
                    print(f"{sample_name} 基因表达量计算失败。")
            except Exception as e:
                print(f"{tmp_label} 基因表达量计算失败：{e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--sample_list_file', type=str, required=True,
                        help='Input sample list file with "fastq_1\tfastq_2\tsample_name" format')
    parser.add_argument('-p', '--fastq_pathway', type=str, default='.',
                        help='Input the pathway of fastq files')
    parser.add_argument('-r', '--ref_gff', type=str, required=True,
                        help='Input the reference gff')
    parser.add_argument('-d', '--index', type=str, required=True,
                        help='Input the index of fasta')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Input the number of samples running simultaneously')
    parser.add_argument('-tc', '--threads_cmd', type=int, default=1, help='Input the needed threads number of single sample')

    args = parser.parse_args()

    run(args.sample_list_file, args.fastq_pathway, args.ref_gff, args.index, args.threads, args.threads_cmd)
    