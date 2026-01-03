# /usr/bin python3

# -*- coding: utf-8 -*-
# @time    : 2026/01/02 15:19
# @author  : HeJunyi
# @file    : auto_fastp.py

import os
import datetime
import time
import subprocess
import concurrent.futures
import argparse

def run_fastp(fastq1: str, fastq2: str, sample_name: str, fastq_pathway: str):
    fastq1_file = os.path.join(fastq_pathway, f"{fastq1}")
    fastq2_file = os.path.join(fastq_pathway, f"{fastq2}")
    fastq1_deal = os.path.join(fastq_pathway, f'{fastq1.replace("fastq.gz", "trim.fastq.gz").replace("fq.gz", "trim.fq.gz")}')
    fastq2_deal = os.path.join(fastq_pathway, f'{fastq2.replace("fastq.gz", "trim.fastq.gz").replace("fq.gz", "trim.fq.gz")}')
    now = datetime.datetime.now()
    print(f"{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}: {sample_name} - {sample_name} 开始运行 fastp ...")
    cmd = f"fastp -j {sample_name}.json -h {sample_name}.html -c -q 30 -u 20 -e 30 -l 40 -w 40 -i {fastq1_file} -I {fastq2_file} -o {fastq1_deal} -O {fastq2_deal} 2> {sample_name}_trim.log"
    print(cmd)
    start = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"fastp deal with {fastq1_file} failed")
        return 1, sample_name
    else:
        end = time.time()
        duration = end - start
        now = datetime.datetime.now()
        print(f"fastp deal with {fastq1_file} success 用时{duration:.1f}s。时间：{now.year}/{now.month}/{now.day} {now.hour}:{now.minute}。")
        return 0, sample_name

def run(input_sample_list_file: str, fastq_pathway: str, threads: int):
    sample_list = []
    with open(input_sample_list_file, "r") as sample_file:
        for line in sample_file:
            line = line.strip("\n").split("\t")
            sample_list.append([line[0], line[1], line[2]])

    futures = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for sample in sample_list:
            future = executor.submit(run_fastp,
                                     sample[0],
                                     sample[1],
                                     sample[2],
                                     fastq_pathway)
            futures[future] = sample

        for future in concurrent.futures.as_completed(futures):
            tmp_label = sample[2]
            try:
                result_code, sample = future.result()
                if result_code == 0:
                    print(f"{sample} was dealed.")
                else:
                    print(f"{sample} failed.")
            except Exception as e:
                print(f"Error in fastp of {tmp_label}: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_sample_list_file', type=str, required=True,
                        help='Input sample list file with "fastq_1\tfastq_2\tsample_name" format')
    parser.add_argument('-p', '--fastq_pathway', type=str, default='.',
                        help='Input the pathway of fastq files')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Input the number of samples running simultaneously')

    args = parser.parse_args()

    run(args.input_sample_list_file, args.fastq_pathway, args.threads)
