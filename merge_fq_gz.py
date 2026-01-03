# /usr/bin python3

# -*- coding: utf-8 -*-
# @time    : 2025/12/24 16:46
# @author  : HeJunyi
# @file    : merge_fq_gz.py

import subprocess
import re
import concurrent.futures
import argparse


def merge_fq_gz(fq_1: str, fq_2: str, sample_id: str):
    sample_fastq = sample_id + ".fq"
    cmd = f"zcat {fq_1} > {sample_fastq} && zcat {fq_2} >> {sample_fastq} && gzip -c {sample_fastq} > {sample_fastq}.gz"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(cmd)
    if result.returncode != 0:
        print(f"{sample_id} failed.")
        return 1, sample_id
    else:
        print(f"{sample_id} success.")
        return 0, sample_id

def run(input_file: str, threads: int):
    sample_list = []
    sample_dict = {}
    with open(input_file, "r") as in_f:
        for line in in_f:
            line = line.strip()  # 去除首尾空白字符

            if line.startswith("冲突"):
                # 匹配冲突行，提取二级和三级编码
                match = re.search(r'二级编码 (\S+) 和三级编码 (\S+)', line)
                if match:
                    secondary_code = match.group(1)
                    third_code = match.group(2)
                    current_sample_id = f"{secondary_code}_{third_code}"

            elif line.startswith("文件1"):
                # 匹配文件1行
                match = re.search(r'文件1: (\S+)（模式）', line)
                if match and current_sample_id:
                    if current_sample_id not in sample_dict:
                        sample_dict[current_sample_id] = []
                    sample_dict[current_sample_id].append(match.group(1))

            elif line.startswith("文件2"):
                # 匹配文件2行
                match = re.search(r'文件2: (\S+)', line)  # 注意：文件中没有"（模式）"
                if match and current_sample_id:
                    if current_sample_id not in sample_dict:
                        sample_dict[current_sample_id] = []
                    sample_dict[current_sample_id].append(match.group(1))
                    # 重置当前sample_id，准备下一组
                    current_sample_id = None

        # 输出结果
        for k, v in sample_dict.items():
            sample_list.append([v[0], v[1], k])

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for sample in sample_list:
            future = executor.submit(merge_fq_gz,
                                     sample[0],
                                     sample[1],
                                     sample[2])
            futures[future] = sample[2]
        for future in concurrent.futures.as_completed(futures):
            tmp_label = futures[future]
            try:
                result_code, sample_id = future.result()
                if result_code == 0:
                    print(f"{sample_id} success.")
                else:
                    print(f"{sample_id} failed.")
            except Exception as e:
                print(f"Error in {tmp_label}: {e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='Input check file')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Input the number of samples running simultaneously')

    args = parser.parse_args()

    run(args.input_file, args.threads)

