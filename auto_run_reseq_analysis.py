# /usr/bin python3

# -*- coding: utf-8 -*-
# @time    : 2025/12/12 16:46
# @author  : HeJunyi
# @file    : auto_run_reseq_analysis.py

import os
import datetime
import time
import subprocess
import concurrent.futures
import argparse


def process_create_index(ref_genome: str, software: str):
    ref_0123 = ref_genome + ".0123"
    ref_amb = ref_genome + ".amb"
    ref_ann = ref_genome + ".ann"
    ref_bwt_2bit_64 = ref_genome + ".bwt.2bit.64"
    ref_pac = ref_genome + ".pac"
    ref_fai = ref_genome + ".fai"
    ref_dict = ref_genome.replace('.fa', '').replace('.fasta', '') + '.dict'
    files = [ref_0123, ref_amb, ref_ann, ref_bwt_2bit_64, ref_pac, ref_fai, ref_dict]
    if all(os.path.exists(f) for f in files):
        print("所有索引文件均存在")
        return 0
    else:
        missing = [f for f in files if not os.path.exists(f)]
        print("缺失文件：", missing)
        for i in files:
            if i not in missing and os.path.isfile(i):
                os.remove(i)
        cmd = f"bwa-mem2 index {ref_genome} && samtools faidx {ref_genome} && {software} CreateSequenceDictionary -R {ref_genome}"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(cmd)
        if proc.returncode != 0:
            print(f"\033[1;31mCreating index failed for {ref_genome}.\ncmd: {cmd}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\033[0m")
            return 1
        else:
            print(f"\033[1;32mCreating index is finished for {ref_genome}.\033[0m")
            return 0

def process_mapping(fastq_1: str, fastq_2: str, ref_genome: str, sample_name: str, pathway: str, threads_cmd: int):

    fastq_1 = os.path.join(pathway, fastq_1)
    fastq_2 = os.path.join(pathway, fastq_2)
    out_prefix = os.path.join(pathway, sample_name)
    mkdup_bam = out_prefix + ".fixmate.sorted.uniq.mkdup.bam"
    mkdup_index = mkdup_bam + ".bai"
    now = datetime.datetime.now()

    # 如果最终索引存在，认为 mapping 已完成
    if os.path.exists(mkdup_index) and os.path.exists(mkdup_bam):
        print(f"\033[1;32mThe mapping of {sample_name} existed ({now.hour}:{now.minute}:{now.second}).\033[0m")
        return 0, mkdup_bam

    # Ensure pathway exists
    os.makedirs(pathway, exist_ok=True)

    start = time.time()

    # 定义中间与输出文件
    fixmate_bam = out_prefix + ".fixmate.bam"
    flagstat_file = out_prefix + ".flagstat.txt"
    tmp_prefix = out_prefix + ".tmp"
    sorted_uniq_bam = out_prefix + ".fixmate.sorted.uniq.bam"

    # 1) bwa -> fixmate -> 输出 BAM (fixmate_bam)
    # note: samtools fixmate does not support -@ threads; remove that option
    cmd1 = (
        f'bwa-mem2 mem -M -t {threads_cmd} -R "@RG\\tID:{sample_name}\\tPL:illumina\\tSM:{sample_name}\\tCN:BERRY" '
        f'{ref_genome} {fastq_1} {fastq_2} '
        f'| samtools fixmate -m - {fixmate_bam}'
    )

    # 2) flagstat on fixmate BAM
    cmd2 = f'samtools flagstat -@ {threads_cmd} {fixmate_bam} > {flagstat_file}'

    # 3) sort and filter (sort -> view) in a pipe to avoid extra full sorted BAM
    cmd3 = (
        f'samtools sort -@ {threads_cmd} -T {tmp_prefix} {fixmate_bam} '
        f'| samtools view -F 4 -q 1 -@ {threads_cmd} -b -o {sorted_uniq_bam} -'
    )

    # 4) markdup
    cmd4 = f'samtools markdup -@ {threads_cmd} {sorted_uniq_bam} {mkdup_bam}'

    # 5) index (generate standard BAI). Do NOT use -c (which creates CSI) unless you need CSI.
    # Let samtools write the index file automatically (sample.bam.bai).
    cmd5 = f'samtools index -@ {threads_cmd} {mkdup_bam}'

    # 按序执行命令，遇到失败立刻返回错误
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5]
    for idx, cmd in enumerate(cmds, start=1):
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(cmd)
        if proc.returncode != 0:
            end = time.time()
            duration = end - start
            print(f"\033[1;31mStage {idx} of mapping failed for {sample_name}.\ncmd: {cmd}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\nTime cost:{duration:.1f}s.\033[0m")
            return 1, ""

    # 清理中间文件，保留最终 mkdup_bam (.bai 已生成) 和 flagstat
    try:
        for f in [fixmate_bam, sorted_uniq_bam]:
            if os.path.exists(f):
                os.remove(f)
    except Exception:
        pass

    end = time.time()
    duration = end - start
    print(f"\033[1;32mThe mapping of {sample_name} is finished. ({now.hour}:{now.minute}:{now.second}\tcost:{duration:.1f}s).\033[0m")
    return 0, mkdup_bam


def process_calling(ref_genome: str, bam_path: str, pathway: str, software: str):
    # bam_path is full path to mkdup BAM
    now = datetime.datetime.now()
    if not os.path.exists(bam_path):
        print(f"\033[1;31mBAM for calling not found: {bam_path}\033[0m")
        return 1, ""

    sample_name = os.path.basename(bam_path).replace('.fixmate.sorted.uniq.mkdup.bam', '')
    out_file = os.path.join(pathway, f"{sample_name}.gatk.raw.gvcf.gz")

    if os.path.exists(out_file):
        print(f"\033[1;32mThe calling of {sample_name} existed ({now.hour}:{now.minute}:{now.second}).\033[0m")
        return 0, out_file

    os.makedirs(pathway, exist_ok=True)
    start = time.time()
    cmd = f'{software} HaplotypeCaller -R {ref_genome} -I {bam_path} -O {out_file} --output-mode EMIT_ALL_CONFIDENT_SITES -ERC GVCF > {sample_name}_gatk_HaplotypeCaller.txt 2>&1'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(cmd)
    end = time.time()
    duration = end - start
    if proc.returncode != 0:
        print(f"\033[1;31mCalling failed for {sample_name}.\ncmd: {cmd}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\nTime cost:{duration:.1f}s\033[0m")
        return 1, ""
    else:
        print(f"\033[1;32mThe calling of {sample_name} is finished. ({now.hour}:{now.minute}:{now.second}\tcost:{duration:.1f}s).\033[0m")
        return 0, out_file


def process_merge_filter(sample_list: list, pathway: str, software: str, ref_genome: str):
    merge_gatk_out = os.path.join(pathway, "merge.gatk.raw.gvcf.gz")
    vcf_gatk_out = os.path.join(pathway, "merge.gatk.raw.vcf.gz")
    extract_snp_out = os.path.join(pathway, "gatk.snp.raw.vcf.gz")
    gatk_filter_snp_out = os.path.join(pathway, "gatk.snp.msk.vcf.gz")
    filter_snp_out = os.path.join(pathway, "gatk.snp.flt.vcf.gz")
    extract_indel_out = os.path.join(pathway, "gatk.indel.raw.vcf.gz")
    gatk_filter_indel_out = os.path.join(pathway, "gatk.indel.msk.vcf.gz")
    filter_indel_out = os.path.join(pathway, "gatk.indel.flt.vcf.gz")
    now = datetime.datetime.now()
    start = time.time()

    with open(os.path.join(pathway, "sample.list"), "w") as out:
        for sample in sample_list:
            if os.path.exists(sample):
                out.write(f"{sample}\n")

    # 1) merge gvcf
    cmd1 = f"{software} CombineGVCFs -R {ref_genome} --variant {os.path.join(pathway, 'sample.list')} -O {merge_gatk_out} > gatk_CombineGVCFs.txt 2>&1"

    # 2) gvcf to vcf
    cmd2 = f"{software} GenotypeGVCFs -R {ref_genome} -stand-call-conf 30.0 --variant {merge_gatk_out} -O {vcf_gatk_out} > gatk_GenotypeGVCFs.txt 2>&1"

    # 3) extract SNPs
    cmd3 = f"{software} SelectVariants -R {ref_genome} --variant {vcf_gatk_out} -O {extract_snp_out} -select-type SNP > gatk_SelectVariants_SNP.txt 2>&1"

    # 4) SNPs VariantFiltration
    cmd4 = (f'{software} VariantFiltration -R {ref_genome} '
        '--filter-name "FilterQD" --filter-expression "QD < 2.0" '
        '--filter-name "FilterMQ" --filter-expression "MQ < 40.0" '
        '--filter-name "FilterFS" --filter-expression "FS > 60.0" '
        '--filter-name "FilterSOR" --filter-expression "SOR > 3.0" '
        '--filter-name "FilterMQRankSum" --filter-expression "MQRankSum < -12.5" '
        '--filter-name "FilterReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" '
        f'--variant {extract_snp_out} -O {gatk_filter_snp_out} > gatk_VariantFiltration_SNP.txt 2>&1')

    # 5) SNPs Filter
    cmd5 = f'zgrep -v Filter {gatk_filter_snp_out} | bgzip -c > {filter_snp_out} && tabix -p vcf {filter_snp_out}'

    # 6) extract Indels
    cmd6 = f'{software} SelectVariants -R {ref_genome} --variant {vcf_gatk_out} -O {extract_indel_out} -select-type INDEL > gatk_SelectVariants_Indel.txt 2>&1'

    # 7) Indel VariantFiltration
    cmd7 = (f'{software} VariantFiltration -R {ref_genome} '
        '--filter-name "FilterQD" --filter-expression "QD < 2.0" '
        '--filter-name "FilterFS" --filter-expression "FS > 200.0" '
        '--filter-name "FilterSOR" --filter-expression "SOR > 10.0" '
        '--filter-name "FilterReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0" '
        f'--variant {extract_indel_out} -O {gatk_filter_indel_out} > gatk_VariantFiltration_Indel.txt 2>&1')

    # 8) Indels Filter
    cmd8 = f'zgrep -v Filter {gatk_filter_indel_out} | bgzip -c > {filter_indel_out} && tabix -p vcf {filter_indel_out}'

    # 按序执行命令，遇到失败立刻返回错误
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8]

    for idx, cmd in enumerate(cmds, start=1):
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(cmd)
        if proc.returncode != 0:
            end = time.time()
            duration = end - start
            print(f"\033[1;31mStage {idx} of merge_and_filter failed.\ncmd: {cmd}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\nTime cost:{duration:.1f}s\033[0m")
            return 1

    end = time.time()
    duration = end - start
    print(f"\033[1;32mThe merge_and_filter is finished. ({now.hour}:{now.minute}:{now.second}\tcost:{duration:.1f}s).\033[0m")
    return 0


def run(input_sample_list_file: str, ref_genome: str, pathway: str, threads: int, threads_cmd: int, software: str):
    sample_list = []
    max_t = threads * threads_cmd
    if max_t > os.cpu_count():
        print(f"\033[1;31mWarning: The total threads number {max_t} exceeds the CPU core count {os.cpu_count()}, adjusting to {os.cpu_count()}.\033[0m")
        max_t = os.cpu_count()

    result_code = process_create_index(ref_genome, software)
    if result_code == 0:
        with open(input_sample_list_file, 'r') as sample_list_file:
            for line in sample_list_file:
                line = line.strip("\n").split("\t")
                fastq_1, fastq_2, sample_name = line[0], line[1], line[2]
                sample_list.append([fastq_1, fastq_2, sample_name])

        access_mapping_list = []
        # Submit mapping jobs, collect as they finish (to keep concurrency),
        # but store results and then rebuild list in original sample_list order.
        futures = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for sample in sample_list:
                future = executor.submit(process_mapping,
                                         sample[0],
                                         sample[1],
                                         ref_genome,
                                         sample[2],
                                         pathway,
                                         threads_cmd)
                futures[future] = sample

            mapping_results = {}  # sample_name -> bam_path
            for future in concurrent.futures.as_completed(futures):
                sample = futures[future]
                tmp_label = sample
                try:
                    result_code, bam_path = future.result()
                    sample_name = sample[2]
                    if result_code == 0 and bam_path and os.path.exists(bam_path):
                        mapping_results[sample_name] = bam_path
                    else:
                        print(f"\033[1;31mFailed in mapping: {tmp_label}\033[0m\n")
                except Exception as e:
                    print(f"Error in mapping of {tmp_label}: {e}")

        # Rebuild access_mapping_list in the original sample_list order
        for sample in sample_list:
            sample_name = sample[2]
            if sample_name in mapping_results:
                access_mapping_list.append(mapping_results[sample_name])

        access_calling_file_list = []
        # Submit calling jobs, collect as they finish, then rebuild ordered list
        futures = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_t) as executor:
            for bam_path in access_mapping_list:
                future = executor.submit(process_calling,
                                         ref_genome,
                                         bam_path,
                                         pathway,
                                         software)
                futures[future] = bam_path

            calling_results = {}  # bam_path -> sample_file
            for future in concurrent.futures.as_completed(futures):
                bam = futures[future]
                try:
                    result_code, sample_file = future.result()
                    if result_code == 0 and sample_file and os.path.exists(sample_file):
                        calling_results[bam] = sample_file
                    else:
                        print(f"\033[1;31mFailed in calling: {bam}\033[0m\n")
                except Exception as e:
                    print(f"Error in calling of {bam}: {e}")

        # Rebuild access_calling_file_list in the original access_mapping_list order
        for bam in access_mapping_list:
            if bam in calling_results:
                access_calling_file_list.append(calling_results[bam])

        try:
            result_code = process_merge_filter(access_calling_file_list, pathway, software, ref_genome)
            if result_code == 0:
                print("\033[1;32mRun completed.\033[0m")
            else:
                print("\033[1;31mRun failed.\033[0m")
        except Exception as e:
            print(f"Error: {e}")
    else:
        print("\033[1;31mRun failed in creating index.\033[0m")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_sample_list_file', type=str, required=True,
                        help=r'Input sample list file with "fastq_1\tfastq_2\tsample_name" format')
    parser.add_argument('-r', '--ref_genome', type=str, required=True,
                        help='Input reference genome')
    parser.add_argument('-p', '--pathway', type=str, default='.',
                        help='Input the pathway of fastq files')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Input the number of samples running simultaneously')
    parser.add_argument('-tc', '--threads_cmd', type=int, default=1, help='Input the needed threads number of single sample in the mapping process')
    parser.add_argument('-s', '--software', type=str, required=True, help='Input the designated version of gatk')

    args = parser.parse_args()

    # 验证参数
    if not os.path.exists(args.input_sample_list_file):
        print(f"\033[1;31mError: Sample list file {args.input_sample_list_file} does not exist\033[0m")
        exit(1)

    if args.threads < 1:
        print(f"\033[1;31mWarning: Threads value {args.threads} is invalid, setting to 1\033[0m")
        args.threads = 1

    run(args.input_sample_list_file, args.ref_genome, args.pathway, args.threads, args.threads_cmd, args.software)
