#!/usr/bin/env python

import os
import sys
import time
import argparse
import datetime
import subprocess

# 原始数据路径
input_dir = '/public/home/changjianye/WGS'

# 输出基础路径
base_output_dir = '/public/home/changjianye/WGS'

# 软件路径
fastqc = '/public/home/changjianye/miniconda3/envs/atac/bin/fastqc'
fastp = '/public/home/changjianye/miniconda3/envs/atac/bin/fastp'
trim_galore = '/public/home/changjianye/miniconda3/envs/atac/bin/trim_galore'
samtools = '/public/home/changjianye/miniconda3/envs/atac/bin/samtools'
bwa = '/public/software/bin/bwa'
tabix = '/public/home/changjianye/miniconda3/envs/atac/bin/tabix'
bgzip = '/public/home/changjianye/miniconda3/envs/atac/bin/bgzip'
gatk = '/public/home/changjianye/miniconda3/envs/atac/bin/gatk'

# CPU数
cores = 30

# 打印时间
print(f"{datetime.datetime.now().strftime('%Y/%m/%d_%H:%M:%S')}")
# 记录所有
START = time.perf_counter()

# 基因组路径
bwa_ref = "/public/home/changjianye/goose_ref/AnsCyg.fa"

# 判断索引文件
if not os.path.isfile(f"{bwa_ref}.fai"):
    print("Genome index not in genome dir ! Now , create it !!!")
    subprocess.call([f"{samtools}", "faidx", f"{bwa_ref}"])
else:
    print("file exists.")

parser = argparse.ArgumentParser(description='Process some files.')
parser.add_argument('--fq1', help='input fastq file 1')
parser.add_argument('--fq2', help='input fastq file 2')
parser.add_argument('--samplename', help='sample name')

args = parser.parse_args()

fq_1 = args.fq1
fq_2 = args.fq2
samplename = args.samplename

def get_samplename(fq_1):
    filename = os.path.basename(fq_1)
    if filename.endswith('.fastq.gz'):
        samplename = os.path.splitext(filename)[0].replace('_1.fastq.gz', '')
    elif filename.endswith('.fq.gz'):
        samplename = os.path.splitext(filename)[0].replace('_1.fq.gz', '')
    else:
        print("Please check the fastq file suffix !!!")
        return False

    return samplename

#samplename = get_samplename(fq_1)
print(f"Sample name is {samplename}")

output_dir = os.path.join(base_output_dir, samplename)
os.makedirs(os.path.join(output_dir, 'fastp'), exist_ok=True)
os.makedirs(os.path.join(output_dir, 'trim'), exist_ok=True)
os.makedirs(os.path.join(output_dir, 'mapping'), exist_ok=True)

trim_out = os.path.join(output_dir, 'trim')
mapping_outdir = os.path.join(output_dir, 'mapping')
fastp_out = os.path.join(output_dir, 'fastp')

# fastp
# subprocess.call([
#     f"{fastp}",
#     "-i", fq_1,
#     "-o", f"{trim_out}/{samplename}_fastp_1.fq.gz",
#     "-I", fq_2,
#     "-O", f"{trim_out}/{samplename}_fastp_2.fq.gz",
#     "--trim_poly_x",
#     "-h", f"{fastp_out}/{samplename}_fastp.html",
#     "-j", f"{fastp_out}/{samplename}_fastp.json",
#     "-w", str(cores)
# ])

# # 去接头
# subprocess.call([
#     f"{trim_galore}", "-output_dir", f"{trim_out}", "--paired", "--length", "25", "--quality", "25",
#     "--stringency", "5", fq_1, fq_2, "--gzip"
# ])

clean_fq1 = f"{trim_out}/{samplename}_fastp_1.fq.gz"
clean_fq2 = f"{trim_out}/{samplename}_fastp_2.fq.gz"

label = f"@RG\\tID:{samplename}\\tPL:MGI\\tSM:{samplename}"  # 可以修改
# 运行BWA
# 创建第一个子进程，执行 bwa mem 命令，并将其输出写入管道
p1 = subprocess.Popen([
    bwa, "mem", "-t", str(cores), "-R", label, bwa_ref, clean_fq1, clean_fq2
], stdout=subprocess.PIPE)

# 创建第二个子进程，执行 samtools view 命令，并将 p1 的输出作为输入
p2 = subprocess.Popen([
    samtools, "view", "-Sb", "-", "-o", f"{mapping_outdir}/{samplename}.bam"
], stdin=p1.stdout)

# 等待第二个子进程完成
p2.communicate()

print("** bwa mapping done **")

bam = f"{mapping_outdir}/{samplename}.bam"

# BAM排序
subprocess.call([
    f"{samtools}", "sort", "-@", str(cores), "-O", "bam", "-o",
    f"{os.path.splitext(bam)[0]}.sorted.bam", f"{bam}"
])
print("** BAM sort done **")
sorted_bam = f"{os.path.splitext(bam)[0]}.sorted.bam"

# 标记PCR重复
subprocess.call([
    f"{gatk}", "MarkDuplicates", "-I", sorted_bam, "-O", f"{os.path.splitext(sorted_bam)[0]}.markdup.bam", "-M",
    f"{os.path.splitext(sorted_bam)[0]}.markdup_metrics.txt", "--REMOVE_DUPLICATES", "true"
])
print("** markdup done **")
markdup_bam = f"{os.path.splitext(sorted_bam)[0]}.markdup.bam"

# 制作index
subprocess.call([f"{samtools}", "index", markdup_bam])
print("** index done **")

# GATK基因组路径
gatk_ref = f"{bwa_ref}.fa"

#制作字典
subprocess.call([
    f"{gatk}", "CreateSequenceDictionary", "-R", gatk_ref
])
print("** dict done **")

# gvcf
subprocess.call([
    f"{gatk}", "HaplotypeCaller", "-R", gatk_ref, "--emit-ref-confidence", "GVCF", "-I", markdup_bam, "-O",
    f"{os.path.splitext(markdup_bam)[0]}.g.vcf"
])
print("** g.vcf done **")
gvcf = f"{os.path.splitext(markdup_bam)[0]}.g.vcf"

# vcf
subprocess.call([
    f"{gatk}", "GenotypeGVCFs", "-R", gatk_ref, "-V", gvcf, "-O", f"{os.path.splitext(gvcf)[0]}.vcf"
])
print("** vcf done **")
vcf = f"{os.path.splitext(gvcf)[0]}.vcf"

# bgzip
subprocess.call([f"{bgzip}", "-f", vcf])
vcf_gz = f"{vcf}.gz"

# tabix
subprocess.call([f"{tabix}", "-p", "vcf", vcf_gz])

# 使用SelectVariants，选出Indel
subprocess.call([
    f"{gatk}", "SelectVariants", "-select-type", "INDEL", "-V", vcf_gz, "-O",
    f"{os.path.splitext(vcf_gz)[0]}.indel.vcf.gz"
])

# 为Indel作过滤
subprocess.call([
    f"{gatk}", "VariantFiltration", "-V", f"{os.path.splitext(vcf_gz)[0]}.indel.vcf.gz", "--filter-expression",
    "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", "--filter-name", "Filter", "-O", f"{os.path.splitext(vcf_gz)[0]}.indel.filter.vcf.gz"
])

# 使用SelectVariants，选出SNP
subprocess.call([
    f"{gatk}", "SelectVariants", "-select-type", "SNP", "-V", vcf_gz, "-O", f"{os.path.splitext(vcf_gz)[0]}.snp.vcf.gz"
])

# 为SNP作硬过滤
subprocess.call([
    f"{gatk}", "VariantFiltration", "-V", f"{os.path.splitext(vcf_gz)[0]}.snp.vcf.gz", "--filter-expression",
    "QD < .0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", "--filter-name", "Filter", "-O", f"{os.path.splitext(vcf_gz)[0]}.snp.filter.vcf.gz"
])

# 打印程序运行时间
END = time.perf_counter()
print(f"Finished in {round(END-START, 2)} seconds")