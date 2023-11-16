#!/usr/bin/python
#coding：utf-8
"""
该脚本用于对单端或双端测序数据进行质量控制、修剪、比对、去重和峰识别等分析流程。
作者：Jianye Chang
日期：2023-07-19
使用方法：
python data_analysis.py [fastq文件1] [fastq文件2（可选）]
参数说明：
- fastq文件1：需要进行分析的第一个fastq文件路径。
- fastq文件2（可选）：需要进行分析的第二个fastq文件路径，仅适用于双端测序数据。
注意事项：
- 请提前安装所需的工具（fastp、samtools、bowtie2、picard、sambamba、alignmentSieve、bamCoverage、MACS2）并设置对应路径。
- 确保参考基因组已经构建并设定正确的路径。
- 分析结果将会保存在指定输出目录中。
"""
import os  # 导入os模块，用来进行文件和目录操作
import sys  # 导入sys模块，用于从命令行获取参数
import subprocess  # 导入subprocess模块，用于执行外部命令

# 定义各个工具的路径
fastp = "/public/home/changjianye/anaconda3/envs/cuttag/bin/fastp"  # fastp工具路径
samtools = "/public/software/env01/bin/samtools"  # samtools工具路径
bowtie2 = "/public/home/changjianye/anaconda3/envs/cuttag/bin/bowtie2"  # bowtie2工具路径
genome = "/public/home/changjianye/project/hzh/Oryza/Oryzabowtie2index/Oryza"  # 参考基因组路径

# 获取输入的fastq文件名
fq_1 = sys.argv[1]  # 第一个fastq文件路径
fq_2 = sys.argv[2]  # 第二个fastq文件路径

# 定义chrM名称
chrM = "NC_011033.1"

# 提取文件名作为输出文件名
samplename = os.path.splitext(os.path.basename(fq_1))[0]

# 指定使用的核心数
cores = 10

# 创建输出目录
output_dir = f"/public/home/changjianye/project/hzh/ATAC/{samplename}/"
os.makedirs(output_dir + "fastp", exist_ok=True)
os.makedirs(output_dir + "bam", exist_ok=True)
os.makedirs(output_dir + "MACS2", exist_ok=True)
os.makedirs(output_dir + "bigwig", exist_ok=True)

# Quality control and read trimming by fastp
print(f"Starting QC and trimming for {samplename}")
subprocess.run([
    fastp,
    "-i", fq_1,
    "-I", fq_2,
    "-o", output_dir + "fastp/trimmed_" + samplename + "_R1.fastq.gz",
    "-O", output_dir + "fastp/trimmed_" + samplename + "_R2.fastq.gz",
    "-h", output_dir + "fastp/" + samplename + "_fastp.html",
    "-j", output_dir + "fastp/" + samplename + "_fastp.json",
    "-w", str(cores),
    "--cut_tail",
    "--cut_tail_window_size=1",
    "--cut_tail_mean_quality=30",
    "--average_qual=30",
    "--length_required=20"
])

# Map reads to reference genome by bowtie2
print(f"Starting mapping for {samplename}")
subprocess.run([
    bowtie2,
    "--very-sensitive",
    "-p", str(cores),
    "-x", genome,
    "-1", output_dir + "fastp/trimmed_" + samplename + "_R1.fastq.gz",
    "-2", output_dir + "fastp/trimmed_" + samplename + "_R2.fastq.gz",
    "-S", output_dir + "bam/" + samplename + ".sam"
])

# Compress sam files to bam files
subprocess.run([
    samtools, "view", "-S", "-@", str(cores), "-b",
    output_dir + "bam/" + samplename + ".sam",
    "-o", output_dir + "bam/" + samplename + ".bam"
])

# Remove mitochondrial reads after alignment
subprocess.run([
    samtools, "view", "-h", output_dir + "bam/" + samplename + ".bam",
    "|", "grep", "-v", chrM,
    "|", samtools, "sort", "-@", str(cores), "-O", "bam",
    "-o", output_dir + "bam/" + samplename + "_sorted_rmChrM.bam"
])

# Remove duplicates using picard
subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/picard", "MarkDuplicates",
    "I=" + output_dir + "bam/" + samplename + "_sorted_rmChrM.bam",
    "O=" + output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup.bam",
    "M=" + output_dir + "bam/" + samplename + "_rmDup_metrics.txt",
    "REMOVE_DUPLICATES=true"
])

# Sort BAM files by genomic coordinates
os.rename(
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup.bam",
    output_dir + "bam/" + samplename + "_rmChrM_rmDup.bam"
)
subprocess.run([
    samtools, "sort", "-@", str(cores),
    "-o", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup.bam",
    output_dir + "bam/" + samplename + "_rmChrM_rmDup.bam"
])

# Filter and keep the uniquely mapped reads
subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/sambamba", "view",
    "-h", "-t", str(cores), "-f", "bam",
    "-F", "[XS] == null and not unmapped",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup.bam",
    ">", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped.bam"
])

# Index BAM files
subprocess.run([
    samtools, "index",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped.bam"
])

# Shift BAM files using deeptools alignmentSieve
subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/alignmentSieve", "-b",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped.bam",
    "-p", str(cores), "--ATACshift",
    "-o", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_shift.bam"
])

# Generate bigwig files
os.rename(
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_shift.bam",
    output_dir + "bam/" + samplename + "_rmChrM_rmDup_mapped_rmbl_shift.bam"
)
subprocess.run([
    samtools, "sort", "-@", "4",
    "-o", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam",
    output_dir + "bam/" + samplename + "_rmChrM_rmDup_mapped_rmbl_shift.bam"
])
subprocess.run([
    samtools, "index",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam"
])

subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/bamCoverage", "--bam",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl.bam",
    "--numberOfProcessors", str(cores), "--binSize", "10",
    "--normalizeUsing", "RPKM",
    "-o", output_dir + "bigwig/" + samplename + "_RPKM_normalized.bw"
])

subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/bamCoverage", "--bam",
    output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam",
    "--numberOfProcessors", str(cores), "--binSize", "10",
    "--normalizeUsing", "RPKM",
    "-o", output_dir + "bigwig/" + samplename + "_RPKM_normalized_shift.bw"
])

# Peak calling by MACS2
genome_size = "1130276682"

subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/macs2", "callpeak",
    "-t", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl.bam",
    "-g", genome_size,
    "-n", samplename + "_noShift",
    "-q", "0.01",
    "--outdir", output_dir + "MACS2",
    "--shift", "-100",
    "--extsize", "200",
    "--nomodel", "-B", "--SPMR",
    "--call-summits"
])

subprocess.run([
    "/public/home/changjianye/anaconda3/envs/cuttag/bin/macs2", "callpeak",
    "-t", output_dir + "bam/" + samplename + "_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam",
    "-g", genome_size,
    "-n", samplename + "_shift",
    "-q", "0.01",
    "--outdir", output_dir + "MACS2",
    "--shift", "0",
    "--extsize", "200",
    "--nomodel", "-B", "--SPMR",
    "--call-summits"
])

print("Analysis complete.")
``