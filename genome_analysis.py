"""
这个脚本用于计算全基因组或甲基化测序过程中参考基因组大小、测序深度。

作者：[常建业]
日期：[2023/07/13]

使用方法：
- 计算参考基因组大小：
  python script.py --reference_fasta_file [参考基因组文件路径] --calculate_genome_size

- 计算测序深度：
  python script.py --reference_fasta_file [参考基因组文件路径] --fastq_file_path_R1 [R1 fastq 文件路径] \
      --fastq_file_path_R2 [R2 fastq 文件路径] --calculate_depth

参数说明：
--reference_fasta_file：参考基因组的FASTA文件路径
--calculate_genome_size：计算参考基因组大小的标志
--fastq_file_path_R1：R1端FASTQ文件的路径
--fastq_file_path_R2：R2端FASTQ文件的路径
--calculate_depth：计算测序深度的标志
"""
import argparse
import gzip
import os
from Bio import SeqIO


def calculate_reference_genome_size(reference_fasta_file):
    total_size = 0

    for record in SeqIO.parse(reference_fasta_file, "fasta"):
        total_size += len(record.seq)

    return total_size


def calculate_gc_content(reference_fasta_file):
    gc_count = 0
    total_count = 0

    for record in SeqIO.parse(reference_fasta_file, "fasta"):
        sequence = record.seq.upper()
        gc_count += sequence.count('G') + sequence.count('C')
        total_count += len(sequence)

    gc_content = (gc_count / total_count) * 100

    return gc_content


def calculate_sequencing_depth(fastq_file_path_R1, fastq_file_path_R2, reference_genome_size):
    total_reads = 0

    with open_fastq_file(fastq_file_path_R1) as fastq_file_R1, open_fastq_file(fastq_file_path_R2) as fastq_file_R2:
        for i, (line_R1, line_R2) in enumerate(zip(fastq_file_R1, fastq_file_R2)):
            if i % 4 == 0:  # 每四行为一个读
                total_reads += 1

    sequencing_depth = (total_reads * 150) / reference_genome_size  # 每条读的长度假设为150

    return sequencing_depth


def open_fastq_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')


def calculate_genome_size(args):
    reference_fasta_file = args.reference_fasta_file
    genome_size = calculate_reference_genome_size(reference_fasta_file)
    print('参考基因组大小：{} bp'.format(genome_size))


def calculate_gc_content_command(args):
    reference_fasta_file = args.reference_fasta_file
    gc_content = calculate_gc_content(reference_fasta_file)
    print('基因组的GC含量：{:.2f}%'.format(gc_content))


def calculate_depth(args):
    fastq_file_path_R1 = args.fastq_file_path_R1
    fastq_file_path_R2 = args.fastq_file_path_R2

    genome_size = calculate_reference_genome_size(args.reference_fasta_file)
    depth = calculate_sequencing_depth(fastq_file_path_R1, fastq_file_path_R2, genome_size)

    file_name_R1, _ = os.path.splitext(os.path.basename(fastq_file_path_R1))
    print('文件名: {}\n测序深度：{:.2f}X'.format(file_name_R1, depth))


if __name__ == '__main__':
    # 创建解析器对象
    parser = argparse.ArgumentParser(
        description='Calculate reference genome size, sequencing depth and GC content.')

    # 添加命令行参数...
    parser.add_argument('--reference_fasta_file', help='Path to the reference FASTA file')
    parser.add_argument('--fastq_file_path_R1', help='Path to the FASTQ file (R1)')
    parser.add_argument('--fastq_file_path_R2', help='Path to the FASTQ file (R2)')
    parser.add_argument('--calculate_all', action='store_true',
                        help='Calculate reference genome size, sequencing depth, and GC content')
    parser.add_argument('--calculate_genome_size', action='store_true', help='Calculate reference genome size')
    parser.add_argument('--calculate_gc_content', action='store_true', help='Calculate GC content')
    parser.add_argument('--calculate_depth', action='store_true', help='Calculate sequencing depth')

    # 解析命令行参数
    args = parser.parse_args()

    if args.calculate_all:  # 如果选择计算所有参数
        # 检查是否提供了所需的文件路径
        if not args.reference_fasta_file or not args.fastq_file_path_R1 or not args.fastq_file_path_R2:
            print('错误：请提供参考基因组文件以及双端FASTQ文件的路径（R1和R2）')
            exit(1)
        calculate_genome_size(args)  # 计算基因组大小
        calculate_gc_content_command(args)  # 计算GC含量
        calculate_depth(args)  # 计算测序深度
    elif args.calculate_genome_size:  # 如果只选择计算基因组大小
        calculate_genome_size(args)  # 计算基因组大小
    elif args.calculate_gc_content:  # 如果只选择计算GC含量
        calculate_gc_content_command(args)  # 计算GC含量
    elif args.calculate_depth:  # 如果只选择计算测序深度
        # 检查是否提供了所需的文件路径
        if not args.fastq_file_path_R1 or not args.fastq_file_path_R2:
            print('错误：请提供双端FASTQ文件的路径（R1和R2）')
            exit(1)
        calculate_depth(args)  # 计算测序深度
    else:
        exit(1)