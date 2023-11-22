"""
Bismark转MethyKit格式脚本

作者：[常建业]
日期：[2023/07/13]

将Bismark格式的甲基化测序结果文件转换为MethyKit格式的文件。

用法:
    python convert_bismark_to_methykit.py <输入文件>

参数:
    <输入文件> (str): Bismark格式的甲基化测序结果文件

输出:
    生成三个MethyKit格式的输出文件，分别对应 CG 位点、CHG 位点和 CHH 位点。

注意:
    1. 输入文件可以是普通文本文件或压缩的 .gz 文件。
    2. 要求输入文件的列顺序如下:
       - 第一列: 染色体号
       - 第二列: 位置
       - 第三列: 碱基链向信息
       - 第四列: 甲基化的读数
       - 第五列: 未甲基化的读数
       - 第六列: CpG上下文类型 (CG/CHG/CHH)

示例:
    python convert_bismark_to_methykit.py input.txt
"""
import sys
import gzip

def convert_bismark_to_methykit(infile):
    # 定义输出文件名
    out1_file = f"{infile}_CG_methykit.txt"
    out2_file = f"{infile}_CHG_methykit.txt"
    out3_file = f"{infile}_CHH_methykit.txt"

    # 打开输入文件和输出文件
    if infile.endswith(".gz"):
        with gzip.open(infile, "rt") as input_file, \
                open(out1_file, "w") as output_file1, \
                open(out2_file, "w") as output_file2, \
                open(out3_file, "w") as output_file3:
            process_lines(input_file, output_file1, output_file2, output_file3)
    else:
        with open(infile, "r") as input_file, \
                open(out1_file, "w") as output_file1, \
                open(out2_file, "w") as output_file2, \
                open(out3_file, "w") as output_file3:
            process_lines(input_file, output_file1, output_file2, output_file3)

def process_lines(input_file, output_file1, output_file2, output_file3):
    # 写入输出文件的表头
    output_file1.write("chrbase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
    output_file2.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
    output_file3.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")

    # 逐行处理输入文件
    for line in input_file:
        line = line.strip()
        columns = line.split("\t")
        chrbase = f"{columns[0]}.{columns[1]}"
        chr_value = columns[0]
        base = columns[1]

        # 根据第三列判断链向
        strand = "R" if "-" in columns[2] else "F"

        # 计算覆盖度和甲基化频率
        coverage = int(columns[3]) + int(columns[4])
        if coverage < 5:
            continue
        freqC = int(columns[3]) / coverage * 100
        freqT = int(columns[4]) / coverage * 100

        # 构建输出行
        output_line = f"{chrbase}\t{chr_value}\t{base}\t{strand}\t{coverage}\t{freqC}\t{freqT}\n"

        # 根据第六列的碱基类型选择写入相应的输出文件
        if "CG" in columns[5]:
            output_file1.write(output_line)
        elif "CHG" in columns[5]:
            output_file2.write(output_line)
        elif "CHH" in columns[5]:
            output_file3.write(output_line)

# 从命令行获取输入文件名
if len(sys.argv) != 2:
    print("请提供输入文件名作为命令行参数。")
    sys.exit(1)
input_file = sys.argv[1]

# 调用函数进行文件转换
convert_bismark_to_methykit(input_file)
