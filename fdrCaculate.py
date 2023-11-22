"""
计算Bismark的FDR脚本

作者：[常建业]
日期：[2023/07/14]

将Bismark格式的甲基化测序结果文件计算FDR并筛选甲基化结果。

用法:
    python fdrCaculate.py input.txt output.txt --conversion_rate 99.98

"""

import argparse
import pandas as pd
import csv
from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests
import gzip
import os


def calculate_qv(input_file, output_file, methy_column=2, unmethy_column=3, conversion_rate):
    # 获取文件扩展名
    file_extension = os.path.splitext(input_file)[1]

    if file_extension == '.gz':
        # 读取.gz文件
        with gzip.open(input_file, 'rt') as f:
            # 读取数据文件
            data = pd.read_csv(f, header=None, sep="\t")
    else:
        # 读取普通文本文件
        data = pd.read_csv(input_file, header=None, sep="\t")

    data.columns = ['chr', 'position', 'methy_column', 'unmethy_column', 'C-context', 'background']

    # 打印处理后的总行数
    print("处理前的总行数: ", len(data))

    # 提取甲基化列和非甲基化列
    N1 = data[methy_column].values
    N2 = data[unmethy_column].values

    # 计算 p 值
    pv = binom_test(N1, N1+N2, p=1-(conversion_rate/100), alternative='greater')

    # 将校正后的 q 值添加到数据中
    data['pv'] = pv

    # 对 p 值进行 FDR 校正
    qv = multipletests(pv, method='fdr_bh')[1]

    # 将校正后的 q 值添加到数据中
    data['qv'] = qv

    # 筛选 q 值小于 0.01 的结果且测序深度大于5
    data = data[(data['qv'] <= 0.01) & (data[methy_column]+data[unmethy_column] >= 10)]

    # 输出合并后的结果到文件
    data.to_csv(output_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

    # 打印处理后的总行数
    print("处理后的总行数: ", len(data))

# 命令行解析
parser = argparse.ArgumentParser(description='Calculate q-values')
parser.add_argument('--input_file', type=str, help='Input data file path')
parser.add_argument('--output_file', type=str, help='Output file path')
parser.add_argument('--conversion_rate', type=float, help='Conversion rate')

args = parser.parse_args()

# 调用函数进行运行
calculate_qv(args.input_file, args.output_file, methy_column=2, unmethy_column=3, conversion_rate=args.conversion_rate)

