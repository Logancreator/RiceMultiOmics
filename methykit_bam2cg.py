import subprocess

# 样本编号列表
samples = ['1C-1', '1C-2', '1C-3', '1D-1', '1D-2', '1D-3', '1E-1', '1E-2', '1E-3', '1F-1', '1F-2', '1F-3']

# BAM文件所在目录
bam_dir = "/public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth3_result/tmp"

# 循环处理每个样本编号
for sample in samples:
    # 构造BAM文件路径
    bam_file = os.path.join(bam_dir, f"{sample}.sort.bam")
    # 构造命令列表
    cmd = ["Rscript", "methykit_bam2cg.R", "-f", bam_file, "-n", sample]
    # 打印信息
    print(f"Processing BAM file: {bam_file} for sample {sample}")
    print(f"Processing command : {cmd}")
    # 执行命令
    subprocess.run(cmd)
    print(f"Processing BAM file: {bam_file} for sample {sample} done !!! ")