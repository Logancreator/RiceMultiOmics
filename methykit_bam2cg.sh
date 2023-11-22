#!/bin/bash
#SBATCH -p smp02,low,big,smp01
#SBATCH -J methykit_bam2cg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o methykit_bam2cg.out

# 激活指定的 Conda 环境
source activate /public/home/changjianye/anaconda3/envs/hic/

# 切换到指定的工作目录
cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth3_result/tmp

echo "开始运行甲基化分析脚本..."

python methykit_bam2cg.py