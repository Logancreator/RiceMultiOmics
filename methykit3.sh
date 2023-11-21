#!/bin/bash
#SBATCH -p smp02,low,big,smp01
#SBATCH -J methykit3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o methykit3.out

# 激活指定的 Conda 环境
source activate /public/home/changjianye/anaconda3/envs/hic/

# 切换到指定的工作目录
cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result

echo "开始运行甲基化分析脚本..."


Rscript methykit3.R -f /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/soapnuke/1F-3/dedup/1F-3_nsort_dedup.deduplicated_sorted.bam -n 1F-3