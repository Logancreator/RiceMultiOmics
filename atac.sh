#!/bin/bash
#SBATCH -p low,big,smp01,smp02
#SBATCH -J atac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o atac.out

# index creating

cd /public/home/changjianye/project/hzh/Oryza/Ref/bwa

bwa index -a bwtsw -p Oryza Oryza.fna

cd /public/home/changjianye/project/hzh/Oryza/Ref/hisat2

hisat2-build Oryza.fna Oryza

cd /public/home/changjianye/project/hzh/Oryza/Ref/star

STAR  \
--runMode genomeGenerate \
--genomeDir Oryza \
--runThreadN 1 \
--genomeFastaFiles  Oryza.fna \
--sjdbGTFfile Oryza.gtf \
--sjdbOverhang 149

cd /public/home/changjianye/project/hzh/Oryza/Ref/bowtie2

bowtie2-build Oryza.fna Oryza 

cd /public/home/changjianye/project/hzh/ATAC

# pipeline running

python /public/home/changjianye/script/ATAC-seq/AccessibleChromatin.py rice_R1_001.fastq.gz rice_R2_001.fastq.gz