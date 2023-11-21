#!/bin/bash

cd /public/home/changjianye/goose_WGBS

BatMeth2='/public/home/changjianye/software/BatMeth2-master/bin/BatMeth2'
calmeth='/public/home/changjianye/software/BatMeth2-master/bin/calmeth'
methyGff='/public/home/changjianye/software/BatMeth2-master/bin/methyGff'

$BatMeth2 pipel -1 1C-1_1.fq.gz -2 1C-1_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1C-1 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1C-2_1.fq.gz -2 1C-2_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1C-2 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1C-3_1.fq.gz -2 1C-3_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1C-3 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf

$BatMeth2 pipel -1 1D-1_1.fq.gz -2 1D-1_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1D-1 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1D-2_1.fq.gz -2 1D-2_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1D-2 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1D-3_1.fq.gz -2 1D-3_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1D-3 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf

$BatMeth2 pipel -1 1E-1_1.fq.gz -2 1E-1_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1E-1 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1E-2_1.fq.gz -2 1E-2_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1E-2 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1E-3_1.fq.gz -2 1E-3_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1E-3 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf

$BatMeth2 pipel -1 1F-1_1.fq.gz -2 1F-1_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1F-1 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1F-2_1.fq.gz -2 1F-2_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1F-2 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf
$BatMeth2 pipel -1 1F-3_1.fq.gz -2 1F-3_2.fq.gz -g /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.fa -o 1F-3 -p 20 --gff /public/home/changjianye/goose_WGBS/batmeth3index/AnsCyg.gtf


# #index
# $BatMeth2 index -g lion-head.fa
# #mapping
# for i in {1C-1,1C-2,1C-3,1D-1,1D-2,1D-3,1E-1,1E-2,1E-3,1F-1,1F-2,1F-3}
# do
# 	mkdir /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/${i}
# 	cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/${i}
# 	$BatMeth2 align \
# 	-g /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/batmeth2index/lion-head.fa \
# 	-1 /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/soapnuke/${i}/data_merge/${i}_1.fq.gz \
# 	-2 /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/soapnuke/${i}/data_merge/${i}_2.fq.gz \
# 	-p 48 -o meth_${i}.sam
# 	$calmeth \
# 	-g /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/batmeth2index/lion-head.fa \
# 	-b /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/${i}/meth_${i}.sam.sort.bam \
# 	-Q 20 \
# 	-c 5 \
# 	-m ${i}.methrario.txt
# 	cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/
# done


# #annoatation
# for i in {1C-1,1C-2,1C-3,1D-1,1D-2,1D-3,1E-1,1E-2,1E-3,1F-1,1F-2,1F-3}
# do
# 	cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/${i}
# 	$methyGff \
# 	-B \
# 	-o gene.meth \
# 	-G /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/batmeth2index/lion-head.fa \
# 	-b /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/goose.gene.bed \
# 	-m /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/${i}/${i}.methrario.txt.methratio.txt
# 	cd /public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth2_result/
# done
