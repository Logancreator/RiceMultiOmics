#!/bin/bash

cd /public/home/changjianye/goose_WGBS

BatMeth2='/public/home/changjianye/software/BatMeth2-master/bin/BatMeth2'
calmeth='/public/home/changjianye/software/BatMeth2-master/bin/calmeth'
methyGff='/public/home/changjianye/software/BatMeth2-master/bin/methyGff'

#$BatMeth2 index -g goose_ncbi.fa

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