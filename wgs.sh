#!/bin/bash
# This script takes two fastq files (paired-end) of WGS data coded by JianyeChang   #2021/03/24
#runs fastp (QC and trimming)
#bwa (mapping trimmed reads)
#gatk (call SNP,Indel,SV)
# USAGE: 
# sh WGSaboutSnpIndel.sh <fq.gz 1 > <fq.gz 2 >
# take consideration of names of fastq files to modify ${samplename}
# take consideration of analyzed species to modify genome
# change ${output_dir} when analyzing different projects

## First To Check about software path (Of course you can module your env to activate brief name.If do it,please modify this script.)

##原始数据路径
input_dir=/public/home/changjianye/project/hzh/rawData

##输出基础路径
base_outputDir=/public/home/changjianye/project/hzh/

##软件路径
fastqc=/public/home/changjianye/anaconda3/envs/env_test/bin/fastqc

fastp=/public/home/changjianye/anaconda3/envs/cuttag/bin/fastp

trim_galore=/public/home/changjianye/anaconda3/envs/env_test/bin/trim_galore

samtools=/public/home/changjianye/anaconda3/envs/cuttag/bin/samtools

bwa=/public/software/env01/bin/bwa

tabix=/public/home/changjianye/anaconda3/envs/cuttag/bin/tabix

bgzip=/public/home/changjianye/anaconda3/envs/cuttag/bin/bgzip

gatk=/public/software/env01/bin/gatk

##CPU数
cores=5

## 打印时间
echo "`date +%Y/%m/%d_%H:%M:%S`"  
## 记录所有
set -x 
START=$(date +%s.%N)

##基因组路径
ref="/public/home/changjianye/project/hzh/HJX74bwaindex/HJX74.fa"

##判断索引文件
if [ ! -f ${ref}.fai ];then
 echo "Genome index not in genome dir ! Now , create it !!!"
 $samtools faidx $ref
 else
 echo "文件存在"
fi

cd $input_dir

fq_1=$1
fq_2=$2

samplename=${fq_1%_1.fq.gz*}
echo "Sample name is $samplename"     

output_dir=${base_outputDir}$samplename/
mkdir -p ${base_outputDir}fastqc
mkdir -p ${output_dir}trim
mkdir -p ${output_dir}mapping
mkdir -p ${output_dir}gatk

fastqc_outdir=${base_outputDir}fastqc/
trim_out=${output_dir}trim/
mapping_outdir=${output_dir}mapping/
gatk_outdir=${output_dir}gatk/

#QC
$fastqc -t $cores -o $fastqc_outdir  $fq_1 $fq_2

#$multiqc  $fastqc_outdir -o $fastqc_outdir -n multiqc_report


##fastp
$fastp \
    -i $fq_1 \
    -o ${fq_1%.fq.gz*}_val_1.fq.gz \
    -I $fq_2 \
    -O ${fq_1%.fq.gz*}_val_2.fq.gz \
    -z 4 \
    -f 5 -t 5 -F 5 -T 5 \
    -5 -W 5 -M 20 \
    -Q \
    -l 50 \
    -c \
    -w 4

#去接头
#$trim_galore -output_dir $trim_out --paired --length 25 --quality 25 --stringency 5 $fq_1 $fq_2 --gzip
clean_fq1=${fq_1%.fq.gz*}_val_1.fq.gz
clean_fq2=${fq_1%.fq.gz*}_val_2.fq.gz


label=@RG\\tID:$samplename\\tPL:MGI\\tSM:$samplename  ##可以修改
##运行BWA
$bwa mem -t $cores -R $label  $ref $clean_fq1 $clean_fq2 | samtools view -Sb - > $mapping_outdir$samplename.bam && echo "** bwa mapping done **"

bam=$mapping_outdir$samplename.bam


##BAM排序
$samtools sort -@ $cores  -O bam -o ${bam%.bam*}.sorted.bam $bam && echo "** BAM sort $samplename"
sorted_bam=${bam%.bam*}.sorted.bam

##标记PCR重复
$gatk MarkDuplicates -I $sorted_bam -O ${bam%.bam*}.sorted.markdup.bam -M ${sorted_bam%.bam*}.markdup_metrics.txt  --REMOVE_DUPLICATES true && echo "** markdup done **"
markdup_bam=${sorted_bam%.bam*}.markdup.bam


##制作index
$samtools index $markdup_bam && echo "** index done **"


##samtools统计信息 
$samtools flagstat $markdup_bam \ 
    > ${markdup_bam%.bam*}.stat


##制作字典
$gatk CreateSequenceDictionary -R $ref -O /public/home/changjianye/project/hzh/HJX74/HJX74.dict && echo "** dict done **"

##gvcf
$gatk  HaplotypeCaller -R $ref  --emit-ref-confidence GVCF -I $markdup_bam  -O ${markdup_bam%.bam*}.g.vcf && echo "** g.vcf done **"
gvcf=${markdup_bam%.bam*}.g.vcf


##vcf
$gatk GenotypeGVCFs \
  -R $ref \
  -V $gvcf \
  -O ${gvcf%.g.vcf*}.vcf && echo "** vcf done **"
vcf=${gvcf%.g.vcf*}.vcf 

##bgzip
$bgzip -f $vcf
vcf_gz=$vcf.gz

##tabix
$tabix  -p vcf $vcf_gz 

# 使用SelectVariants，选出Indel
$gatk SelectVariants \
    -select-type INDEL \
    -V $vcf_gz \
    -O ${vcf_gz%.vcf.gz*}.indel.vcf.gz

# 为Indel作过滤
$gatk VariantFiltration \
    -V ${vcf_gz%.vcf.gz*}.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ${vcf_gz%.vcf.gz*}.indel.filter.vcf.gz
 
# 使用SelectVariants，选出SNP
$gatk SelectVariants \
    -select-type SNP \
    -V $vcf_gz \
    -O ${vcf_gz%.vcf.gz*}.snp.vcf.gz

# 为SNP作硬过滤
$gatk VariantFiltration \
    -V ${vcf_gz%.vcf.gz*}.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ${vcf_gz%.vcf.gz*}.snp.filter.vcf.gz

# 重新合并过滤后的SNP和Indel
$gatk MergeVcfs \
    -I ${vcf_gz%.vcf.gz*}.snp.filter.vcf.gz \
    -I ${vcf_gz%.vcf.gz*}.indel.filter.vcf.gz \
    -O ${vcf_gz%.vcf.gz*}.filter.vcf_includeIndelSNP.gz
    
END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Call SNP Indel completed!"
echo $Duration