#!/bin/bash

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything
START=$(date +%s.%N)


# # initialize a variable with an intuitive name to store the name of the input fastq file
fq_1=$1
fq_2=$2

# grab base of filename for naming outputs
samplename=${fq_1%.R1.fastq.gz*}
echo "Sample name is $samplename"     

# specify the number of cores to use
cores=10

# directory with the bowtie2 genome index
genome=/disk1/liujx/database_genome/genome/NCBI/Sichuan_goose_rename_addN/genome_addN


# make all of the output directories
# The -p option means no error if existing, make parent directories as needed
output_dir=/disk1/liujx/project/cuttag_anser_4stage_pituirary/cjy_cuttag/$samplename/
mkdir -p ${output_dir}fastp
mkdir -p ${output_dir}bam
mkdir -p ${output_dir}bigwig
mkdir -p ${output_dir}MACS2

# set up output directories
fastp_out=${output_dir}fastp/
bam_out=${output_dir}bam/
peak_out=${output_dir}MACS2/
bigwig_out=${output_dir}bigwig/



## Quality control and read trimming by fastp
echo "Starting QC and trimming for $samplename"
/home/liujx/miniconda3/envs/cuttag/bin/fastp -i ${fq_1} \
	    -I ${fq_2} \
		-o ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
		-O ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
		-h ${fastp_out}${samplename}_fastp.html \
		-j ${fastp_out}${samplename}_fastp.json \
		-w ${cores}


## Map reads to reference genome by bowtie2
echo "Starting mapping for $samplename"
/home/liujx/miniconda3/envs/cuttag/bin/bowtie2 \
	--end-to-end --very-sensitive \
	--no-mixed --no-discordant \
	--phred33 -I 10 -X 700  -p ${cores} \
	-x ${genome} \
	-1 ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-2 ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	-S ${bam_out}${samplename}.sam &> ${bam_out}${samplename}_bowtie2_summary.txt

## Compress sam files to bam files and sort bam file by samtools
/home/liujx/miniconda3/envs/cuttag/bin/samtools view -S -b ${bam_out}${samplename}.sam \
	| /home/liujx/miniconda3/envs/cuttag/bin/samtools sort -@ ${cores} -O bam -o ${bam_out}${samplename}_sorted.bam


## Remove duplicates by using picard
echo "Remove duplicates of $samplename"
/home/liujx/miniconda3/envs/cuttag/bin/picard MarkDuplicates \
	I=${bam_out}${samplename}_sorted.bam \
	O=${bam_out}${samplename}_sorted_rmDup.bam \
	M=${bam_out}${samplename}_rmDup_metrics.txt \
	REMOVE_DUPLICATES=true

## Sorting BAM files by genomic coordinates
mv ${bam_out}${samplename}_sorted_rmDup.bam ${bam_out}${samplename}_rmDup.bam 
/home/liujx/miniconda3/envs/cuttag/bin/samtools sort -@ ${cores} -O bam -o ${bam_out}${samplename}_sorted_rmDup.bam ${bam_out}${samplename}_rmDup.bam


## Filter and keep the uniquely mapped reads
# filtered out multimappers by specifying `[XS] == null`
# filtered out unmapped reads by specifying in the filter `not unmapped`
echo "Keep the uniquely mapped reads of $samplename"
/home/liujx/miniconda3/envs/cuttag/bin/sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped" \
	${bam_out}${samplename}_sorted_rmDup.bam > ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam

## index BAM files
cd ${bam_out}
/home/liujx/miniconda3/envs/cuttag/bin/samtools index ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam

## bigwig
/home/liujx/miniconda3/envs/cuttag/bin/bamCoverage --bam ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam \
	--numberOfProcessors ${cores} --binSize 10 \
	--normalizeUsing RPKM \
	-o ${bigwig_out}${samplename}_RPKM_normalized.bw


END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration


