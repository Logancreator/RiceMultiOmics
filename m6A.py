import subprocess
import os
import datetime

# Specify the number of cores to use
cores = 5
# Directory with the hisat2 genome index
genome = "/public/home/changjianye/project/duck/hisat2/geese"
bowtie2_index = "/public/home/changjianye/project/duck/rRNA/rRNA"
STAR_index = /public/home/changjianye/project/duck/STAR_index/
# Software path
fastqc = "/public/home/changjianye/miniconda3/envs/atac/bin/fastqc"
fastp = "/public/home/changjianye/miniconda3/envs/atac/bin/fastp"
bowtie2 = "/public/home/changjianye/miniconda3/envs/atac/bin/bowtie2"
samtools = "/public/home/changjianye/miniconda3/envs/atac/bin/samtools"
STAR = "/public/home/changjianye/miniconda3/envs/atac/bin/STAR"
picard = "/public/home/changjianye/miniconda3/envs/atac/bin/picard"
sambamba = "/public/home/changjianye/miniconda3/envs/atac/bin/sambamba"

def run_command(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode('utf-8'), stderr.decode('utf-8')

def get_current_datetime():
    current_datetime = datetime.datetime.now()
    return current_datetime.strftime("%Y/%m/%d_%H:%M:%S")

def main(fq_1, fq_2,genome, cores, output_prefix):

    # Initialize variables
    samplename = os.path.splitext(os.path.basename(fq_1))[0].split('_')[0]
    print("Sample name is", samplename)

    # Make all of the output directories
    output_dir = f"${output_prefix}/{samplename}/"
    print("Output_Dir is", output_dir)

    # Record the date and time
    date_time = get_current_datetime()
    print("Run start ",date_time)

    # Make Dirs
    os.makedirs(output_dir + "trim", exist_ok=True)
    os.makedirs(output_dir + "bam", exist_ok=True)
    os.makedirs(output_dir + "picard", exist_ok=True)
    os.makedirs(output_dir + "peak", exist_ok=True)

    # Set up output directories
    trim_out = output_dir + "trim/"
    bam_out = output_dir + "bam/"
    picard_out = output_dir + "picard/"
    peak_out = output_dir + "peak/"

    # Quality control and read trimming by fastqc
    run_command(f"${fastqc} -t {cores} -o {trim_out} {fq_1} {fq_2}")

    # Quality control and read trimming by fastp
    print("Starting QC and trimming for", samplename)
    run_command(f"${fastp} -i {fq_1} \
        -I {fq_2} \
        -o {trim_out}{samplename}_1_trimmed.fq.gz \
        -O {trim_out}{samplename}_2_trimmed.fq.gz \
        -h {trim_out}{samplename}_fastp.html \
        -j {trim_out}{samplename}_fastp.json \
        --detect_adapter_for_pe \
        --trim_poly_g \
        -l 20 \
        -w {cores}")

    # Quality control and read trimming by fastqc
    run_command(f"${fastqc} -t {cores} -o {trim_out} {trim_out}{samplename}_1_trimmed.fq.gz {trim_out}{samplename}_2_trimmed.fq.gz")

    # Delete the rRNA by bowtie2
    run_command(f"${bowtie2} --very-sensitive-local --no-unal -I 1 -X 1000 -p {cores} \
        -x ${bowtie2_index} \
        -1 {trim_out}{samplename}_1_trimmed.fq.gz \
        -2 {trim_out}{samplename}_2_trimmed.fq.gz \
        --un-conc-gz {bam_out}{samplename}_rRNAremoved.fq.gz 2>{bam_out}{samplename}_Map2rRNAStat.xls | \
        samtools view -S -b -o {bam_out}{samplename}_Map2rRNA.bam -")

    # Quality control and read trimming by fastqc
    run_command(f"${fastqc} -t {cores} -o {output_dir}trim {bam_out}{samplename}_rRNAremoved.fq.1.gz {bam_out}{samplename}_rRNAremoved.fq.2.gz")

    # Map reads to reference genome by STAR
    run_command(f"${STAR} \
        --runThreadN {cores} \
        --runMode alignReads \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --genomeDir ${STAR_index} \
        --readFilesIn {bam_out}{samplename}_rRNAremoved.fq.1.gz {bam_out}{samplename}_rRNAremoved.fq.2.gz \
        --outFileNamePrefix {bam_out}{samplename}_STAR")

    # Index BAM files
    os.chdir(bam_out)
    run_command(f"${samtools} index {bam_out}{samplename}_STARAligned.sortedByCoord.out.bam")

    # Compress sam files to bam files and sort bam file by samtools
    run_command(f"${samtools}  sort -@ {cores} -O bam -o {bam_out}{samplename}_STAR_sorted.bam {bam_out}{samplename}_STARAligned.sortedByCoord.out.bam")

    # Index BAM files
    os.chdir(bam_out)
    run_command(f"${samtools}  index {bam_out}{samplename}_STAR_sorted.bam")

    # Remove duplicates by using picard
    run_command(f"${picard}  MarkDuplicates \
        I={bam_out}{samplename}_STAR_sorted.bam \
        O={picard_out}{samplename}_STAR_sorted_rmDup.bam \
        M={picard_out}{samplename}_STAR_sorted_rmDup_metrics.txt \
        REMOVE_DUPLICATES=true")

    # Filter and keep the uniquely mapped reads
    run_command(f"${sambamba} view -h -t {cores} -f bam -F \"[XS] == null and not unmapped\" \
        {picard_out}{samplename}_STAR_sorted_rmDup.bam > {picard_out}{samplename}_STAR_sorted_rmDup_mapped.bam")

    # Print completion message
    date_time = get_current_datetime()
    print(date_time, "Run completed")
    # 计算过去的日期和时间
    past_datetime = current_datetime - datetime.timedelta(weeks=2)
    print("Run time", past_datetime)

if __name__ == "__main__":
    fq_1 = "path/to/fq_1.fastq.gz"
    fq_2 = "path/to/fq_2.fastq.gz"
    main(fq_1, fq_2)
