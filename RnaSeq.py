"""
RNA-seq analysis pipeline
Author: Jianye Chang
Date: 2023/07/19

Usage: python pipeline.py -1 <fastq_file_1> -2 <fastq_file_2> -g <genome_dir> -f <genome_fasta> -t <gtf_file> -o <output_dir>

Note:
- This pipeline performs RNA-seq analysis using the fastp and STAR tools.
- Make sure that fastp and STAR are properly installed and accessible in your system.
- The provided input files should be in fastq format (.fastq or .fq) and paired-end reads.
- The genome directory should contain the necessary files for STAR indexing.
- The GTF file should contain gene annotations for the reference genome.
- The output directory will be created if it doesn't exist, and the results will be saved there.

Example:
python pipeline.py -1 sample_1.fastq -2 sample_2.fastq -g genome_dir -f genome.fa -t annotations.gtf -o output_results/
"""

import argparse
import logging
import os
import subprocess

def setup_logging(output_dir, samplename):
    """Set up the logging configuration"""
    log_file = os.path.join(output_dir, samplename+'_pipeline.log')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename=log_file)

def run_fastp(fq_1, fq_2, output_dir):
    """Run fastp for quality control and trimming"""
    logging.info('Running fastp...')
    fastp_out = os.path.join(output_dir, fq_1.split('_')[0], 'fastp')
    os.makedirs(fastp_out, exist_ok=True)

    trimmed_fq_1 = os.path.join(fastp_out, f'trimmed_{fq_1}.gz')
    trimmed_fq_2 = os.path.join(fastp_out, f'trimmed_{fq_2}.gz')
    fastp_html = os.path.join(fastp_out, f'{fq_1.rsplit(".", 1)[0]}_fastp.html')
    fastp_json = os.path.join(fastp_out, f'{fq_1.rsplit(".", 1)[0]}_fastp.json')

    command = [
        'fastp',
        '-i', fq_1,
        '-I', fq_2,
        '-o', trimmed_fq_1,
        '-O', trimmed_fq_2,
        '--trim_poly_x',
        '-h', fastp_html,
        '-j', fastp_json,
        '-w', '5'
    ]
    subprocess.run(command, check=True)
    logging.info('fastp completed.')

def run_star(genome_dir, genome_fasta, gtf_file, fq_1, fq_2,output_dir):
    """Run STAR for alignment"""

    star_out = os.path.join(output_dir, fq_1.split('_')[0], 'STAR')
    os.makedirs(star_out, exist_ok=True)

    fq_1_trimmed = os.path.join(output_dir, fq_1.split('_')[0], 'fastp', f'trimmed_{fq_1}.gz')
    fq_2_trimmed = os.path.join(output_dir, fq_1.split('_')[0], 'fastp', f'trimmed_{fq_2}.gz')
    bam_out = os.path.join(star_out, f'{fq_1.rsplit(".", 1)[0]}_sorted.bam')

    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
        logging.info('Running STAR index generated...')
        command = [
            'STAR',
            '--runThreadN', '5',
            '--runMode', 'genomeGenerate',
            '--genomeDir', genome_dir,
            '--genomeFastaFiles', genome_fasta,
            '--sjdbGTFfile', gtf_file,
            '--sjdbOverhang', '149',
            '--genomeSAindexNbases', '13'
        ]
        subprocess.run(command, check=True)
        logging.info('STAR genome index generated successfully.')

    logging.info('Running STAR alignment...')
    command = [
        'STAR',
        '--runThreadN', '5',
        '--genomeDir', genome_dir,
        '--readFilesIn', fq_1_trimmed, fq_2_trimmed,
        '--readFilesCommand', 'zcat',
        '--outFileNamePrefix', os.path.join(star_out, f'{fq_1.rsplit(".", 1)[0]}_'),
        '--outSAMtype', 'BAM', 'Unsorted'
    ]
    subprocess.run(command, check=True)
    logging.info('STAR alignment completed.')

    samtools_sort_command = [
        'samtools',
        'sort', '-@', '5', '-O', 'BAM',
        '-o', bam_out,
        os.path.join(star_out, f'{fq_1.rsplit(".", 1)[0]}_Aligned.out.bam')
    ]
    subprocess.run(samtools_sort_command, check=True)
    logging.info('BAM file sorted.')

def main():
    parser = argparse.ArgumentParser(description='RNA-seq analysis pipeline')
    parser.add_argument('-1', dest='fq_1', required=True, help='Name of fastq file 1')
    parser.add_argument('-2', dest='fq_2', required=True, help='Name of fastq file 2')
    parser.add_argument('-g', dest='genome_dir', required=True, help='Genome directory')
    parser.add_argument('-f', dest='genome_fasta', required=True, help='Genome FASTA file')
    parser.add_argument('-t', dest='gtf_file', required=True, help='GTF file')
    parser.add_argument('-o', dest='output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    setup_logging(args.output_dir, args.fq_1.split('_')[0])
    logging.info('Pipeline started.')
    run_fastp(args.fq_1, args.fq_2, args.output_dir)
    run_star(args.genome_dir, args.genome_fasta, args.gtf_file, args.fq_1, args.fq_2, args.output_dir)
    logging.info('Pipeline completed.')

if __name__ == '__main__':
    main()

