#!/bin/bash

for i in {CRR154438,CRR154439,CRR154440,CRR154441,CRR154442}
do
	python wgs.py --fq1 ${i}_f1.fastq.gz --fq2 ${i}_r2.fastq.gz --samplename ${i}
done