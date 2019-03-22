#!/bin/bash
srr="SRR7889584"
fastq-dump.2.9.1 --split-files -O FASTQ  SRR7889584
numLines=$(fastq-dump.2.9.1 --gzip -X 1 -Z --split-spot $srr | wc -l)
if [ $numLines -eq 8 ]; then cat FASTQ/${srr}_1.fastq FASTQ/${srr}_2.fastq > FASTQ/$srr.fastq && rm FASTQ/${srr}_1.fastq FASTQ/${srr}_2.fastq; fi
if [ -f FASTQ/${srr}_1.fastq ]; then mv FASTQ/${srr}_1.fastq FASTQ/${srr}.fastq ; elif [ -f FASTQ/${srr}_2.fastq ]; then mv FASTQ/${srr}_2.fastq FASTQ/${srr}.fastq; fi
