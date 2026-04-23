#!/usr/bin/env bash

# LIFE748 Assessment 2 - Part 1
# Genome assembly, annotation, and figure generation workflow

# Activate environment
conda activate genomics_env

# Quality control
fastqc GN9_hifix30.fastq -o qc/

# Assembly with Flye
flye --pacbio-hifi GN9_hifix30.fastq --out-dir assemblies/flye --threads 8

# Assembly with SPAdes
spades.py --bio GN9_hifix30.fastq -o assemblies/spades

# QUAST benchmarking
quast.py assemblies/flye/assembly.fasta assemblies/spades/contigs.fasta -o quast/

# Prokka annotation
prokka assemblies/flye/assembly.fasta --outdir annotations/prokka_flye --prefix GN9_flye

# Bakta annotation
bakta --db /path/to/bakta_db --output annotations/bakta_flye assemblies/flye/assembly.fasta

# Extract annotation summaries
python extract_annotation_summary.py

# Generate assembly figures
python make_assembly_figures.py

# Generate annotation figures
python make_annotation_figures.py

# Generate any additional figures
python make_figures.py
