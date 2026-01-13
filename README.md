# FASTQ → BAM → Variants Pipeline (BWA / Samtools / BCFtools)

This repo contains a simple bash pipeline for processing paired-end FASTQ files through alignment, variant calling, and consensus generation.  
It is designed to run locally on MAC OS for now

## Overview

For each paired-end sample (`*_R1*.fastq` / `*_R2*.fastq`), the pipeline performs:

1. Read alignment with **BWA MEM**
2. BAM conversion, sorting, and indexing (**samtools**)
3. Alignment statistics (`samtools flagstat`)
4. Variant calling (**bcftools mpileup + call**, haploid)
5. Generation of:
   - BCF and indexed VCF
   - CSV summary of variants
   - Consensus FASTA
   - BigWig coverage file (for IGV)


## Requirements

The following tools must be installed and available in `$PATH`:

- **bwa**
- **samtools**
- **bcftools**
- **htslib**
- **deeptools** (for `bamCoverage`)

### Recommended installation (Homebrew)

```bash
brew install bwa samtools bcftools htslib tabix deeptools
```


## Usage

Add info about the pipeline in config.txt:

```text
# Sample name
SAMPLENAME=sample1

# Input directory for FASTQ files
INPUT_DIR=/path/to/fastq/files

# Output directory
OUTPUT_DIR=/path/to/output/directory

# Reference genome FASTA
REFERENCE=/path/to/reference/genome.fasta

# Number of threads
THREADS=4

# Minimum mapping quality
MIN_MAPQ=20

# Minimum base quality
MIN_BASEQ=20

# Coverage threshold for consensus
COVERAGE_THRESHOLD=10
```

Then run:
```bash
bash ./run_pipeline.sh
``` 