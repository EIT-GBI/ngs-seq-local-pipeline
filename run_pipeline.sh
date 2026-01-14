#!/bin/bash

# Exit immediately if a command doesn't run successfully
set -e

# Load config
source config.txt

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# hanldle case where no files match
shopt -s nullglob

# Find all R1 and R2 files
R1_FILES=("$INPUT_DIR"/*_R1*.fastq)
R2_FILES=("$INPUT_DIR"/*_R2*.fastq)

echo "Found ${#R1_FILES[@]} R1 files."
echo "Found ${#R2_FILES[@]} R2 files."

# Loop through each pair of R1 and R2 files
for R1_FILE in "${R1_FILES[@]}"; do

    BASENAME=$(basename "$R1_FILE")
    SAMPLENAME=${BASENAME%%_R1*}

    R2_FILE="${R1_FILE/_R1/_R2}"

    echo "Processing $R1_FILE"

    # QC block with FastQC
    echo "Running FastQC for $SAMPLENAME"
    mkdir -p "$OUTPUT_DIR/qc"
    fastqc -t "$THREADS" -o "$OUTPUT_DIR/qc" "$R1_FILE" "$R2_FILE"

    # Align reads with BWA 
    echo "Aligning reads for $SAMPLENAME"
    
    # bwa mem
    # -t specifies number of threads
    #
    # samtools view
    # -@ specifies number of threads
    # -bS specifies input is SAM and output is BAM
    bwa mem -t "$THREADS" "$REFERENCE_GENOME" "$R1_FILE" "$R2_FILE" | \
        samtools view -@ "$THREADS" -bS > "$OUTPUT_DIR/${SAMPLENAME}.bam"


    # Sort and index the BAM file
    echo "Sorting and indexing BAM file for $SAMPLENAME"
    samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/$SAMPLENAME.sorted.bam" "$OUTPUT_DIR/$SAMPLENAME.bam"
    samtools index "$OUTPUT_DIR/$SAMPLENAME.sorted.bam"

    # Remove unsorted BAM file to save space
    rm "$OUTPUT_DIR/${SAMPLENAME}.bam"

    ## Stats
    samtools flagstat "$OUTPUT_DIR/$SAMPLENAME.sorted.bam" > "$OUTPUT_DIR/$SAMPLENAME.flagstat.txt"

    ## Variant calling with bcftools
    echo "Calling variants for $SAMPLENAME"
    
    # bcftools mpileup 
    # -O u output uncompressed BCF
    # -f reference genome
    # -q minimum mapping quality
    # -a AD,DP to annotate allele depth and total depth

    # bcftools call
    # --ploidy 1 for haploid genomes
    # -mv for multiallelic calling and variant sites only
    # -Ob output in compressed BCF format
    
    # bcftools filter
    # -e expression to filter variants
    # -Ob output in compressed BCF format

    bcftools mpileup \
        -Ou \
        -f "$REFERENCE_GENOME" \
        -q "$MIN_MAPQ" \
        -a AD,DP \
        "$OUTPUT_DIR/$SAMPLENAME.sorted.bam" \
        | bcftools call \
            --ploidy 1 \
            -mv \
            -Ou \
        | bcftools filter \
            -e "FMT/DP<$MIN_DEPTH || QUAL<$MIN_QUAL" \
            -Ob \
            -o "$OUTPUT_DIR/$SAMPLENAME.calls.bcf"
    
    bcftools index "$OUTPUT_DIR/$SAMPLENAME.calls.bcf"
    

    ## IGV friendly compressed VCF
    # -Oz output in compressed VCF format
    # -p vcf for tabix indexing
    bcftools view -Oz -o "$OUTPUT_DIR/$SAMPLENAME.calls.vcf.gz" "$OUTPUT_DIR/$SAMPLENAME.calls.bcf"
    tabix -p vcf "$OUTPUT_DIR/$SAMPLENAME.calls.vcf.gz"


    bcftools query \
    -f '%CHROM,%POS,%REF,%ALT,%DP,[%AD],%QUAL,%FILTER\n' \
    "$OUTPUT_DIR/$SAMPLENAME.calls.bcf" \
    > "$OUTPUT_DIR/$SAMPLENAME.calls.csv"


    sed -i.bak '1i\
    CHROM,POS,REF,ALT,DP,AD(ref,alt),QUAL,FILTER
    ' "$OUTPUT_DIR/$SAMPLENAME.calls.csv" && rm -f "$OUTPUT_DIR/$SAMPLENAME.calls.csv.bak"


    # Consensus file
    echo "Generating consensus sequence for $SAMPLENAME"
    bcftools consensus -f "$REFERENCE_GENOME" "$OUTPUT_DIR/$SAMPLENAME.calls.bcf" > "$OUTPUT_DIR/$SAMPLENAME.consensus.fna"

    # No more bigWig files, as coverage can be loaded into IGV directly from sorted BAM files
    # Can be re-added if needed

done


