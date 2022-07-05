#!/bin/bash

### The following code extracts reads from Arriba output "fusions.tsv" and determines which cell they are associated with

# Extract desired reads
echo "Extracting reads"
grep "NC_004107" fusions.tsv | grep "WCK01" > WHV_WCK_fusions.tsv
#grep "NC_004107" fusions.discarded.tsv | grep "WCK01" > WHV_WCK_fusions.discarded.tsv
head -n 1 fusions.tsv > header.tsv
cat header.tsv WHV_WCK_fusions.tsv > filtered_fusions.tsv

# Extract read identifiers associated with any fusion, making into newlines and unique
cut -f30 filtered_fusions.tsv | tr , "\n" | tail -n +2 | sort | uniq > fusion_reads.txt

# Extract cell barcodes and UMIs from unique list
echo "Identifying cells with viral integration events"
seqkit grep \
--pattern-file fusion_reads.txt \
/data/zoe_fastqs/L9647/TLH/*R1* \
> R1_fusions.fastq

# Convert barcode fastq file to fasta, then remove fastq files
echo "Converting fastq files"
seqkit fq2fa R1_fusions.fastq -o R1_fusions.fasta
rm R1_fusions.fastq

# Trim to only have cell barcode in fasta file (First 16 bp for scRNA-seq)
fastx_trimmer -f 1 -l 16 -i R1_fusions.fasta -o R1_fusions_cellBarcode.fasta

# Trim to only have read name in fasta header
cut -d ' ' -f1 R1_fusions_cellBarcode.fasta | sponge R1_fusions_cellBarcode.fasta

# Convert fasta to table
seqkit fx2tab R1_fusions_cellBarcode.fasta > R1_fusions_cellBarcode.tsv

# Isolate reads identifiers from WHV_WCK_fusions.tsv and replace them with corresponding cell barcode
cut -f30 WHV_WCK_fusions.tsv > fusion_barcodes.tsv
while IFS=$'\t' read -r header barcode
do
    echo "Replacing" $header "with" ${barcode}-1
    sed -i "s/$header/${barcode}-1/g" fusion_barcodes.tsv
done < R1_fusions_cellBarcode.tsv
# Add "barcodes" header to this file
sed -i 's/read_identifiers/barcodes/g' fusion_barcodes.tsv

# Using R, make a TSV file of cluster numbers and cell IDs
