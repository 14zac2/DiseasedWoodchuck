#!/bin/bash

### The following code extracts reads from Arriba output "fusions.tsv" and determines which cell they are associated with
### The product is a matrix of barcodes by fusion expression

# Set sample variables
woodchuck="L9647"
tissue="TLH"

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
/data/zoe_fastqs/"$woodchuck"/"$tissue"/*R1* \
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
sed -i '1s/^/barcodes\n/g' fusion_barcodes.tsv

# Add this to filtered_fusions.tsv
paste -d "\t" filtered_fusions.tsv fusion_barcodes.tsv > filtered_fusions_withBarcodes.tsv
# First extract barcodes and associated genes from filtered_fusions_withBarcodes.tsv
cut -f1,2,31 filtered_fusions_withBarcodes.tsv > all_fusion_barcodes.tsv
# Combine both gene columns
awk '{print $1 "," $2 "\t" $3}' all_fusion_barcodes.tsv | sponge all_fusion_barcodes.tsv
# Remove super long parts of gene names
sed -i 's/mikado\.WCK01_AAH20201022_F8-//g' all_fusion_barcodes.tsv
# Remove WhBvgp
sed -e 's/WhBvgp//g' -e 's/,\t/\t/g' -e 's/^,//g' -e 's/,,/,/g' -i all_fusion_barcodes.tsv
# Remove lines that contain periods (as these are empty gene names at this point)
grep -v "\." all_fusion_barcodes.tsv | sponge all_fusion_barcodes.tsv
# Expand the comma separated barcodes, while repeated gene fusion information
sed -Ee 's/^((\S+\t)[^,]+),/\1\n\2/;P;D' -i all_fusion_barcodes.tsv
# Remove header
tail -n +2 all_fusion_barcodes.tsv | sponge all_fusion_barcodes.tsv
# NOTE - could make a contingency table in R at this point, but I'm going to make one in bash
# Count occurances of each fusion
sort all_fusion_barcodes.tsv | uniq -c | sponge all_fusion_barcodes.tsv
# Make contingency table, using 0s instead of NAs for empty values
datamash -sW crosstab 3,2 unique 1 --filler=0 < all_fusion_barcodes.tsv > fusion_countMatrix.tsv
