#!/bin/bash

### The following code extracts reads from Arriba output "fusions.tsv" and determines which cell they are associated with

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

# Using R, make a TSV file of cluster numbers and cell IDs
# This file is called cluster_barcode_IDs.tsv

# Make a TSV file formatted like fusion_barcodes.tsv, but replace barcodes with clusters
cp fusion_barcodes.tsv fusion_clusters.tsv
while IFS=$'\t' read -r barcode cluster
do
    echo "Replacing" $barcode "with" $cluster
    sed -i "s/$barcode/$cluster/g" fusion_clusters.tsv
done < cluster_barcode_IDs.tsv
# Add "clusters" header to this file
sed -i 's/barcodes/clusters/' fusion_clusters.tsv

# Sometimes there are hundreds or thousands of fusions with some duplicates,
# so it's easiest to assign them a unique ID
# Make sure to delete old labeled_unique_integrations.txt 
rm labeled_unique_integrations.txt
cut -f30 WHV_WCK_fusions.tsv > fusion_reads.tsv
n=0
oldP=0
# Loop through fusion reads
while read p; do
    #echo "Labeling integration event" $n
    if [ "$p" == "$oldP" ]; then
        echo "The integration event is duplicated!"
        echo "Current read" $p
        echo "Last read" $oldP
        #echo -e "Int_${n}" >> labeled_unique_integrations.txt
    else
        #echo "The integration event is unique!"
        #echo -e "Int_${n}" >> labeled_unique_integrations.txt
        n=$((n+1))
    fi
    echo -e "Int_${n}" >> labeled_unique_integrations.txt
    oldP=$p
done < fusion_reads.tsv

# Added header to new file
sed -i '1s/^/intSite\n/g' labeled_unique_integrations.txt

# Add columns to filtered_fusions.tsv
paste -d "\t" labeled_unique_integrations.txt filtered_fusions.tsv fusion_barcodes.tsv fusion_clusters.tsv > filtered_fusions_new.tsv

# Find integration events for Seurat object barcodes
# Also have to delete any old integration_sites_per_cell.txt before running this loop
rm integration_sites_per_cell.txt
while IFS=$'\t' read -r barcode cluster
do
    echo "Looking for barcode" $barcode
    if grep -q "$barcode" filtered_fusions_new.tsv; then
        grep "$barcode" filtered_fusions_new.tsv | awk '{print $1}' | uniq | awk -vORS=, '{print}' | sed 's/,$/\n/' >> integration_sites_per_cell.txt
    else
        echo "NA" >> integration_sites_per_cell.txt
    fi
done < cluster_barcode_IDs.tsv

# Add columns to cluster_barcode_IDs.tsv
paste -d "\t" cluster_barcode_IDs.tsv integration_sites_per_cell.txt > integration_sites_appended_to_cells.tsv

# Looking to try and paste integration site-associated gene names into row beside cell ID
cp integration_sites_per_cell.txt impacted_genes_per_cell.txt
# Loop through integration sites
while read line
do
    intSite=$(awk '{print $1}' <<< "$line")
    gene1=$(awk '{print $2}' <<< "$line")
    gene2=$(awk '{print $3}' <<< "$line")
    genes=$(echo ${gene1},$gene2) # Just in case I don't want to just replace the integration site with gene2
    echo "Replacing $intSite with $genes" ;
    sed -i "s/\b${intSite}\b/$genes/g" impacted_genes_per_cell.txt
done < filtered_fusions_new.tsv
# Remove WhBvgp
sed -e 's/WhBvgp//g' -e 's/,$//g' -e 's/,,/,/g' -i impacted_genes_per_cell.txt

# Append this to per cell info
paste -d "\t" integration_sites_appended_to_cells.tsv impacted_genes_per_cell.txt > int_sites_plus_impacted_genes_per_cell.tsv
# Add header
header=$(echo -e 'barcode\tcluster\tintSite\tgenes')
sed -i -e "1 i\\$header" int_sites_plus_impacted_genes_per_cell.tsv
