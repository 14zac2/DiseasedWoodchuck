#!/bin/bash

# The purpose of this script is to create fasta of flanking sites around integration sites
# of various sizes to search for motifs
# Uses output from Arriba (fusions.tsv) that has already been filtered to only have WHV reads

# Get all coordinates
echo "Grabbing integration coordinates"
cut -f5-6 WHV_WCK_fusions.tsv > integration_coords.tsv

# Get rid of virus coordinates and remove tabs
# (might be interested in these later for distribution of viral int coords)
echo "Keep and sort woodchuck coordinates"
sed -e 's/\S*NC_\S*//ig' -e 's/\t//g' -i integration_coords.tsv
# Sort coordinates by ascending number
sort -n integration_coords.tsv | uniq | sponge integration_coords.tsv
# Get rid of any blank lines from the file
awk NF integration_coords.tsv | sponge integration_coords.tsv

# STOP HERE. Collect from all samples, then do the rest.
# Concatenate integration coordinates
# All samples begin with L (L*/*/arriba* representing the file path where integration_coords.csv files are stored)
# Place in motif_analysis folder
cat L*/*/arriba*/integration_coords.tsv > ./motif_analysis/integration_coords_allSamples.tsv
# Sort and uniq this file to only keep unique coordinates
sort integration_coords_allSamples.tsv | uniq | sponge integration_coords_allSamples.tsv

# Create temp file used to find flank sites
echo "Creating temporary flank site file"
# Replace colon of coordinate with tab
sed -e 's/:/\t/g' integration_coords_allSamples.tsv > flank_coords.tmp
# Duplicate coordinate column
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' flank_coords.tmp | sponge flank_coords.tmp
# Replace any spaces with tabs
sed -i -e 's/ /\t/g' flank_coords.tmp

# Create different flank site bed files
echo "Creating flanking files"
# Create 10bp flank site bed file
awk -v s=10 'BEGIN { FS="\t"; OFS="\t" } {print $1, $2-s, $3+s}' flank_coords.tmp > flank_coords_10bp.bed
# Create 50bp flank site bed file
awk -v s=50 'BEGIN { FS="\t"; OFS="\t" } {print $1, $2-s, $3+s}' flank_coords.tmp > flank_coords_50bp.bed
# Create 100bp flank site bed file
awk -v s=100 'BEGIN { FS="\t"; OFS="\t" } {print $1, $2-s, $3+s}' flank_coords.tmp > flank_coords_100bp.bed
# Create 100bp flank site bed file
awk -v s=500 'BEGIN { FS="\t"; OFS="\t" } {print $1, $2-s, $3+s}' flank_coords.tmp > flank_coords_500bp.bed
# Create 1000bp flank site bed file
awk -v s=1000 'BEGIN { FS="\t"; OFS="\t" } {print $1, $2-s, $3+s}' flank_coords.tmp > flank_coords_1000bp.bed

# Replace negative values with zero if necessary
sed -i 's/-[0-9][0-9]*/0/' flank_coords_500bp.bed
sed -i 's/-[0-9][0-9]*/0/' flank_coords_1000bp.bed

# Get fasta sequences from files
bedtools getfasta -fi ~/Dropbox/Zoe/scf_version/genome/sc2_ortho_mito_virus.fa \
-bed flank_coords_10bp.bed -fo flank_seq_10bp.fasta
bedtools getfasta -fi ~/Dropbox/Zoe/scf_version/genome/sc2_ortho_mito_virus.fa \
-bed flank_coords_50bp.bed -fo flank_seq_50bp.fasta
bedtools getfasta -fi ~/Dropbox/Zoe/scf_version/genome/sc2_ortho_mito_virus.fa \
-bed flank_coords_100bp.bed -fo flank_seq_100bp.fasta
bedtools getfasta -fi ~/Dropbox/Zoe/scf_version/genome/sc2_ortho_mito_virus.fa \
-bed flank_coords_500bp.bed -fo flank_seq_500bp.fasta
bedtools getfasta -fi ~/Dropbox/Zoe/scf_version/genome/sc2_ortho_mito_virus.fa \
-bed flank_coords_1000bp.bed -fo flank_seq_1000bp.fasta
