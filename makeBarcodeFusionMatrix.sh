#!/bin/bash

# Making a matrix of barcodes by fusion expression

# Use output from intSite_cellAssociation_scRNAseq.sh
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
