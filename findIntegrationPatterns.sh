# Pulling out important information from fusions.tsv to compare across samples
# Note that it uses the output WHV_WCK_fusions.tsv created in intSite_cellAssociation_scRNAseq.sh
# But here is the prep code for convenience
grep "NC_004107" fusions.tsv | grep "WCK01" > WHV_WCK_fusions.tsv
head -n 1 fusions.tsv > header.tsv
cat header.tsv WHV_WCK_fusions.tsv > filtered_fusions.tsv

# Grab a unique list of genes that integrations occur in
# Note that when an integration is intergenic, this code gets rid of the specific base coordinates but treats the
# two comma separated genes as one gene
cut -f1,2 WHV_WCK_fusions.tsv | sed 's/WhBvgp//g' | 
sed 's/\.//g' | sed 's/\t//g' | sed '/^$/d' | 
sed 's/([^()]*)//g' | tr -d ' ' | 
sort | uniq > integrated_genes.tsv
# Add sample ID and tissue type
sampleID="L9647_TLH"
tissue="Infected"
awk -v tissue="$tissue" -v sampleID="$sampleID" -F'\t' 'BEGIN {OFS = FS} {print sampleID, $0, tissue}' integrated_genes.tsv | 
sponge integrated_genes.tsv

# Collect more extensive information: woodchuck genes involved (with coordinate), woodchuck chromosome, 
# gene structure at integration, integration confidence determined by Arriba, number of supporting reads
# Keeps different integration sites with different coordinates separate from each other
cut -f1,2,5,6,7,8,15,30 WHV_WCK_fusions.tsv | 
awk -v col=8 -F '\t' '{$col=gsub(",", "", $col)+1; print}' | 
sed 's/WhBvgp //g' | sed 's/\. //g' | sed 's/NC_004107[^ ]* //g' | 
sed 's/:[0-9]* / /g' | sed 's/ /\t/g' | 
sed 's/\tCDS\t/\t/g' > int_genes_chroms.tsv
# Add sample ID and tissue type
sampleID="L9647_TLH"
tissue="Infected"
awk -v tissue="$tissue" -v sampleID="$sampleID" -F'\t' 'BEGIN {OFS = FS} {print sampleID, $0, tissue}' int_genes_chroms.tsv | 
sponge int_genes_chroms.tsv
