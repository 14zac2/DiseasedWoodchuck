# Pulling out important information from fusions.tsv to compare across samples

# Grab a unique list of genes that integrations occur in
# Note that when an integration is intergenic, this code gets rid of the specific base coordinates but treats the
# two comma separated genes as one gene
cut -f1,2 WHV_WCK_fusions.tsv | sed 's/WhBvgp//g' | 
sed 's/\.//g' | sed 's/\t//g' | sed '/^$/d' | 
sed 's/([^()]*)//g' | tr -d ' ' | 
sort | uniq > integrated_genes.txt

# Collect more extensive information: woodchuck genes involved, woodchuck chromosome, 
# gene structure at integration, integration confidence determined by Arriba, number of supporting reads
# Keeps different integration sites with different coordinates separate from each other
cut -f1,2,5,6,7,8,15,30 WHV_WCK_fusions.tsv | 
awk -v col=8 -F '\t' '{$col=gsub(",", "", $col)+1; print}' | 
sed 's/WhBvgp //g' | sed 's/\. //g' | sed 's/NC_004107[^ ]* //g' | 
sed 's/:[0-9]* / /g' | sed 's/([^()]*)//g' | sed 's/ /\t/g' | 
sed 's/\tCDS\t/\t/g' > int_genes_chroms.tsv
