# WHV Analysis
### Analysis of woodchucks infected with the Woodchuck Hepatitis Virus (WHV)

## File list

intSite_cellAssociation_scRNAseq.sh - Process output from Arriba to associate integration sites with specific cells. This requires the following files as input: Arriba output `fusions.tsv`, a TSV file of cell barcodes and their associated cluster from processed scRNA-seq data (`cluster_barcode_IDs.tsv`), and R1 barcodes FASTQ files from CellRanger (`/file/path/to/fastqs/*R1*`).

findIntegrationPatterns.sh - Pulling out important information from fusions.tsv to compare across samples.
