# WHV Analysis
### Analysis of woodchucks infected with the Woodchuck Hepatitis Virus (WHV)

## File list

run_arriba_starsolo_scRNAseq_piped.sh - Run STARSolo and Arriba on scRNA-seq FASTQ files to find integration sites. Makes two directories: `starSolo_whitelist_stringent_output` and `arriba_output_solo_stringent` which contain STAR and Arriba outputs respectively.

intSite_cellAssociation_scRNAseq.sh - Process output from Arriba to associate integration sites with specific cells. This requires the following files as input: Arriba output `fusions.tsv`, a TSV file of cell barcodes and their associated cluster from processed scRNA-seq data (`cluster_barcode_IDs.tsv`), and R1 barcodes FASTQ files from CellRanger (`/file/path/to/fastqs/*R1*`).

findIntegrationPatterns.sh - Pulling out important information from fusions.tsv to compare across samples.

makeBarcodeFusionMatrix.sh - Make a matrix of the number of fusion events for each barcode based on Arriba output.
