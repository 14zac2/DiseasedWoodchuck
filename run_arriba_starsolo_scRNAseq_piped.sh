#!/bin/bash

mkdir starSolo_whitelist_stringent_output
mkdir arriba_output_solo_stringent
STAR \
--genomeDir /data/zoe_analysis/virus_integration/index_STAR_sc2_ortho_mito_virus_chrom \
--runThreadN 20 \
--readFilesIn /data/zoe_fastqs/L191A/S4/L191A_S4_L001_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L002_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L003_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L004_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L005_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L006_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L007_R2_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L008_R2_001.fastq.gz \
/data/zoe_fastqs/L191A/S4/L191A_S4_L001_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L002_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L003_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L004_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L005_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L006_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L007_R1_001.fastq.gz,\
/data/zoe_fastqs/L191A/S4/L191A_S4_L008_R1_001.fastq.gz \
--soloType CB_UMI_Simple \
--clipAdapterType CellRanger4 \
--outFilterScoreMin 30 \
--soloCBwhitelist /data/zoe_analysis/virus_integration/737K-august-2016.txt \
--outSAMattributes sM CR UR \
--genomeSAsparseD 3 \
--soloCBmatchWLtype Exact \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloCellFilter  None \
--readFilesCommand gunzip -c \
--outFileNamePrefix ./starSolo_whitelist_stringent_output/ \
--outStd BAM_Unsorted \
--genomeLoad NoSharedMemory \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outBAMcompression 0 \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 3 \
--chimSegmentMin 15 \
--chimOutType WithinBAM SoftClip \
--chimJunctionOverhangMin 5 \
--chimScoreMin 1 \
--chimFilter banGenomicN \
--chimScoreDropMax 37 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentReadGapMax 3 |
arriba \
-x /dev/stdin \
-f "blacklist intronic relative_support homologs read_through short_anchor in_vitro intragenic_exonic low_coverage_viral_contigs hairpin end_to_end" \
-g /data/zoe_analysis/sc2_ortho_mito_virus_chrom/genes/genes.gtf \
-a /data/zoe_analysis/sc2_ortho_mito_virus_chrom/fasta/genome.fa \
-o ./arriba_output_solo_stringent/fusions.tsv \
-O ./arriba_output_solo_stringent/fusions.discarded.tsv \
-i "*" -I -T 10 -u -X -S 1 -A 10 -U 32767
