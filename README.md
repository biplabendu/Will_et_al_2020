# Will_et_al_2020
This repository contains annotated code for the function that performs enrichment analysis of GO/PFAM terms using a hypergeometric test. This function was used to perform functional enrichments in Will et al. 2020. 

- Will_et_al_2020
  - Enrichment_analysis
    - tutorial || contains a md file that walks you through the analysis and code
    - results_csv || contains results of the enrichment analysis performed in Will et al. 2020
    - enrichment_fungal_rnaseq_IW.R || R script for customizing the analysis for re-use with other organsisms
    - FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv || gene annotation file for Ophiocodyceps camponoti-floridani used in Will et al. 2020 (replace this with your organism's gene annotation file)
    - enrichtest_fungal.csv || a subset of genes to perform functional enrichment on (replace this with your test geneset)
---------
Citation:

Genetic Underpinnings of Host Manipulation by Ophiocordyceps as Revealed by Comparative Transcriptomics.
I Will, B Das, T Trinh, A Brachmann, RA Ohm, C de Bekker. 
G3: Genes, Genomes, Genetics 10 (7), 2275-2296
