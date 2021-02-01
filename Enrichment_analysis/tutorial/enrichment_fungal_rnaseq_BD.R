#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = T, warning = F, message = F)

#'# Enrichment analysis (hypergeometric test) for the fungal DEGs

#' This code runs the enrichment_analysis function for the fungal DEGs
#' as a sanity check to see if the function works as expected and gives us 
#' the same results as using any other program for enrichment analysis.

#+
library(tidyverse)

#+
#'# Initial loading and restructuring the data files
#' Set your working directory to the folder with the csv files
# /* setwd("~/OneDrive - University of Central Florida/loop_test/results") */

#' Read the file with Ophio gene annotations
#' [Replace the following file with your organism's gene annotation; ensure the column names GO, PFAM, etc are the same]
all_genes <- read.csv("./FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv", header = T, stringsAsFactors = FALSE)
#+ eval=FALSE
names(all_genes)

#+ 
#' Select only the columns that we need.
#' 
all_genes <- all_genes[,c("arb2_gene","GO","PFAM","signalp","TMHMM")]
head(all_genes)

#' There seem to be duplicate rows for the same gene. Keep only one row per gene.
all_genes <- dplyr::distinct((all_genes))
head(all_genes)


#'# The following code performs GO enrichment only.
#' Repeat the same protocol for any other enrichment (PFAM, signalp, etc)
#'
#' Is there any NAs in the GO column, meaning, are there fungal genes with no GO annotations?
summary(is.na(as.factor(all_genes$GO)))

#' Now, let's replace the NAs in GOs and PFAMs with the term "no_annot"
all_genes[is.na(all_genes)] <- "no_annot"
#' Check for NAs again
summary(is.na(as.factor(all_genes$GO)))
#' There is no more NAs since we replaced them with the category "no annot".

#' Now, let's flatten the file (aim: for each gene, one GO term per line)

all_genes_gos <- all_genes %>%
  dplyr::mutate(go_split = str_split(GO, "; ")) %>% 
  unnest() %>% 
  #dplyr::select(-GO) %>%
  separate(go_split, c("GO","GO_desc"), sep = "([\\|])", extra = "drop") %>% 
  dplyr::select(gene_name = "arb2_gene", GO, GO_desc)

#' Let's take a look to see if it worked.
head(all_genes_gos)

unique(all_genes_gos$gene_name) %>% 
  length()  ## sanity check to see all starting genes have made it to the flattened file

#' How many unique GO terms are there in the Ophio genome?
unique(all_genes_gos$GO) %>% 
  length() 

#+
#'# Now let’s make the function that will do the enrichment analysis

enrichment_analysis <- function(background, geneset) {
  
  # this is going to be your test geneset
  genes <- geneset
  
  # background genes to run enrichment against 
  df.whole <- background # (the flattened 'all_genes_gos' file above)
  
  go_to_desc <- dplyr::distinct(as.data.frame(background[-1]))
  
  # subset the background file to keep only genes in the test set
  df.test <- background %>% 
    filter(gene_name %in% genes)
  
  ## all unique GO terms in the background data set
  annot_terms <- unique(df.whole[[2]])  
  
  ## Create an empty list to save results from the enrichment test
  df.list <- list()
  
  
  for (i in 1:length(annot_terms)) {
    
    # assign the annotation term to test for enrichment
    annot <- annot_terms[i]
    
    # number of genes in the test set
    k <- length(genes)
    
    # number of genes that have the annotation term in the background geneset
    m <- df.whole %>% 
      filter(GO == annot) %>% 
      nrow()
    
    # number of genes that DON'T have the annotation term in the background geneset
    n <- length(unique(df.whole[[1]])) - m
    
    # number of genes that have the annotation term in the test geneset
    x <- df.test %>% 
      filter(GO == annot) %>% 
      nrow()
    
    ### Perform the hypergeometric test
    pval <- dhyper(x, m, n, k, log=F)
    
    # Save the results into the list
    df.list[[i]] <- data.frame(GO = annot_terms[i],
                               x = x,
                               k=k,
                               m=m,
                               n=n,
                               tot_genes = sum(m,n),
                               pVal = pval)
    
    
  }
  
  
  # Row-bind the results for all tested annotation terms
  df.enriched <- bind_rows(df.list, .id = "column_label") %>% 
    select(-column_label) %>% 
    arrange(pVal) %>% 
    ## Adjust the p-values for multiple hypothesis testing (BH-correction)
    mutate(adj_pVal = p.adjust(pVal, "BH")) %>% 
    ## Filter to keep annotation terms that have adjusted p-value < 0.1
    filter(adj_pVal < 0.1) %>%
    ## add a column with the descriptions for each of the GO terms
    left_join(go_to_desc, by="GO")
  
  
  # Return the results as a data.frame
  return(df.enriched);
  
}

#+
#'# Read a test geneset

#' read file
ophio_rnaseq <- read.csv("./enrichtest_fungal.csv", stringsAsFactors = F, header = T)

#' I just need the gene names (first column), let's only keep that.
ophio_rnaseq <- ophio_rnaseq[,1, drop = T]   ## drop = T is necessary for the function to work. 
head(ophio_rnaseq)

#'# Let's run the enrichment analysis
#+
ophio_genes_enrich <- enrichment_analysis(background = all_genes_gos, geneset = ophio_rnaseq)

#' What are top most hits?
knitr::kable(head(ophio_genes_enrich))
#' What about the tail?
knitr::kable(tail(ophio_genes_enrich))

#' Let's look at only the GO terms that have a BH-adjusted pVal < 0.06

enriched_ophio <- ophio_genes_enrich %>% 
  filter(adj_pVal < 0.06)

knitr::kable(enriched_ophio)

#+
#'# Converting the R script to a markdown file.
#'
#' Type the following in the R console:
#' knitr::spin('./tutorial/enrichment_fungal_rnaseq_IW.R')
