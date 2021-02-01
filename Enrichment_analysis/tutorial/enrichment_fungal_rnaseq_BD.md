


# Enrichment analysis (hypergeometric test) for the fungal DEGs
This code runs the enrichment_analysis function for the fungal DEGs
as a sanity check to see if the function works as expected and gives us 
the same results as using any other program for enrichment analysis.


```r
library(tidyverse)
```


# Initial loading and restructuring the data files
Set your working directory to the folder with the csv files
Read the file with Ophio gene annotations
[Replace the following file with your organism's gene annotation; ensure the column names GO, PFAM, etc are the same]


```r
all_genes <- read.csv("./FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv", header = T, stringsAsFactors = FALSE)
```

```r
names(all_genes)

```


Select only the columns that we need.



```r
all_genes <- all_genes[,c("arb2_gene","GO","PFAM","signalp","TMHMM")]
head(all_genes)
##      arb2_gene                                                            GO
## 1 Ophcf2|00001 GO:0003824|catalytic activity; GO:0003674|molecular_function;
## 2 Ophcf2|00001 GO:0003824|catalytic activity; GO:0003674|molecular_function;
## 3 Ophcf2|00001 GO:0003824|catalytic activity; GO:0003674|molecular_function;
## 4 Ophcf2|00002                                                          <NA>
## 5 Ophcf2|00002                                                          <NA>
## 6 Ophcf2|00002                                                          <NA>
##                                PFAM signalp TMHMM
## 1              PF00501|AMP-binding;    <NA>  <NA>
## 2              PF00501|AMP-binding;    <NA>  <NA>
## 3              PF00501|AMP-binding;    <NA>  <NA>
## 4 PF13415|Kelch_3; PF13418|Kelch_4;    <NA>  <NA>
## 5 PF13415|Kelch_3; PF13418|Kelch_4;    <NA>  <NA>
## 6 PF13415|Kelch_3; PF13418|Kelch_4;    <NA>  <NA>
```

There seem to be duplicate rows for the same gene. Keep only one row per gene.


```r
all_genes <- dplyr::distinct((all_genes))
head(all_genes)
##      arb2_gene
## 1 Ophcf2|00001
## 2 Ophcf2|00002
## 3 Ophcf2|00003
## 4 Ophcf2|00004
## 5 Ophcf2|00005
## 6 Ophcf2|00006
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      GO
## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         GO:0003824|catalytic activity; GO:0003674|molecular_function;
## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  <NA>
## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        GO:0003824|catalytic activity; GO:0051920|peroxiredoxin activity; GO:0055114|oxidation-reduction process; GO:0003674|molecular_function; GO:0016491|oxidoreductase activity; GO:0008150|biological_process; GO:0016684|oxidoreductase activity, acting on peroxide as acceptor; GO:0016209|antioxidant activity; GO:0008152|metabolic process;
## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        GO:0016298|lipase activity; GO:0003824|catalytic activity; GO:0016787|hydrolase activity; GO:0052689|carboxylic ester hydrolase activity; GO:0003674|molecular_function; GO:0004806|triglyceride lipase activity; GO:0016788|hydrolase activity, acting on ester bonds; GO:0008150|biological_process; GO:0071704|organic substance metabolic process; GO:0006629|lipid metabolic process; GO:0044238|primary metabolic process; GO:0008152|metabolic process;
## 5 GO:0006694|steroid biosynthetic process; GO:0016614|oxidoreductase activity, acting on CH-OH group of donors; GO:0003674|molecular_function; GO:0071704|organic substance metabolic process; GO:0006629|lipid metabolic process; GO:0003854|3-beta-hydroxy-delta5-steroid dehydrogenase activity; GO:0016229|steroid dehydrogenase activity; GO:0055114|oxidation-reduction process; GO:0003824|catalytic activity; GO:0008610|lipid biosynthetic process; GO:0016491|oxidoreductase activity; GO:0044238|primary metabolic process; GO:0033764|steroid dehydrogenase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor; GO:0008202|steroid metabolic process; GO:0009058|biosynthetic process; GO:1901362|organic cyclic compound biosynthetic process; GO:0008150|biological_process; GO:1901360|organic cyclic compound metabolic process; GO:1901576|organic substance biosynthetic process; GO:0008152|metabolic process; GO:0016616|oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor;
## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      GO:0003824|catalytic activity; GO:0140096|catalytic activity, acting on a protein; GO:0016787|hydrolase activity; GO:0005488|binding; GO:0008270|zinc ion binding; GO:0046872|metal ion binding; GO:0003674|molecular_function; GO:0008233|peptidase activity; GO:0043167|ion binding; GO:0043169|cation binding; GO:0046914|transition metal ion binding; GO:0070011|peptidase activity, acting on L-amino acid peptides; GO:0008237|metallopeptidase activity;
##                                                             PFAM signalp TMHMM
## 1                                           PF00501|AMP-binding;    <NA>  <NA>
## 2                              PF13415|Kelch_3; PF13418|Kelch_4;    <NA>  <NA>
## 3         PF00578|AhpC-TSA; PF08534|Redoxin; PF10417|1-cysPrx_C;    <NA>  <NA>
## 4                              PF01734|Patatin; PF11815|DUF3336;    <NA>  <NA>
## 5       PF05368|NmrA; PF13460|NAD_binding_10; PF01073|3Beta_HSD;    <NA>  <NA>
## 6 PF17900|Peptidase_M1_N; PF01433|Peptidase_M1; PF11838|ERAP1_C;    <NA>  <NA>
```

# The following code performs GO enrichment only.
Repeat the same protocol for any other enrichment (PFAM, signalp, etc)

Is there any NAs in the GO column, meaning, are there fungal genes with no GO annotations?


```r
summary(is.na(as.factor(all_genes$GO)))
##    Mode   FALSE    TRUE 
## logical    3515    3940
```

Now, let's replace the NAs in GOs and PFAMs with the term "no_annot"


```r
all_genes[is.na(all_genes)] <- "no_annot"
```

Check for NAs again


```r
summary(is.na(as.factor(all_genes$GO)))
##    Mode   FALSE 
## logical    7455
```

There is no more NAs since we replaced them with the category "no annot".
Now, let's flatten the file (aim: for each gene, one GO term per line)


```r
all_genes_gos <- all_genes %>%
  dplyr::mutate(go_split = str_split(GO, "; ")) %>% 
  unnest() %>% 
  #dplyr::select(-GO) %>%
  separate(go_split, c("GO","GO_desc"), sep = "([\\|])", extra = "drop") %>% 
  dplyr::select(gene_name = "arb2_gene", GO, GO_desc)
```

Let's take a look to see if it worked.


```r
head(all_genes_gos)
## # A tibble: 6 x 3
##   gene_name    GO         GO_desc                    
##   <chr>        <chr>      <chr>                      
## 1 Ophcf2|00001 GO:0003824 catalytic activity         
## 2 Ophcf2|00001 GO:0003674 molecular_function;        
## 3 Ophcf2|00002 no_annot   <NA>                       
## 4 Ophcf2|00003 GO:0003824 catalytic activity         
## 5 Ophcf2|00003 GO:0051920 peroxiredoxin activity     
## 6 Ophcf2|00003 GO:0055114 oxidation-reduction process

unique(all_genes_gos$gene_name) %>% 
  length()  ## sanity check to see all starting genes have made it to the flattened file
## [1] 7455
```

How many unique GO terms are there in the Ophio genome?


```r
unique(all_genes_gos$GO) %>% 
  length() 
## [1] 2400
```


# Now let’s make the function that will do the enrichment analysis


```r
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
```


# Read a test geneset
read file


```r
ophio_rnaseq <- read.csv("./enrichtest_fungal.csv", stringsAsFactors = F, header = T)
```

I just need the gene names (first column), let's only keep that.


```r
ophio_rnaseq <- ophio_rnaseq[,1, drop = T]   ## drop = T is necessary for the function to work. 
head(ophio_rnaseq)
## [1] "Ophcf2|00036" "Ophcf2|00037" "Ophcf2|00145" "Ophcf2|00188" "Ophcf2|00402" "Ophcf2|00450"
```

# Let's run the enrichment analysis


```r
ophio_genes_enrich <- enrichment_analysis(background = all_genes_gos, geneset = ophio_rnaseq)
```

What are top most hits?


```r
knitr::kable(head(ophio_genes_enrich))
```



|GO         |  x|   k|  m|    n| tot_genes|      pVal|  adj_pVal|GO_desc                    |
|:----------|--:|---:|--:|----:|---------:|---------:|---------:|:--------------------------|
|GO:0006270 |  5| 168|  7| 7448|      7455| 0.0000001| 0.0002643|DNA replication initiation |
|GO:0044421 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular region part  |
|GO:0090729 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|toxin activity             |
|GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space;       |
|GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space        |
|GO:0009405 |  6| 168| 35| 7420|      7455| 0.0001027| 0.0493157|pathogenesis               |

What about the tail?


```r
knitr::kable(tail(ophio_genes_enrich))
```



|   |GO         |  x|   k|  m|    n| tot_genes|      pVal|  adj_pVal|GO_desc                                    |
|:--|:----------|--:|---:|--:|----:|---------:|---------:|---------:|:------------------------------------------|
|3  |GO:0090729 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|toxin activity                             |
|4  |GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space;                       |
|5  |GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space                        |
|6  |GO:0009405 |  6| 168| 35| 7420|      7455| 0.0001027| 0.0493157|pathogenesis                               |
|7  |GO:0044419 |  6| 168| 38| 7417|      7455| 0.0001635| 0.0560722|interspecies interaction between organisms |
|8  |GO:0051704 |  6| 168| 38| 7417|      7455| 0.0001635| 0.0560722|multi-organism process                     |

Let's look at only the GO terms that have a BH-adjusted pVal < 0.06


```r
enriched_ophio <- ophio_genes_enrich %>% 
  filter(adj_pVal < 0.06)

knitr::kable(enriched_ophio)
```



|GO         |  x|   k|  m|    n| tot_genes|      pVal|  adj_pVal|GO_desc                                    |
|:----------|--:|---:|--:|----:|---------:|---------:|---------:|:------------------------------------------|
|GO:0006270 |  5| 168|  7| 7448|      7455| 0.0000001| 0.0002643|DNA replication initiation                 |
|GO:0044421 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular region part                  |
|GO:0090729 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|toxin activity                             |
|GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space;                       |
|GO:0005615 |  6| 168| 34| 7421|      7455| 0.0000870| 0.0493157|extracellular space                        |
|GO:0009405 |  6| 168| 35| 7420|      7455| 0.0001027| 0.0493157|pathogenesis                               |
|GO:0044419 |  6| 168| 38| 7417|      7455| 0.0001635| 0.0560722|interspecies interaction between organisms |
|GO:0051704 |  6| 168| 38| 7417|      7455| 0.0001635| 0.0560722|multi-organism process                     |


# Converting the R script to a markdown file.

Type the following in the R console:
knitr::spin('./tutorial/enrichment_fungal_rnaseq_IW.R')
