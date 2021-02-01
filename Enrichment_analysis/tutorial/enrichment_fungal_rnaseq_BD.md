


# Enrichment analysis (hypergeometric test) for the fungal DEGs
This code runs the enrichment_analysis function for the fungal DEGs
as a sanity check to see if the function works as expected and gives us 
the same results as using any other program for enrichment analysis.


```r
library(goseq)
library(tidyverse)
```


# Initial loading and restructuring the data files
Set your working directory to the folder with the csv files
Read the complete Cflo file with ant annotations and the time course rnaseq data


```r
all_genes <- read.csv("./Enrichment_analysis/Ophio/FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv", header = T, stringsAsFactors = FALSE)
## Error in file(file, "rt"): cannot open the connection
```

```r
names(all_genes)

```


Select only the columns that we need.



```r
all_genes <- all_genes[,c("arb2_gene","GO","PFAM","signalp","TMHMM")]
```

head(all_genes)
There seem to be duplicate rows for the same gene. Keep only one row per gene.


```r
all_genes <- dplyr::distinct((all_genes))
```

head(all_genes)


```r
summary(all_genes)
##   arb2_gene              GO                PFAM          
##  Length:7455        Length:7455        Length:7455       
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##    signalp             TMHMM          
##  Length:7455        Length:7455       
##  Class :character   Class :character  
##  Mode  :character   Mode  :character
```

# GO enrichment only.

Is there any NAs in the GO column, meaning, are there fungal genes with no GO annotations


```r
summary(is.na(as.factor(all_genes$GO)))
##    Mode   FALSE 
## logical    7455
```

Need a filtering step here. 
Basically, which genes should I include in the background? For now, I am using everything.
Now, let's replace the NAs in GOs and pfams with "no_annot"


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
Now, let's flatten the file - for each gene, one GO term per line.


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

How many unique GO terms are there in the Ophio genome


```r
unique(all_genes_gos$GO) %>% 
  length() 
## [1] 2400
```


# Now let’s make the function that will do the enrichment analysis


```r
enrichment_analysis <- function(background, geneset) {
  
  genes <- geneset
  
  df.whole <- background
  
  go_to_desc <- dplyr::distinct(as.data.frame(background[-1]))
  
  df.test <- background %>% 
    filter(gene_name %in% genes)
  annot_terms <- unique(df.whole[[2]])  ## all unique GO terms in the background data set
  
  df.list <- list()
  
  for (i in 1:length(annot_terms)) {
    annot <- annot_terms[i]
    k <- length(genes)
    m <- df.whole %>% 
      filter(GO == annot) %>% 
      nrow()
    n <- length(unique(df.whole[[1]])) - m
    x <- df.test %>% 
      filter(GO == annot) %>% 
      nrow()
    
    pval <- dhyper(x, m, n, k, log=F)
    
    df.list[[i]] <- data.frame(GO = annot_terms[i],
                               x = x,
                               k=k,
                               m=m,
                               n=n,
                               tot_genes = sum(m,n),
                               pVal = pval)
    
    
  }
  
  df.enriched <- bind_rows(df.list, .id = "column_label") %>% 
    select(-column_label) %>% 
    arrange(pVal) %>% 
    mutate(adj_pVal = p.adjust(pVal, "BH")) %>% 
    filter(pVal < 0.1 | adj_pVal < 0.1) %>%
    left_join(go_to_desc, by="GO")
  return(df.enriched);
  
}
```


# Read and format the DEGs from IW's rnaseq data
read the RNASeq file for Ophio DEGs


```r
ophio_rnaseq <- read.csv("./Enrichment_analysis/Ophio/enrichtest_fungal.csv", stringsAsFactors = F, header = T)
## Error in file(file, "rt"): cannot open the connection
```

I just need the gene names, let's only keep that.


```r
ophio_rnaseq <- ophio_rnaseq[,1, drop = T]   ## drop = T is necessary for the function to work. 
## Error in ophio_rnaseq[, 1, drop = T]: incorrect number of dimensions
head(ophio_rnaseq)
## [1] "Ophcf2|00036" "Ophcf2|00037" "Ophcf2|00145" "Ophcf2|00188"
## [5] "Ophcf2|00402" "Ophcf2|00450"
```

# Let's run the enrichment analysis


```r
ophio_genes_enrich <- enrichment_analysis(background = all_genes_gos, geneset = ophio_rnaseq)
```

What are top most hits?


```r
knitr::kable(head(ophio_genes_enrich))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> GO </th>
   <th style="text-align:right;"> x </th>
   <th style="text-align:right;"> k </th>
   <th style="text-align:right;"> m </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> tot_genes </th>
   <th style="text-align:right;"> pVal </th>
   <th style="text-align:right;"> adj_pVal </th>
   <th style="text-align:left;"> GO_desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0006270 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7448 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000001 </td>
   <td style="text-align:right;"> 0.0002643 </td>
   <td style="text-align:left;"> DNA replication initiation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0044421 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular region part </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0090729 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> toxin activity </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0005615 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular space; </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0005615 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular space </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009405 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 7420 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0001027 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> pathogenesis </td>
  </tr>
</tbody>
</table>

What about the tail?


```r
knitr::kable(tail(ophio_genes_enrich))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> GO </th>
   <th style="text-align:right;"> x </th>
   <th style="text-align:right;"> k </th>
   <th style="text-align:right;"> m </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> tot_genes </th>
   <th style="text-align:right;"> pVal </th>
   <th style="text-align:right;"> adj_pVal </th>
   <th style="text-align:left;"> GO_desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 179 </td>
   <td style="text-align:left;"> GO:0019538 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 351 </td>
   <td style="text-align:right;"> 7104 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0931512 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> protein metabolic process </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 180 </td>
   <td style="text-align:left;"> GO:0043603 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 103 </td>
   <td style="text-align:right;"> 7352 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0940360 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> cellular amide metabolic process </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 181 </td>
   <td style="text-align:left;"> GO:0043169 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 348 </td>
   <td style="text-align:right;"> 7107 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0956480 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> cation binding </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 182 </td>
   <td style="text-align:left;"> GO:0043169 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 348 </td>
   <td style="text-align:right;"> 7107 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0956480 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> cation binding; </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 183 </td>
   <td style="text-align:left;"> GO:0046872 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 347 </td>
   <td style="text-align:right;"> 7108 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0964863 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> metal ion binding </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 184 </td>
   <td style="text-align:left;"> GO:0004672 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 101 </td>
   <td style="text-align:right;"> 7354 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0984846 </td>
   <td style="text-align:right;"> 0.9774648 </td>
   <td style="text-align:left;"> protein kinase activity </td>
  </tr>
</tbody>
</table>

Let's look at only the GO terms that have a BH-adjusted pVal < 0.06


```r
enriched_ophio <- ophio_genes_enrich %>% 
  filter(adj_pVal < 0.06)

knitr::kable(enriched_ophio)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> GO </th>
   <th style="text-align:right;"> x </th>
   <th style="text-align:right;"> k </th>
   <th style="text-align:right;"> m </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> tot_genes </th>
   <th style="text-align:right;"> pVal </th>
   <th style="text-align:right;"> adj_pVal </th>
   <th style="text-align:left;"> GO_desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0006270 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7448 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000001 </td>
   <td style="text-align:right;"> 0.0002643 </td>
   <td style="text-align:left;"> DNA replication initiation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0044421 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular region part </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0090729 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> toxin activity </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0005615 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular space; </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0005615 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7421 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0000870 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> extracellular space </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009405 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 7420 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0001027 </td>
   <td style="text-align:right;"> 0.0493157 </td>
   <td style="text-align:left;"> pathogenesis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0044419 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 7417 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0001635 </td>
   <td style="text-align:right;"> 0.0560722 </td>
   <td style="text-align:left;"> interspecies interaction between organisms </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0051704 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 7417 </td>
   <td style="text-align:right;"> 7455 </td>
   <td style="text-align:right;"> 0.0001635 </td>
   <td style="text-align:right;"> 0.0560722 </td>
   <td style="text-align:left;"> multi-organism process </td>
  </tr>
</tbody>
</table>


# Converting the R script to a markdown file.

Type the following in the R console:
knitr::spin('./Enrichment_analysis/Ophio/enrichment_fungal_rnaseq_IW.R')
