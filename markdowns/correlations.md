Correlations
================
Timothy Yu

This notebook generates correlations between different datasets.

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets 
    ##  [9] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.15.4       WGCNA_1.70-3                fastcluster_1.2.3          
    ##  [4] dynamicTreeCut_1.63-1       forcats_0.5.1               stringr_1.4.0              
    ##  [7] purrr_0.3.4                 readr_2.1.1                 tibble_3.1.8               
    ## [10] tidyverse_1.3.1             EnhancedVolcano_1.6.0       ggrepel_0.9.1              
    ## [13] apeglm_1.10.0               GenomicFeatures_1.40.1      AnnotationDbi_1.50.3       
    ## [16] DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
    ## [19] matrixStats_0.61.0          Biobase_2.48.0              GenomicRanges_1.40.0       
    ## [22] GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1           
    ## [25] BiocGenerics_0.34.0         RColorBrewer_1.1-3          viridis_0.6.2              
    ## [28] viridisLite_0.4.0           cowplot_1.1.1               ggplot2_3.3.6              
    ## [31] tidyr_1.2.0                 dplyr_1.0.9                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] circlize_0.4.15          readxl_1.3.1             backports_1.4.0         
    ##   [4] Hmisc_4.6-0              BiocFileCache_1.12.1     systemfonts_1.0.4       
    ##   [7] plyr_1.8.6               splines_4.0.2            BiocParallel_1.22.0     
    ##  [10] digest_0.6.29            invgamma_1.1             foreach_1.5.1           
    ##  [13] htmltools_0.5.2          GO.db_3.11.4             SQUAREM_2021.1          
    ##  [16] fansi_1.0.3              checkmate_2.0.0          magrittr_2.0.3          
    ##  [19] memoise_2.0.1            cluster_2.1.2            doParallel_1.0.16       
    ##  [22] tzdb_0.2.0               Biostrings_2.56.0        annotate_1.66.0         
    ##  [25] modelr_0.1.8             askpass_1.1              bdsmatrix_1.3-4         
    ##  [28] prettyunits_1.1.1        jpeg_0.1-9               colorspace_2.0-3        
    ##  [31] blob_1.2.2               rvest_1.0.2              rappdirs_0.3.3          
    ##  [34] textshaping_0.3.6        haven_2.4.3              xfun_0.28               
    ##  [37] crayon_1.4.2             RCurl_1.98-1.5           jsonlite_1.7.2          
    ##  [40] genefilter_1.70.0        impute_1.62.0            survival_3.3-1          
    ##  [43] iterators_1.0.13         glue_1.6.2               gtable_0.3.0            
    ##  [46] zlibbioc_1.34.0          XVector_0.28.0           GetoptLong_1.0.5        
    ##  [49] shape_1.4.6              scales_1.2.0             mvtnorm_1.1-3           
    ##  [52] DBI_1.1.1                Rcpp_1.0.7               htmlTable_2.3.0         
    ##  [55] xtable_1.8-4             progress_1.2.2           emdbook_1.3.12          
    ##  [58] clue_0.3-60              foreign_0.8-81           bit_4.0.4               
    ##  [61] preprocessCore_1.50.0    Formula_1.2-4            truncnorm_1.0-8         
    ##  [64] htmlwidgets_1.5.4        httr_1.4.2               ellipsis_0.3.2          
    ##  [67] pkgconfig_2.0.3          XML_3.99-0.8             farver_2.1.1            
    ##  [70] nnet_7.3-16              dbplyr_2.1.1             locfit_1.5-9.4          
    ##  [73] utf8_1.2.2               tidyselect_1.1.2         labeling_0.4.2          
    ##  [76] rlang_1.0.4              reshape2_1.4.4           munsell_0.5.0           
    ##  [79] cellranger_1.1.0         tools_4.0.2              cachem_1.0.6            
    ##  [82] cli_3.3.0                generics_0.1.3           RSQLite_2.2.8           
    ##  [85] broom_0.7.10             evaluate_0.14            fastmap_1.1.0           
    ##  [88] yaml_2.2.1               ragg_1.2.2               knitr_1.36              
    ##  [91] bit64_4.0.5              fs_1.5.2                 xml2_1.3.3              
    ##  [94] biomaRt_2.44.4           compiler_4.0.2           rstudioapi_0.13         
    ##  [97] png_0.1-7                curl_4.3.2               reprex_2.0.1            
    ## [100] geneplotter_1.66.0       stringi_1.7.6            highr_0.9               
    ## [103] lattice_0.20-45          Matrix_1.3-4             vctrs_0.4.1             
    ## [106] pillar_1.8.0             lifecycle_1.0.1          GlobalOptions_0.1.2     
    ## [109] data.table_1.14.2        bitops_1.0-7             irlba_2.3.3             
    ## [112] rtracklayer_1.48.0       latticeExtra_0.6-29      R6_2.5.1                
    ## [115] gridExtra_2.3            codetools_0.2-18         MASS_7.3-54             
    ## [118] assertthat_0.2.1         rjson_0.2.21             openssl_1.4.5           
    ## [121] withr_2.5.0              GenomicAlignments_1.24.0 Rsamtools_2.4.0         
    ## [124] GenomeInfoDbData_1.2.3   hms_1.1.1                rpart_4.1-15            
    ## [127] coda_0.19-4              rmarkdown_2.11           ashr_2.2-47             
    ## [130] mixsqp_0.3-43            bbmle_1.0.24             base64enc_0.1-3         
    ## [133] numDeriv_2016.8-1.1      lubridate_1.8.0

``` r
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))
```

Load in datasets.

``` r
# load in liver lipid data
liver_lipids = read.csv('../processed_data/datasets/Lipidomics_liver_normliverweight.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in plasma lipid data
plasma_lipids = read.csv('../processed_data/datasets/Lipidomics_plasma.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in RNA-seq FPKM data
RNAseq = read.csv('../processed_data/datasets/RNAseq_FPKM.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in intestinal bile acids data
intestinal_ba = read.csv('../processed_data/datasets/BA_intestine.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in liver bile acids data
liver_ba = read.csv('../processed_data/datasets/BA_liver.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in biliary bile acids data
biliary_ba = read.csv('../processed_data/datasets/BA_biliary.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in total bile acids data
total_ba = read.csv('../processed_data/datasets/BA_totalpool.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# load in microbiome (rarefied) data
microbiome = read.csv('../processed_data/datasets/Microbiome_otu_rarefied.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
```

``` r
createCorrTable <- function(corr_object, data1_name, data2_name) {
  "
  Create correlation / BH-adjusted p-value table from bicor object.
  
  Args ---
  corr_object : bicor object that contains bicor and p 
  data1_name : data type of first argument in bicorAndPValue()
  data2_name : data type of second argument in bicorAndPValue()
  "
  bicor <- corr_object$bicor
  p <- corr_object$p

  x <- bicor %>% as.data.frame() %>%
  tibble::rownames_to_column(var = data1_name) %>%
  gather(!!data2_name, bicor, 2:ncol(.))
    
  y <- p %>% as.data.frame() %>%
  tibble::rownames_to_column(var = data1_name) %>%
  gather(!!data2_name, pval, 2:ncol(.)) %>%
  mutate(adj_pval = p.adjust(pval, method = 'BH', n = length(pval)))

  table <- left_join(x, y, by = c(data1_name, data2_name)) %>%
    dplyr::select(!!data1_name, !!data2_name, bicor, adj_pval)
  return(table)
}
```

Correlate liver lipids -> RNA-seq

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(liver_lipids$ID, RNAseq$ID)
liver_lipids_mat = liver_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
RNAseq_mat = RNAseq %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, RNAseq_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((RNAseq_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
gene_x_liver_lipid = createCorrTable(corr, 'Gene', 'Liver_lipid') %>%
  filter(adj_pval < 0.05)
write.csv(gene_x_liver_lipid, '../processed_data/correlations/gene_x_liver_lipid_adjpval0.05.csv', row.names = FALSE)

top = (gene_x_liver_lipid %>% 
  arrange(adj_pval) %>% 
  dplyr::select(Gene) %>% 
  distinct() %>% 
  head(150))$Gene
```

Correlate liver lipids -> top 150 RNAseq

``` r
top_RNAseq_mat = RNAseq_mat %>% dplyr::select(ID, Treatment, all_of(top))

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, top_RNAseq_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((top_RNAseq_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
topgene_x_liver_lipid = createCorrTable(corr, 'Gene', 'Liver_lipid')
write.csv(topgene_x_liver_lipid, '../processed_data/correlations/top150gene_x_liver_lipid.csv', row.names = FALSE)
rm(gene_x_liver_lipid, liver_lipids_mat, RNAseq_mat, top_RNAseq_mat, topgene_x_liver_lipid)
```

Correlate liver lipids -> intestinal BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(liver_lipids$ID, intestinal_ba$ID)
liver_lipids_mat = liver_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
intestinal_ba_mat = intestinal_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, intestinal_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (intestinal_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
liver_lipid_x_intestinal_ba = createCorrTable(corr, 'Liver_lipid', 'Intestinal_ba')
write.csv(liver_lipid_x_intestinal_ba, '../processed_data/correlations/liver_lipid_x_intestinal_ba.csv',
          row.names = FALSE)
rm(liver_lipid_x_intestinal_ba, liver_lipids_mat, intestinal_ba_mat)
```

Correlate liver lipids -> liver BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(liver_lipids$ID, liver_ba$ID)
liver_lipids_mat = liver_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
liver_ba_mat = liver_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, liver_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (liver_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
liver_lipid_x_liver_ba = createCorrTable(corr, 'Liver_lipid', 'Liver_ba')
write.csv(liver_lipid_x_liver_ba, '../processed_data/correlations/liver_lipid_x_liver_ba.csv',
          row.names = FALSE)
rm(liver_lipid_x_liver_ba, liver_lipids_mat, liver_ba_mat)
```

Correlate liver lipids -> biliary BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(liver_lipids$ID, biliary_ba$ID)
liver_lipids_mat = liver_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
biliary_ba_mat = biliary_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, biliary_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (biliary_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
liver_lipid_x_biliary_ba = createCorrTable(corr, 'Liver_lipid', 'Biliary_ba')
write.csv(liver_lipid_x_biliary_ba, '../processed_data/correlations/liver_lipid_x_biliary_ba.csv',
          row.names = FALSE)
rm(liver_lipid_x_biliary_ba, liver_lipids_mat, biliary_ba_mat)
```

Correlate liver lipids -> microbiome

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(liver_lipids$ID, microbiome$ID)
liver_lipids_mat = liver_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(liver_lipids_mat$ID, microbiome_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((liver_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
liver_lipid_x_microbiome = createCorrTable(corr, 'Liver_lipid', 'Microbiome_otu')
write.csv(liver_lipid_x_microbiome, '../processed_data/correlations/liver_lipid_x_microbiome.csv',
          row.names = FALSE)
rm(liver_lipid_x_microbiome, liver_lipids_mat, microbiome_mat)
```

Correlate plasma lipids -> RNA-seq

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(plasma_lipids$ID, RNAseq$ID)
plasma_lipids_mat = plasma_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
RNAseq_mat = RNAseq %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(plasma_lipids_mat$ID, RNAseq_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((RNAseq_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (plasma_lipids_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
gene_x_plasma_lipid = createCorrTable(corr, 'Gene', 'Plasma_lipid') %>%
  filter(adj_pval < 0.05)
write.csv(gene_x_plasma_lipid, '../processed_data/correlations/gene_x_plasma_lipid_adjpval0.05.csv', 
          row.names = FALSE)
rm(gene_x_plasma_lipid, plasma_lipids_mat, RNAseq_mat)
```

Correlate plasma lipids -> intestinal BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(plasma_lipids$ID, intestinal_ba$ID)
plasma_lipids_mat = plasma_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
intestinal_ba_mat = intestinal_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(plasma_lipids_mat$ID, intestinal_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((plasma_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (intestinal_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
plasma_lipid_x_intestinal_ba = createCorrTable(corr, 'Plasma_lipid', 'Intestinal_ba')
write.csv(plasma_lipid_x_intestinal_ba,'../processed_data/correlations/plasma_lipid_x_intestinal_ba.csv',
          row.names = FALSE)
rm(plasma_lipid_x_intestinal_ba, plasma_lipids_mat, intestinal_ba_mat)
```

Correlate plasma lipids -> liver BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(plasma_lipids$ID, liver_ba$ID)
plasma_lipids_mat = plasma_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
liver_ba_mat = liver_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(plasma_lipids_mat$ID, liver_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((plasma_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (liver_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
plasma_lipid_x_liver_ba = createCorrTable(corr, 'Plasma_lipid', 'Liver_ba')
write.csv(plasma_lipid_x_liver_ba,'../processed_data/correlations/plasma_lipid_x_liver_ba.csv',
          row.names = FALSE)
rm(plasma_lipid_x_liver_ba, plasma_lipids_mat, liver_ba_mat)
```

Correlate plasma lipids -> biliary BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(plasma_lipids$ID, biliary_ba$ID)
plasma_lipids_mat = plasma_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
biliary_ba_mat = biliary_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(plasma_lipids_mat$ID, biliary_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((plasma_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (biliary_ba_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
plasma_lipid_x_biliary_ba = createCorrTable(corr, 'Plasma_lipid', 'Biliary_ba')
write.csv(plasma_lipid_x_biliary_ba,'../processed_data/correlations/plasma_lipid_x_biliary_ba.csv',
          row.names = FALSE)
rm(plasma_lipid_x_biliary_ba, plasma_lipids_mat, biliary_ba_mat)
```

Correlate plasma lipids -> microbiome

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(plasma_lipids$ID, microbiome$ID)
plasma_lipids_mat = plasma_lipids %>% filter(ID %in% common_ids) %>% arrange(ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(plasma_lipids_mat$ID, microbiome_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((plasma_lipids_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
plasma_lipid_x_microbiome = createCorrTable(corr, 'Plasma_lipid', 'Microbiome_otu')
write.csv(plasma_lipid_x_microbiome,'../processed_data/correlations/plasma_lipid_x_microbiome.csv',
          row.names = FALSE)
rm(plasma_lipid_x_microbiome, plasma_lipids_mat, microbiome_mat)
```

Correlate microbiome -> RNA-seq

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(microbiome$ID, RNAseq$ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)
RNAseq_mat = RNAseq %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(microbiome_mat$ID, RNAseq_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((RNAseq_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
gene_x_microbiome = createCorrTable(corr, 'Gene', 'Microbiome_otu') %>%
  filter(adj_pval < 0.05)
write.csv(gene_x_microbiome, '../processed_data/correlations/gene_x_microbiome_adjpval0.05.csv', row.names = FALSE)
rm(gene_x_microbiome, microbiome_mat, RNAseq_mat)
```

Correlate microbiome -> intestinal BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(microbiome$ID, intestinal_ba$ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)
intestinal_ba_mat = intestinal_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(microbiome_mat$ID, intestinal_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((intestinal_ba_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
intestinal_ba_x_microbiome = createCorrTable(corr, 'Intestinal_ba', 'Microbiome_otu')
write.csv(intestinal_ba_x_microbiome, '../processed_data/correlations/intestinal_ba_x_microbiome.csv',
          row.names = FALSE)
rm(intestinal_ba_x_microbiome, microbiome_mat, intestinal_ba_mat)
```

Correlate microbiome -> liver BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(microbiome$ID, liver_ba$ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)
liver_ba_mat = liver_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(microbiome_mat$ID, liver_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((liver_ba_mat %>% dplyr::select(-c('ID', 'Treatment'))), 
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
liver_ba_x_microbiome = createCorrTable(corr, 'Liver_ba', 'Microbiome_otu')
write.csv(liver_ba_x_microbiome, '../processed_data/correlations/liver_ba_x_microbiome.csv',
          row.names = FALSE)
rm(liver_ba_x_microbiome, microbiome_mat, liver_ba_mat)
```

Correlate microbiome -> biliary BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(microbiome$ID, biliary_ba$ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)
biliary_ba_mat = biliary_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(microbiome_mat$ID, biliary_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((biliary_ba_mat %>% dplyr::select(-c('ID', 'Treatment'))),
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
biliary_ba_x_microbiome = createCorrTable(corr, 'Biliary_ba', 'Microbiome_otu')
write.csv(biliary_ba_x_microbiome, '../processed_data/correlations/biliary_ba_x_microbiome.csv',
          row.names = FALSE)
rm(biliary_ba_x_microbiome, microbiome_mat, biliary_ba_mat)
```

Correlate microbiome -> total BA

``` r
# only correlate if ID's are in both data sets
common_ids = intersect(microbiome$ID, total_ba$ID)
microbiome_mat = microbiome %>% filter(ID %in% common_ids) %>% arrange(ID)
total_ba_mat = total_ba %>% filter(ID %in% common_ids) %>% arrange(ID)

# check that order of ID's are the same in both data sets
assertthat::are_equal(microbiome_mat$ID, total_ba_mat$ID)
```

    ## [1] TRUE

``` r
# correlate (NA's and missing values exist, so warnings expected)
corr = bicorAndPvalue((total_ba_mat %>% dplyr::select(-c('ID', 'Treatment'))),
                      (microbiome_mat %>% dplyr::select(-c('ID', 'Treatment')))
       )
```

``` r
total_ba_x_microbiome = createCorrTable(corr, 'Total_ba', 'Microbiome_otu')
write.csv(total_ba_x_microbiome, '../processed_data/correlations/total_ba_x_microbiome.csv',
          row.names = FALSE)
rm(total_ba_x_microbiome, microbiome_mat, total_ba_mat)
```

``` r
# rmarkdown::render("correlations.Rmd")
# mv correlations.md ../markdowns/
```
