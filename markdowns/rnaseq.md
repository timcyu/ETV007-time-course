Analysis of RNA-seq data
================
Timothy Yu

This notebook analyzes the RNA-seq dataset.

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
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] EnhancedVolcano_1.6.0       ggrepel_0.9.1               apeglm_1.10.0              
    ##  [4] GenomicFeatures_1.40.1      AnnotationDbi_1.50.3        DESeq2_1.28.1              
    ##  [7] SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.61.0         
    ## [10] Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
    ## [13] IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0        
    ## [16] WGCNA_1.70-3                fastcluster_1.2.3           dynamicTreeCut_1.63-1      
    ## [19] rmarkdown_2.11              cowplot_1.1.1               RColorBrewer_1.1-3         
    ## [22] DescTools_0.99.44           viridis_0.6.2               viridisLite_0.4.0          
    ## [25] bnstruct_1.0.11             igraph_1.2.9                bitops_1.0-7               
    ## [28] ggfortify_0.4.13            forcats_0.5.1               stringr_1.4.0              
    ## [31] dplyr_1.0.9                 purrr_0.3.4                 readr_2.1.1                
    ## [34] tidyr_1.2.0                 tibble_3.1.8                ggplot2_3.3.6              
    ## [37] tidyverse_1.3.1            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2               tidyselect_1.1.2         RSQLite_2.2.8           
    ##   [4] htmlwidgets_1.5.4        grid_4.0.2               BiocParallel_1.22.0     
    ##   [7] munsell_0.5.0            codetools_0.2-18         ragg_1.2.2              
    ##  [10] preprocessCore_1.50.0    withr_2.5.0              colorspace_2.0-3        
    ##  [13] highr_0.9                knitr_1.36               rstudioapi_0.13         
    ##  [16] labeling_0.4.2           bbmle_1.0.24             GenomeInfoDbData_1.2.3  
    ##  [19] mixsqp_0.3-43            bit64_4.0.5              farver_2.1.1            
    ##  [22] coda_0.19-4              vctrs_0.4.1              generics_0.1.3          
    ##  [25] xfun_0.28                BiocFileCache_1.12.1     R6_2.5.1                
    ##  [28] doParallel_1.0.16        invgamma_1.1             locfit_1.5-9.4          
    ##  [31] cachem_1.0.6             assertthat_0.2.1         scales_1.2.0            
    ##  [34] nnet_7.3-16              rootSolve_1.8.2.3        gtable_0.3.0            
    ##  [37] lmom_2.8                 rlang_1.0.4              genefilter_1.70.0       
    ##  [40] systemfonts_1.0.4        splines_4.0.2            rtracklayer_1.48.0      
    ##  [43] impute_1.62.0            broom_0.7.10             checkmate_2.0.0         
    ##  [46] yaml_2.2.1               reshape2_1.4.4           modelr_0.1.8            
    ##  [49] backports_1.4.0          Hmisc_4.6-0              tools_4.0.2             
    ##  [52] ellipsis_0.3.2           proxy_0.4-26             Rcpp_1.0.7              
    ##  [55] plyr_1.8.6               base64enc_0.1-3          progress_1.2.2          
    ##  [58] zlibbioc_1.34.0          RCurl_1.98-1.5           prettyunits_1.1.1       
    ##  [61] rpart_4.1-15             openssl_1.4.5            ashr_2.2-47             
    ##  [64] haven_2.4.3              cluster_2.1.2            fs_1.5.2                
    ##  [67] magrittr_2.0.3           data.table_1.14.2        reprex_2.0.1            
    ##  [70] truncnorm_1.0-8          mvtnorm_1.1-3            SQUAREM_2021.1          
    ##  [73] hms_1.1.1                evaluate_0.14            xtable_1.8-4            
    ##  [76] XML_3.99-0.8             emdbook_1.3.12           jpeg_0.1-9              
    ##  [79] readxl_1.3.1             gridExtra_2.3            compiler_4.0.2          
    ##  [82] biomaRt_2.44.4           bdsmatrix_1.3-4          crayon_1.4.2            
    ##  [85] htmltools_0.5.2          mgcv_1.8-38              tzdb_0.2.0              
    ##  [88] Formula_1.2-4            geneplotter_1.66.0       expm_0.999-6            
    ##  [91] Exact_3.1                lubridate_1.8.0          DBI_1.1.1               
    ##  [94] dbplyr_2.1.1             MASS_7.3-54              rappdirs_0.3.3          
    ##  [97] boot_1.3-28              Matrix_1.3-4             cli_3.3.0               
    ## [100] pkgconfig_2.0.3          GenomicAlignments_1.24.0 numDeriv_2016.8-1.1     
    ## [103] foreign_0.8-81           xml2_1.3.3               foreach_1.5.1           
    ## [106] annotate_1.66.0          XVector_0.28.0           rvest_1.0.2             
    ## [109] digest_0.6.29            Biostrings_2.56.0        cellranger_1.1.0        
    ## [112] htmlTable_2.3.0          gld_2.6.3                curl_4.3.2              
    ## [115] Rsamtools_2.4.0          lifecycle_1.0.1          nlme_3.1-153            
    ## [118] jsonlite_1.7.2           askpass_1.1              fansi_1.0.3             
    ## [121] pillar_1.8.0             lattice_0.20-45          fastmap_1.1.0           
    ## [124] httr_1.4.2               survival_3.3-1           GO.db_3.11.4            
    ## [127] glue_1.6.2               png_0.1-7                iterators_1.0.13        
    ## [130] bit_4.0.4                class_7.3-19             stringi_1.7.6           
    ## [133] blob_1.2.2               textshaping_0.3.6        latticeExtra_0.6-29     
    ## [136] memoise_2.0.1            irlba_2.3.3              e1071_1.7-9

``` r
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))
```

Load in all sample files. Note: one outlier `7B_Index_09_control8` was
removed after visualizing PCA. If want to visualize it without that
sample, do not run line of code where it removes file name from
`sampleFiles`.

``` r
directory <- "../raw_data/ETV007_RNAseq_counts"
sampleFiles <- (list.files(directory))
sampleFiles <- sampleFiles[-33]

# load in the gene key information
key = read.table('../raw_data/Mappings/ensembl84_info.txt', header = TRUE, sep = '\t', quote="", fill=FALSE, stringsAsFactors = FALSE)
```

Get gene lengths for downstream FPKM normalization.

``` r
gene_info = key %>% mutate(gene_length = Gene.End..bp. - Gene.Start..bp. + 1) %>% dplyr::select(Ensembl.Gene.ID, Associated.Gene.Name, gene_length) %>% distinct()
```

Create `DESeq` object.

``` r
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleFiles)

sampleTable <- sampleTable %>% 
  separate(col = 'condition', into = c("exp", "index", "condition"), sep = "_", remove = F)
sampleTable <- subset(sampleTable, select = -c(exp, index))
sampleTable$condition <- gsub(".count", "", sampleTable$condition)
sampleTable$condition <- gsub("control8", "0", sampleTable$condition)
sampleTable$condition <- gsub("treated", "", sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~condition)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula
    ## are characters, converting to factors

``` r
#Set the reference group to 0 week
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "0")

#Filter out reads whose row sums are less than 10
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=10,]

# Run DESeq
dds <- DESeq(ddsHTSeq)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 7 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

PCA of variance-stabilized transformed data.

``` r
vsd = vst(dds)
levels(vsd$condition) = c("CHOW8", "WD2", "WD4", "WD6", "WD8")
plotPCA(vsd) + scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'), name = "Treatment")
```

<img src="rnaseq_files/figure-gfm/rna1-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/rnaseq_pca.pdf', height = 4.5, width = 5)
```

Check normalized data.

``` r
raw <- reshape2::melt(assay(dds))
colnames(raw) <- c("Gene","Sample","Expression")
ggplot(raw, aes(x=Sample, y=Expression)) + geom_violin() + theme(axis.text.x = element_text(angle=45, hjust=1))
```

![](rnaseq_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
matrix <- reshape2::melt(assay(vsd))
colnames(matrix) <- c("Gene","Sample","Expression")
ggplot(matrix, aes(x=Sample, y=Expression)) + geom_violin() + theme(axis.text.x = element_text(angle=45, hjust=1))
```

![](rnaseq_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

Generate all volcano plots. We use `ashr` to shrink fold-changes.
Vertical lines indicate \|log2FC\| \> 0.5. Horizontal lines indicate
padj \< 0.05.

``` r
mapping = key %>% dplyr::select(Ensembl.Gene.ID, Associated.Gene.Name) %>% distinct()

for(i in c(0, 2, 4, 6)) {
  for (j in c(2, 4, 6, 8)) {
    if(j > i) {
      res = lfcShrink(dds, contrast = c("condition", j, i), type = 'ashr')
      annotation = data.frame("Ensembl.Gene.ID" = rownames(res), stringsAsFactors = FALSE) %>% 
      left_join(., mapping, by = "Ensembl.Gene.ID") 
      
      # create custom key-value pairs for 'up-regulated', 'down-regulated', 'non-significant'
      keyvals <- ifelse(
        abs(res$log2FoldChange) > 0.5 & res$padj < 0.05, 'red', ifelse(
          res$padj < 0.05, 'blue', 'gray40'
        )
      )
      keyvals[is.na(keyvals)] <- 'gray40'
      names(keyvals)[keyvals == 'blue'] <- 'padj < 0.05'
      names(keyvals)[keyvals == 'red'] <- 'padj < 0.05 and log2FC > 0.5'
      names(keyvals)[keyvals == 'gray40'] <- 'NS'
      
      EnhancedVolcano(res, 
                lab = annotation$Associated.Gene.Name, 
                labSize = 6,
                x = 'log2FoldChange', 
                y = 'padj', 
                title = paste(i, 'vs.', j), 
                pCutoff = 0.05, 
                FCcutoff = 0.5, 
                colCustom = keyvals,
                colAlpha = 0.9, 
                pointSize = 3.0)
      ggsave(paste('../figures/volcano/', i, 'vs', j, '.pdf', sep=''), height = 7.5, width = 6.2)
    }
  }
}
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

Obtain log2FC, padj, and direction of FC for each gene between each
time-point comparison. Calculate the number of DE genes between each
time-point comparison. DE means padj \< 0.05 and \|log2FC\| \> 0.5.

``` r
mapping = key %>% dplyr::select(Ensembl.Gene.ID, Associated.Gene.Name) %>% distinct()
class_df = data.frame(matrix(NA, ncol=1, nrow=nrow(annotation)))[-1]
summary_df = data.frame("t1"=integer(), "t2"=integer(), "Total" =integer(), 
                        "Up"=integer(), "Down"=integer(), stringsAsFactors=FALSE)

for(i in c(0, 2, 4, 6)) {
  for (j in c(2, 4, 6, 8)) {
    if(j > i) {
      res = lfcShrink(dds, contrast = c("condition", j, i), type = 'ashr')
      fcname = paste('log2FC_', i, j, sep='')
      dename = paste('DE_', i, j, sep='')
      dirname = paste('Direct_', i, j, sep='')
  
      genes = as.data.frame(res) %>% 
        tibble::rownames_to_column(var = "Ensembl.Gene.ID") %>% 
        left_join(., mapping, by = "Ensembl.Gene.ID") %>%
        mutate(
          !!fcname := log2FoldChange,
          !!dename := ifelse(padj < 0.05 & abs(log2FoldChange) > 0.5, "yes", "no"),
          !!dirname := ifelse(log2FoldChange > 0, "up", "down")
        ) %>%
        dplyr::select(Associated.Gene.Name, !!as.symbol(fcname), 
                      !!as.symbol(dename), !!as.symbol(dirname))
  
      if(ncol(class_df) == 0) {
        class_df = cbind(class_df, genes)
      } else {
        stopifnot(identical(class_df$Associated.Gene.Name, genes$Associated.Gene.Name))
        class_df = cbind(class_df, genes %>% dplyr::select(-Associated.Gene.Name))
      }
      
      num_total = nrow(filter(genes, !!as.symbol(dename) == "yes"))
      num_up = nrow(filter(genes, !!as.symbol(dename) == "yes", !!as.symbol(dirname) == "up"))
      num_down = nrow(filter(genes, !!as.symbol(dename) == "yes", !!as.symbol(dirname) == "down"))
      num_category = data.frame("t1"=i, "t2"=j, "Total"=num_total, "Up"=num_up, 
                              "Down"=num_down, stringsAsFactors=FALSE)
      summary_df = rbind(summary_df, num_category)
    }
  }
}
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

Plot heat maps of significantly DE genes (padj \< 0.05, \|log2FC\| \>
0.5).

``` r
a <- summary_df %>% 
  ggplot(aes(x=t1, y=t2)) + 
  geom_tile(aes(fill=Total), color='black', size=0.5, alpha = 0.9) + 
  scale_fill_gradient(low="#f1f1f1", high="#696969") + 
  scale_x_continuous(breaks=seq(0,8,2)) + 
  scale_y_continuous(breaks=seq(0,8,2)) + 
  labs(x="Weeks on Western Diet", y="Weeks on Western Diet", title="Total") + 
  geom_text(aes(label=Total), size=6) + 
  theme(legend.position="none")

b <- summary_df %>% 
  ggplot(aes(x=t1, y=t2)) + 
  geom_tile(aes(fill=Up), color='black', size=0.5, alpha=0.9) +
  scale_fill_gradient(low="#e5e5ff", high="blue") + 
  scale_x_continuous(breaks=seq(0,8,2)) + 
  scale_y_continuous(breaks=seq(0,8,2)) + 
  labs(x="Weeks on Western Diet", y="Weeks on Western Diet", title="Up-regulated") +
  geom_text(aes(label=Up), size=6) + 
  theme(legend.position="none")

c <- summary_df %>% 
  ggplot(aes(x=t1, y=t2)) + 
  geom_tile(aes(fill=Down), color='black', size=0.5, alpha=0.9) +
  scale_fill_gradient(low="#ffe5e5", high="red") +
  scale_x_continuous(breaks=seq(0,8,2)) + 
  scale_y_continuous(breaks=seq(0,8,2)) + 
  labs(x="Weeks on Western Diet", y="Weeks on Western Diet", title="Down-regulated") +
  geom_text(aes(label=Down), size=6) + 
  theme(legend.position="none")

plot_grid(a,b,c, ncol = 3)
```

<img src="rnaseq_files/figure-gfm/rna2-1.png" style="display: block; margin: auto;" />

``` r
# ggsave('../figures/rnaseq_pairwise_de_genes_heatmap_fccutoff.pdf', height = 3.7, width = 10)
rm(a,b,c)
```

The RNA-seq profiles have a largely `biphasic` response. We can
categorize genes accordingly and count the number in each category.

``` r
class_df = class_df %>% 
  mutate(biphasic_class = ifelse(DE_02 == 'yes' & DE_28 == 'no', 'uniq_02',
                 ifelse(DE_02 == 'no' & DE_28 == 'yes', 'uniq_28',
                 ifelse(DE_02 == 'yes' & DE_28 == 'yes', 'both', 'none'))))

print(nrow(filter(class_df, biphasic_class == 'uniq_02')))
```

    ## [1] 659

``` r
print(nrow(filter(class_df, biphasic_class == 'uniq_28')))
```

    ## [1] 494

``` r
print(nrow(filter(class_df, biphasic_class == 'both'))) 
```

    ## [1] 85

``` r
print(nrow(filter(class_df, biphasic_class == 'none'))) 
```

    ## [1] 21494

Save lists of DE genes that fall into each category for downstream KEGG
pathway enrichment analysis.

``` r
save_list <- function(l, fname) {
  write.table(
    l,
    paste0('../processed_data/degenes/', fname, '.txt'),
    quote = F, col.names = F, row.names = F
  )
}

# DE exclusively between 02, up
up_02_list = (
  class_df %>% 
    filter(biphasic_class == 'uniq_02', Direct_02 == 'up') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(up_02_list, 'up_02')

# DE exclusively between 02, down
down_02_list = (
  class_df %>% 
    filter(biphasic_class == 'uniq_02', Direct_02 == 'down') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(down_02_list, 'down_02')

# DE exclusively between 28, up
up_28_list = (
  class_df %>% 
    filter(biphasic_class == 'uniq_28', Direct_28 == 'up') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(up_28_list, 'up_28')

# DE exclusively between 28, down
down_28_list = (
  class_df %>% 
    filter(biphasic_class == 'uniq_28', Direct_28 == 'down') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(down_28_list, 'down_28')

# DE between 02 (up) and 28 (up) 
up02_up28_list = (
  class_df %>% 
    filter(biphasic_class == 'both', Direct_02 == 'up', Direct_28 == 'up') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(up02_up28_list, 'up02_up28')

# DE between 02 (up) and 28 (down) 
up02_down28_list = (
  class_df %>% 
    filter(biphasic_class == 'both', Direct_02 == 'up', Direct_28 == 'down') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(up02_down28_list, 'up02_down28')

# DE between 02 (down) and 28 (up) 
down02_up28_list = (
  class_df %>% 
    filter(biphasic_class == 'both', Direct_02 == 'down', Direct_28 == 'up') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(down02_up28_list, 'down02_up28')

# DE between 02 (down) and 28 (down) 
down02_down28_list = (
  class_df %>% 
    filter(biphasic_class == 'both', Direct_02 == 'down', Direct_28 == 'down') %>% 
    dplyr::select(Associated.Gene.Name)
)
save_list(down02_down28_list, 'down02_down28')
```

We input the above gene lists into KEGG pathway enrichment via
`https://biit.cs.ut.ee/gprofiler/gost`. We set organism to Mouse.

We then filter out the outputted KEGG pathways by: 1. any obvious
irrelevant pathways (i.e., cancer, viral infection, etc.)

Ran on 5/22/23. These are in files `kegg_02_up.csv`, `kegg_02_down.csv`,
`kegg_28_up.csv`, and `kegg_28_down.csv`.

This chunk calculates and writes FPKM measurements for all genes.

``` r
# no need to run this if already saved dataframe.
ordered_geneLength = data.frame("Ensembl.Gene.ID" = rownames(assay(dds)), stringsAsFactors = FALSE) %>% 
    left_join(., gene_info, by = "Ensembl.Gene.ID")
mcols(dds)$basepairs = ordered_geneLength$gene_length 
fpkm_dds <- fpkm(dds)

ID_mapping = read.table('../raw_data/Mappings/rnaseq_mapping_list.txt', header = TRUE)

all_genes <- as.data.frame(fpkm_dds) %>% 
  tibble::rownames_to_column(var = "Ensembl.Gene.ID") %>%
  left_join(., gene_info, by = "Ensembl.Gene.ID") %>%
  dplyr::select(-c("Ensembl.Gene.ID", "gene_length")) %>%
  distinct(Associated.Gene.Name, .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "Associated.Gene.Name") %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "RNA_seq_ID") %>%
  left_join(., ID_mapping, by = "RNA_seq_ID") %>%
  relocate(Lipidomics_ID) %>%
  separate(Lipidomics_ID, c("Treatment", "ID"), sep = "_") %>%
  mutate(Treatment = ifelse(Treatment == 'CHOW8', 'Chow8', Treatment)) %>%
  dplyr::select(-RNA_seq_ID) %>%
  relocate(ID)

# write.csv(all_genes, '../processed_data/datasets/RNAseq_FPKM.csv', row.names = F)
```

``` r
# rmarkdown::render("rnaseq.Rmd")
# mv rnaseq.md ../markdowns/
# mv rnaseq_files ../markdowns/
```
