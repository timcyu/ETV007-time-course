Analysis of lipidomics data
================
Timothy Yu

This notebook analyzes the liver and plasma lipidomics datasets.

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] rmarkdown_2.11     cowplot_1.1.1      RColorBrewer_1.1-3 DescTools_0.99.44 
    ##  [5] viridis_0.6.2      viridisLite_0.4.0  bnstruct_1.0.11    igraph_1.2.9      
    ##  [9] bitops_1.0-7       ggfortify_0.4.13   forcats_0.5.1      stringr_1.4.0     
    ## [13] dplyr_1.0.9        purrr_0.3.4        readr_2.1.1        tidyr_1.2.0       
    ## [17] tibble_3.1.8       ggplot2_3.3.6      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fs_1.5.2          lubridate_1.8.0   httr_1.4.2        tools_4.0.2      
    ##  [5] backports_1.4.0   utf8_1.2.2        R6_2.5.1          DBI_1.1.1        
    ##  [9] colorspace_2.0-3  withr_2.5.0       tidyselect_1.1.2  gridExtra_2.3    
    ## [13] Exact_3.1         compiler_4.0.2    textshaping_0.3.6 cli_3.3.0        
    ## [17] rvest_1.0.2       expm_0.999-6      xml2_1.3.3        labeling_0.4.2   
    ## [21] scales_1.2.0      mvtnorm_1.1-3     proxy_0.4-26      systemfonts_1.0.4
    ## [25] digest_0.6.29     pkgconfig_2.0.3   htmltools_0.5.2   dbplyr_2.1.1     
    ## [29] fastmap_1.1.0     highr_0.9         rlang_1.0.4       readxl_1.3.1     
    ## [33] rstudioapi_0.13   generics_0.1.3    farver_2.1.1      jsonlite_1.7.2   
    ## [37] magrittr_2.0.3    Matrix_1.3-4      Rcpp_1.0.7        munsell_0.5.0    
    ## [41] fansi_1.0.3       lifecycle_1.0.1   stringi_1.7.6     yaml_2.2.1       
    ## [45] MASS_7.3-54       rootSolve_1.8.2.3 plyr_1.8.6        grid_4.0.2       
    ## [49] crayon_1.4.2      lmom_2.8          lattice_0.20-45   haven_2.4.3      
    ## [53] hms_1.1.1         knitr_1.36        pillar_1.8.0      boot_1.3-28      
    ## [57] gld_2.6.3         reshape2_1.4.4    reprex_2.0.1      glue_1.6.2       
    ## [61] evaluate_0.14     data.table_1.14.2 modelr_0.1.8      vctrs_0.4.1      
    ## [65] tzdb_0.2.0        cellranger_1.1.0  gtable_0.3.0      assertthat_0.2.1 
    ## [69] xfun_0.28         broom_0.7.10      e1071_1.7-9       ragg_1.2.2       
    ## [73] class_7.3-19      ellipsis_0.3.2

``` r
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))
```

Read in liver lipidomics dataset.

``` r
liverData = read.csv('../processed_data/datasets/Lipidomics_liver_normliverweight.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
```

By 90% confidence ellipse, Mice 53 and 39 are outliers. We will remove
them and re-plot the PCA. They are excluded from the remaining analysis.

``` r
autoplot(prcomp((liverData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = liverData, colour = 'Treatment', 
         size = 0.1, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm', 
         frame.level = 0.90) + 
         scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41')) +
         scale_fill_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

<img src="lipidomics_files/figure-gfm/lipid1-1.png" style="display: block; margin: auto;" />

PCA plot of liver lipidomics data.

``` r
liverData = liverData %>% filter(!ID %in% c(39,53))
autoplot(prcomp((liverData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = liverData, colour = 'Treatment', size = 2.5) + 
  scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

<img src="lipidomics_files/figure-gfm/lipid2-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_liver_pca.pdf', height = 3.5, width = 5)
```

Bubble plot showing relative fold-change in liver lipid abundance.

``` r
colors <- c("#A00000", "#00C000", "#5757F9", "#FF6000", "#0000C0", "#C0C000", "#D7B0B0",
            "#9F044D", "#077E97", "#C5944E", "#034E61", "#FFA040", "#606060", "#fcfc81")

bubdata = makeBubbleData(liverData)
bubp = getBubbleSignificance(bubdata)
makeBubblePlot(bubp, abundance=FALSE, colors=colors, bubble_range=c(3,10))
```

<img src="lipidomics_files/figure-gfm/lipid3-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_liver_bubble_plot.pdf', height = 5.5, width = 15)
#write.csv(bubp, '../processed_data/misc/lipidomics_liver_bubble_stats.csv')
```

Abundance of individual hepatic TG species split by SFA, MUFA, and PUFA.

``` r
scaled_liverData = cbind(liverData[,c(1,2)], scale(liverData[,-c(1,2)])) # mean 0, sd 1

mean_TG_data = scaled_liverData %>%
  select_if(~ !any(is.na(.))) %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class", "Chain"), sep = "\\.", 
          remove = FALSE, extra = 'drop') %>%
  filter(Class == 'TG') %>%
  group_by(Treatment, Feature, Chain) %>%
  mutate(Mean_abundance = mean(Value)) %>%
  select(Treatment, Feature, Chain, Mean_abundance) %>%
  distinct() %>%
  mutate(Saturation = case_when(
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() == 0 ~ "SFA",
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() == 1 ~ "MUFA",
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() > 1 ~ "PUFA",
    TRUE ~ NA_character_
  )) %>% mutate(Treatment = as.factor(Treatment)) %>% mutate(
    Treatment = fct_recode(
      Treatment,
      "0" = "Chow8",
      "2" = "WD2",
      "4" = "WD4",
      "6" = "WD6",
      "8" = "WD8"
  )
)

mean_TG_data %>% 
  ggplot(aes(x = Treatment, y = Mean_abundance, color = Saturation, group = Feature)) +
  geom_line(size=0.5, alpha=0.7) +
  scale_color_manual(values = c("SFA" = "#3969AC", "MUFA" = "#F2B701", "PUFA" = "#E73F74")) +
  labs(x = "Weeks on Western Diet",
       y = "Mean abundance (z-score)",
       color = "Saturation") +
  theme(legend.position = "right") +
  facet_wrap(~Saturation, scales = 'free')
```

<img src="lipidomics_files/figure-gfm/lipid4-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_liver_TG_saturation_change.pdf', height = 3, width = 9)
```

Read in plasma lipidomics dataset.

``` r
plasmaData = read.csv('../processed_data/datasets/Lipidomics_plasma.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
```

Plot PCA. There are no outliers.

``` r
autoplot(prcomp((plasmaData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = plasmaData, colour = 'Treatment', 
         size = 0.1, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm', 
         frame.level = 0.95) + 
         scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41')) +
         scale_fill_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

<img src="lipidomics_files/figure-gfm/lipid5-1.png" style="display: block; margin: auto;" />
PCA plot of plasma lipidomics data

``` r
autoplot(prcomp((plasmaData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), data = plasmaData, colour = 'Treatment', size = 2.5) + scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

<img src="lipidomics_files/figure-gfm/lipid6-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_plasma_pca.pdf', height = 3.5, width = 5)
```

Bubble plot showing relative fold-change in plasma lipid abundance.

``` r
colors <- c("#A00000", "#00C000", "#5757F9", "#FF6000", "#0000C0", "#C0C000", "#D7B0B0",
            "#9F044D", "#077E97", "#606060", "#fcfc81")

bubdata = makeBubbleData(plasmaData) %>% filter(!Feat %in% c("PE", "PE_O", "PE_P"))
bubp = getBubbleSignificance(bubdata)
makeBubblePlot(bubp, abundance=FALSE, colors=colors, bubble_range=c(3,10))
```

<img src="lipidomics_files/figure-gfm/lipid7-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_plasma_bubble_plot.pdf', height = 5, width = 15)
#write.csv(bubp, '../processed_data/misc/lipidomics_plasma_bubble_stats.csv')
```

Plot abundance of plasma CE species.

``` r
abundant_CE = plasmaData %>% select(grep("CE", colnames(.))) %>%
  colMeans() %>% as.data.frame() %>%
  filter(. > 50)

scaled_plasmaData = cbind(plasmaData[,c(1,2)], 
                    scale(plasmaData %>% select(rownames(abundant_CE)))) # mean 0, sd 1

CE_data = scaled_plasmaData %>%
  gather(Feature, Value, 3:ncol(.))

order = c(unique((CE_data %>% filter(Treatment == 'Chow8'))$ID),
          unique((CE_data %>% filter(Treatment == 'WD2'))$ID),
          unique((CE_data %>% filter(Treatment == 'WD4'))$ID),
          unique((CE_data %>% filter(Treatment == 'WD6'))$ID),
          unique((CE_data %>% filter(Treatment == 'WD8'))$ID))

CE_data %>% ggplot(aes(x = Feature, y = as.factor(ID))) + 
  geom_tile(aes(fill = Value), alpha = 1) + 
  scale_fill_distiller(palette='RdYlBu') +
  scale_y_discrete(limits = rev(factor(order))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank())
```

<img src="lipidomics_files/figure-gfm/lipid8-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_plasma_CE_heatmap.pdf', height = 14, width = 9)
```

Abundance of individual plasma TG species split by SFA, MUFA, and PUFA.

``` r
scaled_plasmaData = cbind(plasmaData[,c(1,2)], scale(plasmaData[,-c(1,2)])) # mean 0, sd 1

mean_plasma_TG_data = scaled_plasmaData %>%
  select_if(~ !any(is.na(.))) %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class", "Chain"), sep = "\\.", 
          remove = FALSE, extra = 'drop') %>%
  filter(Class == 'TG') %>%
  group_by(Treatment, Feature, Chain) %>%
  mutate(Mean_abundance = mean(Value)) %>%
  select(Treatment, Feature, Chain, Mean_abundance) %>%
  distinct() %>%
  mutate(Saturation = case_when(
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() == 0 ~ "SFA",
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() == 1 ~ "MUFA",
    str_extract(Feature, "TG\\.\\d+\\.\\d+") %>% 
      str_extract("\\.\\d+$") %>%
      str_replace("\\.", "") %>%
      as.numeric() > 1 ~ "PUFA",
    TRUE ~ NA_character_
  )) %>% mutate(Treatment = as.factor(Treatment)) %>% mutate(
    Treatment = fct_recode(
      Treatment,
      "0" = "Chow8",
      "2" = "WD2",
      "4" = "WD4",
      "6" = "WD6",
      "8" = "WD8"
  )
)

mean_plasma_TG_data %>% 
  ggplot(aes(x = Treatment, y = Mean_abundance, color = Saturation, group = Feature)) +
  geom_line(size=0.5, alpha=0.7) +
  scale_color_manual(values = c("SFA" = "#3969AC", "MUFA" = "#F2B701", "PUFA" = "#E73F74")) +
  labs(x = "Weeks on Western Diet",
       y = "Mean abundance (z-score)",
       color = "Saturation",
       title = "Plasma TG") +
  theme(legend.position = "right") +
  facet_wrap(~Saturation, scales = 'free')
```

<img src="lipidomics_files/figure-gfm/lipid9-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_plasma_TG_saturation_change.pdf', height = 3, width = 9)
```

Spearman correlation between hepatic and plasma lipids.

``` r
# correlate total lipid classes between liver and plasma
totalLipidbyClass_liver = liverData %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class"), sep = "\\.", remove = FALSE, extra = 'drop') %>%
  mutate(RefinedClass = ifelse(grepl("PE[.]O", Feature), "PE_O",
                        ifelse(grepl("PE[.]P", Feature), "PE_P", Class))) %>%
  select(-Class) %>%
  group_by(ID, RefinedClass) %>%
  drop_na() %>%
  mutate(totalLipid = sum(Value)) %>%
  ungroup() %>%
  select(ID, Treatment, RefinedClass, totalLipid) %>%
  distinct() %>%
  group_by(RefinedClass) %>%
  mutate(n_groups = length(unique(Treatment))) %>%
  filter(n_groups == 5) %>%
  ungroup() %>%
  select(-n_groups)

totalLipidbyClass_plasma = plasmaData %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class"), sep = "\\.", remove = FALSE, extra = 'drop') %>%
  mutate(RefinedClass = ifelse(grepl("PE[.]O", Feature), "PE_O",
                        ifelse(grepl("PE[.]P", Feature), "PE_P", Class))) %>%
  select(-Class) %>%
  group_by(ID, RefinedClass) %>%
  drop_na() %>%
  mutate(totalLipid = sum(Value)) %>%
  ungroup() %>%
  select(ID, Treatment, RefinedClass, totalLipid) %>%
  distinct() %>%
  group_by(RefinedClass) %>%
  mutate(n_groups = length(unique(Treatment))) %>%
  filter(n_groups == 5) %>%
  ungroup() %>%
  select(-n_groups)
 
lipids_in_common = intersect(totalLipidbyClass_liver$RefinedClass, totalLipidbyClass_plasma$RefinedClass)
 
df = data.frame("LipidClass" = character(), "Spearman_corr" = numeric(), "p_val" = numeric())
for (lipid1 in lipids_in_common) {
  lipid_in_liver = totalLipidbyClass_liver %>% filter(RefinedClass == lipid1)
  for (lipid2 in lipids_in_common) {
    lipid_in_plasma = totalLipidbyClass_plasma %>% filter(RefinedClass == lipid2)
    combined = left_join(lipid_in_liver, lipid_in_plasma, by = 'ID')
    corr = cor.test(combined$totalLipid.x, 
                    combined$totalLipid.y, 
                    method=c("pearson", "kendall", "spearman"))$estimate
    p = cor.test(combined$totalLipid.x, 
                 combined$totalLipid.y, 
                 method=c("pearson", "kendall", "spearman"))$p.value
    df = rbind(df, data.frame("Lipid_liver" = lipid1,
                              "Lipid_plasma" = lipid2,
                              "Spearman_corr" = corr,
                              "p_val" = p))
  }
}

# save dataframe
write.csv(df, '../processed_data/misc/liver_plasma_lipid_class_pairwise_corr.csv', row.names=FALSE)

df %>%
  ggplot(aes(x=Lipid_liver, y=Lipid_plasma)) +
  geom_tile(aes(fill=Spearman_corr), alpha=1, color='black', size=0.1) +
  scale_fill_gradient2(low = "#0000FF", mid = "#FFFFFF", high ="#FF0000", midpoint=0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

<img src="lipidomics_files/figure-gfm/lipid10-1.png" style="display: block; margin: auto;" />

``` r
#ggsave('../figures/lipidomics_plasma_liver_total_corr_heatmap.pdf', height = 4, width = 6)
```

``` r
# rmarkdown::render("lipidomics.Rmd")
# mv lipidomics.Rmd ../markdowns/
# mv lipidomics_files ../markdowns/
```
