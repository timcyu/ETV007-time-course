Lipidomics analysis
================
Timothy Yu

``` r
knitr::opts_chunk$set(echo = TRUE)
```

Liver lipidomics

``` r
liverData = read.csv('../processed_data/datasets/Lipidomics_liver_normliverweight.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
```

``` r
autoplot(prcomp((liverData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = liverData, colour = 'Treatment', 
         size = 0.1, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm', 
         frame.level = 0.90) + 
         scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41')) +
         scale_fill_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# By 90% confidence ellipse, Mice 53 and 39 are outliers. We will remove them and re-plot the PCA. They are also excluded from the remaining analysis.
```

``` r
# PCA plot of liver lipidomics data
liverData = liverData %>% filter(!ID %in% c(39,53))
autoplot(prcomp((liverData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = liverData, colour = 'Treatment', size = 2.5) + 
  scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_pca.pdf', height = 3.5, width = 5)
```

``` r
# make bubble plot
colors <- c("#A00000", "#00C000", "#5757F9", "#FF6000", "#0000C0", "#C0C000", "#D7B0B0",
            "#9F044D", "#077E97", "#C5944E", "#034E61", "#FFA040", "#606060", "#fcfc81")

bubdata = makeBubbleData(liverData)
bubp = getBubbleSignificance(bubdata)
makeBubblePlot(bubp, abundance=FALSE, colors=colors, bubble_range=c(3,10))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_bubble_plot.pdf', height = 5.5, width = 15)
#write.csv(bubp, '../processed_data/misc/lipidomics_liver_bubble_stats.csv')
```

``` r
# triglyceride abundance heat map for a specific TG class (TG48)
scaled_liverData = cbind(liverData[,c(1,2)], scale(liverData[,-c(1,2)])) # mean 0, sd 1

TG48_data = scaled_liverData %>%
  select_if(~ !any(is.na(.))) %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class", "Chain"), sep = "\\.", 
          remove = FALSE, extra = 'drop') %>%
  filter(Class == 'TG', Chain == 48)

order = c(unique((TG48_data %>% filter(Treatment == 'Chow8'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD2'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD4'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD6'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD8'))$ID))

TG48_data %>%
  ggplot(aes(x = Feature, y = as.factor(ID))) + 
  geom_tile(aes(fill = Value), alpha = 1) + 
  scale_fill_distiller(palette='RdYlBu') +
  scale_y_discrete(limits = rev(factor(order))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank())
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_TG48_normliverweight_heatmap.pdf', height = 14, width = 9)
```

``` r
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
  geom_line(data = . %>% filter(Saturation == "PUFA"), size=0.5, alpha=0.7) +
  geom_line(data = . %>% filter(Saturation == "MUFA"), size=0.5, alpha=0.7) +
  geom_line(data = . %>% filter(Saturation == "SFA"), size=0.5, alpha=0.7) +
  scale_color_manual(values = c("SFA" = "#3969AC", "MUFA" = "#F2B701", "PUFA" = "#E73F74")) +
  labs(x = "Weeks on Western Diet",
       y = "Mean abundance (z-score)",
       color = "Saturation") +
  theme(legend.position = "right") +
  facet_wrap(~Chain, scales = 'free')
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_TG_saturation_by_chain.pdf', height = 7, width = 11)
```

``` r
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

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_TG_saturation_change.pdf', height = 3, width = 9)
```

``` r
TG_data = scaled_liverData %>%
  select_if(~ !any(is.na(.))) %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class", "Chain"), sep = "\\.", 
          remove = FALSE, extra = 'drop') %>%
  filter(Class == 'TG') %>%
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
  ))

order = c(unique((TG48_data %>% filter(Treatment == 'Chow8'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD2'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD4'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD6'))$ID),
          unique((TG48_data %>% filter(Treatment == 'WD8'))$ID))

TG_data %>%
  ggplot(aes(x = Feature, y = as.factor(ID))) + 
  geom_tile(aes(fill = Value), alpha = 1) + 
  scale_fill_distiller(palette='RdYlBu', limits = c(NA, 3), oob = scales::squish) +
  scale_y_discrete(limits = rev(factor(order))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  facet_wrap(~Saturation, scales = 'free')
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_allTG_heatmap.pdf', height = 15, width = 30)
```

``` r
TG48_data %>% 
  group_by(Treatment, Feature) %>%
  mutate(Mean_abundance = mean(Value)) %>%
  select(Treatment, Feature, Mean_abundance) %>%
  distinct() %>%
  mutate(Saturation = case_when(
    str_extract(Feature, "FA\\d+\\.\\d+$") %>% 
      str_extract("\\d+$") %>% 
      as.numeric() == 0 ~ "SFA",
    str_extract(Feature, "FA\\d+\\.\\d+$") %>% 
      str_extract("\\d+$") %>% 
      as.numeric() == 1 ~ "MUFA",
    str_extract(Feature, "FA\\d+\\.\\d+$") %>% 
      str_extract("\\d+$") %>% 
      as.numeric() > 1 ~ "PUFA",
    TRUE ~ NA_character_
  )) %>%
  mutate(Treatment = fct_recode(Treatment,
                             "0" = "Chow8",
                             "2" = "WD2",
                             "4" = "WD4",
                             "6" = "WD6",
                             "8" = "WD8")) %>%
  ggplot(aes(x = Treatment, y = Mean_abundance, color = Saturation, group = Feature)) +
  geom_line(size=1, alpha=0.7) +
  scale_color_manual(values = c("SFA" = "#3969AC", "MUFA" = "#F2B701", "PUFA" = "#E73F74")) +
  labs(x = "Weeks on Western Diet",
       y = "Mean abundance (z-score)",
       color = "Saturation") +
  theme(legend.position = "right")
```

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: WD2, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD4, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD6, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD8

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

    ## Warning: Unknown levels in `f`: Chow8, WD2, WD4, WD6

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_liver_TG48_saturation_change.pdf', height = 5, width = 7)
```

Plasma lipidomics

``` r
plasmaData = read.csv('../processed_data/datasets/Lipidomics_plasma.csv', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
```

``` r
autoplot(prcomp((plasmaData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), 
         data = plasmaData, colour = 'Treatment', 
         size = 0.1, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm', 
         frame.level = 0.95) + 
         scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41')) +
         scale_fill_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# PCA plot of plasma lipidomics data
autoplot(prcomp((plasmaData[-c(1,2)] %>% select_if(~ !any(is.na(.)))), scale. = TRUE), data = plasmaData, colour = 'Treatment', size = 2.5) + scale_color_manual(values = c('black', '#3b9bb3', '#75af4f', '#e39225', '#c13b41'))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_plasma_pca.pdf', height = 3.5, width = 5)
```

``` r
# make bubble plot
colors <- c("#A00000", "#00C000", "#5757F9", "#FF6000", "#0000C0", "#C0C000", "#D7B0B0",
            "#9F044D", "#077E97", "#606060", "#fcfc81")

bubdata = makeBubbleData(plasmaData) %>% filter(!Feat %in% c("PE", "PE_O", "PE_P"))
bubp = getBubbleSignificance(bubdata)
makeBubblePlot(bubp, abundance=FALSE, colors=colors, bubble_range=c(3,10))
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_plasma_bubble_plot.pdf', height = 5, width = 15)
#write.csv(bubp, '../processed_data/misc/lipidomics_plasma_bubble_stats.csv')
```

``` r
# CE abundance heat map for a specific CE (CE(18:1))
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

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_plasma_CE_heatmap.pdf', height = 14, width = 9)
```

``` r
# corresponding plot for plasma -- triglyceride abundance heat map for a specific TG class (TG48)
scaled_plasmaData = cbind(plasmaData[,c(1,2)], scale(plasmaData[,-c(1,2)])) # mean 0, sd 1

plasma_TG48_data = scaled_plasmaData %>%
  select_if(~ !any(is.na(.))) %>%
  gather(Feature, Value, 3:ncol(.)) %>%
  separate(col = 'Feature', into = c("Class", "Chain"), sep = "\\.", 
          remove = FALSE, extra = 'drop') %>%
  filter(Class == 'TG', Chain == 48) %>%
  filter(!ID %in% c(55,56)) # 2 outliers

order = c(unique((plasma_TG48_data %>% filter(Treatment == 'Chow8'))$ID),
          unique((plasma_TG48_data %>% filter(Treatment == 'WD2'))$ID),
          unique((plasma_TG48_data %>% filter(Treatment == 'WD4'))$ID),
          unique((plasma_TG48_data %>% filter(Treatment == 'WD6'))$ID),
          unique((plasma_TG48_data %>% filter(Treatment == 'WD8'))$ID))

plasma_TG48_data %>%
  ggplot(aes(x = Feature, y = as.factor(ID))) + 
  geom_tile(aes(fill = Value), alpha = 1) + 
  scale_fill_distiller(palette='RdYlBu') +
  scale_y_discrete(limits = rev(factor(order))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank())
```

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
#ggsave('../figures/lipidomics_plasma_TG48_heatmap.pdf', height = 14, width = 9)
```

``` r
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

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

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

![](/Users/timyu/Desktop/ETV007-time-course/markdowns/lipidomics_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggsave('../figures/lipidomics_plasma_liver_total_corr_heatmap.pdf', height = 4, width = 6)
```
