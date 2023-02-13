library(tidyverse)
library(ggfortify)
library(bnstruct)
library(viridis)
library(DescTools)
library(RColorBrewer)
require(cowplot)
theme_set(theme_cowplot())

"utility functions"

makeBubbleData <- function(df, islipidomics=TRUE) {
  "
  builds a summary dataframe for making bubble plots.
  
  args ---
  df : dataframe
    an omics dataset. Must contain the columns `ID` and
    `treatment`, while all other columns are a measured
    lipid, bile acid, etc.
  islipidomics : boolean
    denotes if the dataset is a lipidomics dataset.
  "
  if(islipidomics == TRUE) {
    bubble_df = df %>%
        gather(f, v, -c("ID", "Treatment")) %>%
        separate(col="f", into="class", sep="\\.", 
                 remove=FALSE, extra='drop') %>%
        mutate(Feat = ifelse(
          grepl("PE[.]O", f), "PE_O",
          ifelse(grepl("PE[.]P", f), "PE_P", 
                 class)
          )
        ) %>%
        select(-class) %>%
        group_by(ID, Feat) %>%
        drop_na() %>%
        mutate(Val = sum(v)) %>%
        ungroup() %>%
        select(ID, Treatment, Feat, Val) %>%
        distinct()
  }
  else {
    bubble_df = df %>% 
      gather(Feat, Val, -c("ID", "Treatment"))
  }
  return(bubble_df)
}


getBubbleSignificance <- function(bubble_df, thresh=0.05) {
  '
  performs Dunnetts test to obtain significance of difference
  in each feature at each time point (WD2,4,6,8) relative to 
  the control treatment (CHOW8).
  
  args ---
  bubble_df : dataframe
    dataframe made by makeBubbleData()
  thresh : float
    threshold for significance
  '
  bubble_df$Treatment <- as.factor(bubble_df$Treatment)
  featList = levels(factor(bubble_df$Feat))

  pvals = data.frame(
    "Comparison"=character(), "pval"=double(), 
    "Median"=double(), "FoldChange"=double(), 
    "Feat"=character(), stringsAsFactors=FALSE)
  
  for(i in featList) {
    res = DunnettTest(
      Val ~ Treatment, 
      data=filter(bubble_df, Feat == i),
      control = 'Chow8'
    )
    
    # fold-changes by median across mice of same treatment
    fc = bubble_df %>% 
      filter(Feat == i) %>% 
      group_by(Treatment) %>% 
      mutate(Median = median(Val)) %>% 
      ungroup() %>%
      select(Treatment, Feat, Median) %>% 
      distinct() %>% 
      mutate(
        FoldChange = Median / filter(., Treatment == 'Chow8')$Median
      ) %>%
      filter(Treatment != 'Chow8') %>% 
      mutate(
        Comparison = paste(Treatment, "-", "Chow8", sep = "")
      ) %>% 
      select(Comparison, Median, FoldChange, Feat)
    
    # combine p-values and fold-change into one data frame
    row = as.data.frame(res$Chow8) %>% 
      select(pval) %>% 
      tibble::rownames_to_column(var="Comparison") %>% 
      left_join(., fc, by="Comparison")
    
    pvals = rbind(pvals, row)
  }
  pvals = pvals %>% mutate(
    Significant = ifelse(pval < thresh, "yes", "no")
  )
  return(pvals)
}


makeBubblePlot <- function(
  pval_df, 
  colors=list(),
  y_lines=list(),
  bubble_range=list(),
  bubble_breaks=list(),
  ylim=list(),
  xlim=list()
  ) {
  '
  makes a bubble plot.
  
  args ---
    pval_df : dataframe made by getBubbleSignificance()
    colors : list of hex values
    y_lines : list. ex. c(0,2,4,6)
    bubble_range : list with two floats. ex. c(5,20)
    bubble_breaks : list. ex. c(0.01,0.1,1,5,9)
    ylim : list with two floats. ex. c(0,7)
    xlim : list with two floats. ex. c(-0.5,5.5)
  '
  if(length(colors) == 0) {
    colors = brewer.pal(12, "Set3")
    if(length(unique(pval_df$Feat)) > length(colors)) {
      return("Too many features. Specify custom color list")
    }
    else {
      colors = colors[0:length(unique(pval_df$Feat))]
    }
  }
  
  if(length(y_lines) == 0) {
    min_y = floor(min(log2(pval_df$FoldChange)))
    max_y = ceiling(max(log2(pval_df$FoldChange))) + 1
    y_lines = seq(min_y, max_y, by=2)
  }
  
  if(length(bubble_range) == 0) {
    min_size = ceiling(-log10(max(pval_df$pval)))
    max_size = ceiling(-log10(min(pval_df$pval)))
    bubble_range = c(min_size, max_size)
  }
  
  if(length(bubble_breaks) == 0) {
    if(bubble_range[1] > 10) {
      bubble_breaks = c(0.01, 0.1, 1, 10, bubble_range[1])
    } 
    else {
      bubble_breaks = c(0.01, 0.1, 1, 10)
    }
  }
  
  if(length(xlim) == 0) {
    xlim = c(
      floor(min(log2(pval_df$FoldChange))) - 1,
      ceiling(max(log2(pval_df$FoldChange))) + 1
    )
  }
  
  if(length(ylim) == 0) {
    ylim = c(
      floor(min(log10(pval_df$Median)))-1,
      ceiling(max(log10(pval_df$Median)))+1
    )
  }
  
  plot = pval_df %>%
    mutate(log2FoldChange = log2(FoldChange)) %>%
    mutate(log10pval = -log10(pval)) %>%
    mutate(log10Abundance = log10(Median)) %>%
    arrange(desc(log10Abundance)) %>%
    ggplot(
      aes(x=log2FoldChange, 
          y=log10Abundance, 
          size=log10pval,
          fill=Feat)
    ) +
    geom_vline(
      xintercept=0, 
      linetype='solid', 
      color='black'
    ) +
    geom_hline(
      yintercept=y_lines, 
      linetype='dotted',
      color='grey70',
      size=1
    ) +
    geom_point(
      aes(fill=Feat), 
      alpha=0.9,
      color='black',
      pch=21, 
      stroke=1
    ) +
    scale_size(
      range=bubble_range,
      breaks=bubble_breaks,
      name="Significance -log10(P)", 
    ) +
    ylim(ylim) +
    xlim(xlim) +
    scale_fill_manual(values=colors) +
    facet_grid(~Comparison, scales='fixed') +
    labs(x="log2(Fold-change)", y="log10(Abundance)")
  return(plot)
}

