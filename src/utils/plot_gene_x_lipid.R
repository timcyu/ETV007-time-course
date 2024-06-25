# heatmaps of gene x lipid correlations
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyverse)

# set working directory here
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# obtain top 50 most abundant TG species
TAG = read.csv(
  '../../processed_data/datasets/Lipidomics_liver_normliverweight.csv', 
  header = TRUE, 
  stringsAsFactors = FALSE, 
  check.names = FALSE
) %>% gather(
  lipid,
  value,
  3:ncol(.)
) %>% filter(
  grepl("TG", lipid)
) %>% group_by(
  lipid
) %>% mutate(
  med = median(value)
) %>% ungroup() %>% dplyr::select(
  lipid, med
) %>% distinct() %>% top_n(50, wt=med)

# extract remaining species that are not TAG
other_species = read.csv(
  '../../processed_data/datasets/Lipidomics_liver_normliverweight.csv', 
  header = TRUE, 
  stringsAsFactors = FALSE, 
  check.names = FALSE
) %>% gather(
  lipid,
  value,
  3:ncol(.)
) %>% filter(
  !grepl("TG", lipid)
) %>% dplyr::select(lipid) %>% distinct()

# list of lipids to (potentially) include in the heatmap
lipid_list = append(TAG$lipid, other_species$lipid)

# remove any correlations that don't contain lipids in `lipid_list`
top50gene_x_liver_lipid = read.csv(
  '../../processed_data/correlations/top50gene_x_liver_lipid.csv'
) %>% filter(
  Liver_lipid %in% lipid_list
) %>% drop_na()

# create gene x liver lipid matrix to plot
top50gene_x_liver_lipid_mat = top50gene_x_liver_lipid %>% dplyr::select(-adj_pval) %>%
  pivot_wider(names_from=Liver_lipid, values_from=bicor) %>%
  tibble::column_to_rownames(var="Gene") %>% as.matrix()

# plot heat map
colors <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))

lipidclasses = (data.frame(
  "lipid" = colnames(top50gene_x_liver_lipid_mat)
) %>% separate(
  col="lipid", 
  into=c("class"), 
  sep="\\.", 
  remove=FALSE, 
  extra="drop"
) %>% mutate(
  lipid_class = ifelse(
    grepl("PE[.]O", lipid), "PE_O", ifelse(
      grepl("PE[.]P", lipid), "PE_P", class
    )
  )
))$lipid_class

column_ha = HeatmapAnnotation(
  lipid_class = lipidclasses,
  col = list(lipid_class=c(
    "CE"="#A00000", "Cer"="#00C000", "DG"="#5757F9", 
    "FA"="#FF6000", "HexCER"="#0000C0", "LacCER"="#C0C000", 
    "LPC"="#D7B0B0", "LPE"="#9F044D", "PC"="#077E97", 
    "PE"="#C5944E", "PE_O"="#034E61", "PE_P"="#FFA040", 
    "SM"="#606060", "TG"="#fcfc81"
  )
  ),
  border = TRUE
)
set.seed(123)
ht = Heatmap(
  top50gene_x_liver_lipid_mat, 
  top_annotation=column_ha,
  cluster_rows=TRUE,
  cluster_columns=TRUE,
  show_row_dend=FALSE,
  show_column_dend=FALSE,
  show_column_names=FALSE,
  col=colors,
  column_split=lipidclasses,
  border=TRUE
)

# overall heatmap
pdf('../../figures/gene_heatmaps/gene_x_all_normliverweight.pdf', height=9, width=18)
ht = draw(ht)
dev.off()

# only plot DG
DG_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'DG')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='DG'
)
pdf('../../figures/gene_heatmaps/gene_x_DG.pdf', height=9, width=13)
DG_ht
dev.off()

# only plot Cer
Cer_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'Cer')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='Cer'
)
pdf('../../figures/gene_heatmaps/gene_x_Cer.pdf', height=9, width=7)
Cer_ht
dev.off()

# only plot FA
FA_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'FA')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='FA'
)
pdf('../../figures/gene_heatmaps/gene_x_FA.pdf', height=9, width=7)
FA_ht
dev.off()

# only plot LPC
LPC_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'LPC')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='LPC'
)
pdf('../../figures/gene_heatmaps/gene_x_LPC.pdf', height=9, width=7)
LPC_ht
dev.off()

# only plot PE
PE_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'PE')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='PE'
)
pdf('../../figures/gene_heatmaps/gene_x_PE.pdf', height=9, width=14)
PE_ht
dev.off()

# only plot PC
PC_ht = Heatmap(
  top50gene_x_liver_lipid_mat[,which(lipidclasses == 'PC')], 
  row_order=rownames(top50gene_x_liver_lipid_mat)[row_order(ht)],
  cluster_columns=TRUE,
  cluster_rows=FALSE,
  show_row_dend=FALSE,
  show_column_dend=FALSE, 
  show_column_names=TRUE,
  border=TRUE,
  col=colors,
  column_title ='PC'
)
pdf('../../figures/gene_heatmaps/gene_x_PC.pdf', height=9, width=16)
PC_ht
dev.off()