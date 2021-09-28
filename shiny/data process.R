#global for Bethan data shiny app
setwd("C:/Users/baolongsu/Desktop/Projects/BethanApp/test")

library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(data.table)
library(viridis)

#prepare data and save as .rdata
gene_lipid.corr <- read.csv("data/GENE-LIVERLIPID-corrTable_FDR_0.05.csv")
rnaseq <- read.csv("data/RNA-seq-FPKM-with-gene-names.csv")
rnaseq<- rnaseq  %>%
  separate(X, c("AB", "Index", "Group"), "_") 
rnaseq$Group <-  sub(".count", "", rnaseq$Group)

#test data.table
DT.gene_lipid.corr <- as.data.table(gene_lipid.corr)
DT.gene_lipid.corr %>% select(Gene)

#convert gene_lipid.corr to data.table
gene_lipid.corr <- as.data.table(gene_lipid.corr)
rm(DT.gene_lipid.corr)

#this merged data is processed by "C:\Users\baolongsu\Desktop\Projects\BethanApp\ETV007-time-course\src\correlations.Rmd"
rna_livlipid_merge <- read.csv("data/nraseq_liverlipid_merge_ori.csv")
#add sample groupname column
#rna_livlipid_merge["Group"] <- str_extract(rna_livlipid_merge$sampleID, "^[^_]*")
#rna_livlipid_merge$Group
#rna_livlipid_merge<-rna_livlipid_merge[c("Group", colnames(rna_livlipid_merge)[-c(1,23308)])]
write.csv(rna_livlipid_merge, "data/nraseq_liverlipid_merge.csv") #23306 columns, 46 rows

#get gene list
genelist <- colnames(rna_livlipid_merge)[c(3:22730)]


save(gene_lipid.corr, rna_livlipid_merge,genelist, #rnaseq,
     file="data/mydata.rdata")








#load from .rdata
load("data/mydata.rdata")

gene.box <- rnaseq %>% select("Group", "Apoh", "AB", "Index") %>%
  ggplot(aes(x = Group, y = !!as.symbol("Apoh"), color=Group,
                    text = paste(AB, Index))) +
  geom_boxplot(outlier.shape = NA, color="black", outlier.color = NA,outlier.size = 0) +
  geom_jitter(size = 2, alpha = 0.7,
              position=position_jitter(width=0.1, height=0.1))+
  scale_colour_viridis(discrete = T)+
  theme_classic()

#p <- ggplot(data=data, aes_string(x='Month', y = input$Measure)) +geom_line(size = 1.5) + theme_minimal()
#gene.box <- ggplotly(gene.box)
gene.box <- ggplotly(gene.box) %>% layout(legend = list(orientation = "h", x = 0, y = -0.1))

#remove outlier from boxplot
gene.box$x$data[1] <- lapply(gene.box$x$data[1], 
                          FUN = function(x){
                            x$marker = list(opacity = 0)
                            return(x)
                            })
gene.box



#lipid box

rna_livlipid_merge[c(1,2,2+22730)] %>%
  ggplot(aes(x = Group, y = !!as.symbol(names(rna_livlipid_merge)[2+22730]), color=Group,
             text = paste(sampleID))) +
  geom_boxplot(outlier.shape = NA, color="black", outlier.color = NA,outlier.size = 0)+
  geom_jitter(size = 2, alpha = 0.7,
              position=position_jitter(width=0.1, height=0,seed=123),
              na.rm=TRUE)+
  scale_colour_viridis(discrete = T)+
  theme_classic()






p <- rna_livlipid_merge %>% select("Group", input$gene, "sampleID") %>%
  ggplot(aes(x = Group, y = !!as.symbol(input$gene), color=Group,
             text=paste("ID ",sampleID))) +
  geom_boxplot(color="black") +
  geom_jitter(size = 2, alpha = 0.7,
              position=position_jitter(width=0.1, height=0,seed=123),
              na.rm=TRUE)+
  scale_colour_viridis(discrete = T)+
  theme_classic()+
  theme(legend.position="bottom")
p <- ggplotly(p, height = 600, width = 900)