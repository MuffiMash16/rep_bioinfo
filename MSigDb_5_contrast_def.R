# MSigDb library 
# Script to run enricher( ) on Categories H, C2, C5
library(msigdbr)
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)

# 1. Load and prepare data set
data_file <- readRDS("~/data14Features.rds")
colnames(rowData(data_file))

analysis_data <- metadata(data_file)$limma[,"Sirt_vs_Mock",] %>%
  as.data.frame()%>%
  dplyr::select(effect,p)%>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)

analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>% 
  stringr::str_split(.,";",simplify = TRUE)%>%
  .[,1]

analysis_data <- analysis_data %>% mutate(p_GRP_ID = ifelse(p <= 0.05,1,
                         ifelse(p <= 0.05,2,3)))

# Step 3 - create Gene List 
geneList <- analysis_data %>% dplyr::filter(ENTREZ_ID != "" & effect != "")%>% 
  dplyr::select(ENTREZ_ID,effect)

# Step 4 - Numeric vector - effect value
gene_num_vec <- geneList [,2]
#        - Character vector - ENTREZ IDs
names(gene_num_vec) <- as.character(geneList [,1])
#        - Ordered Vector
gene_num_vec <- sort(gene_num_vec,decreasing = TRUE, na.last = TRUE)

# Step 5 - Remove duplicate values
any(duplicated(names(gene_num_vec)))
gene_num_vec <- gene_num_vec[!duplicated(names(gene_num_vec))]
geneID <- names(gene_num_vec)

# ----------------------------------------------------------
# 2.Extract from MSigDb Analysis 
msigdbr_data <- msigdbr(species = "Mus musculus")%>%
  dplyr::select(gs_name, entrez_gene, gs_cat)

any(msigdbr_data$gs_cat== "H")

cat_vector <- c("H","C2","C5")

musculus_cat_data <- cat_vector %>% lapply(.,function(x){
  cat_based_dataset <- msigdbr_data %>% dplyr::filter(gs_cat == x)
  enrichr_cat_msigDB <- enricher(geneID, TERM2GENE = cat_based_dataset)
  enrichr_cat_msigDB@result %>%
    mutate(gs_cat = x) %>%
    dplyr::mutate(p_log = -1*log(pvalue, 10))%>%
    dplyr::filter(p_log > 1.3) 
  })

musculus_cat_df <- musculus_cat_data %>% do.call(rbind, .)

# Verify results
H1_data <- musculus_cat_df %>% filter(gs_cat == 'H')

