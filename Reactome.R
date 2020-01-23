# Reactome pathway Analysis (local dataset)- using enrichPathway & gsePathway

library(clusterProfiler)
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)
library(ReactomePA)

data_file <- readRDS("data14Features.rds")
analysis_data <-metadata(data_file)$limma[,"Sirt_vs_Mock",]%>%
  as.data.frame()%>%
  dplyr::select(effect,p)%>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)

analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>%
  stringr::str_split(.,";",simplify=TRUE)%>%
  .[,1]

geneList <- analysis_data %>% dplyr::filter(ENTREZ_ID != "" & effect != "")%>% 
  dplyr::select(ENTREZ_ID,effect)

gene_num_vec <- geneList [,2]
names(gene_num_vec) <- as.character(geneList [,1])
gene_num_vec <- sort(gene_num_vec,decreasing = TRUE, na.last = TRUE)

# Remove duplicate values
any(duplicated(names(gene_num_vec)))
gene_num_vec <- gene_num_vec[!duplicated(names(gene_num_vec))]

# Gene ID
gene <- names(gene_num_vec)[abs(gene_num_vec)> 2]

# --- Reactome Analysis -------
# - enrichPathway () - Over Representation Test 


reactome_enrich <- enrichPathway(gene = gene,
                                 organism = "mouse",
                                 pvalueCutoff = 0.05)

Reactome_enrich_result <- reactome_enrich@result

# - gsePathway () - Gene Set Enrichment Analysis -----

reactome_gsea <- gsePathway(geneList = gene_num_vec, 
                           organism = "mouse",
                           nPerm = 1000,
                           pvalueCutoff = 0.05,
                           verbose = TRUE)

reactome_gsea_result <- reactome_gsea@result

