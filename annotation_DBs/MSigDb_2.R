# Comments :-d default "p value cutoff" in the GSEA function was modified
# Category 6 (contains cellualr pathways often disregulated in cancer) error when perfomring gsea()

# MSigDb Analysis - 14 Feature Data set - Mus musculus
library(msigdbr)

# Step 1 - from "MSigDb" select a specific collection use Category (6) 
m_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

# Step 2 - Local DataSet
data_file <- readRDS("data14Features.rds")

analysis_data <- metadata(data_file)$limma[,"Sirt_vs_Mock",] %>%
  as.data.frame()%>%
  dplyr::select(effect)%>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)
  
analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>% 
  stringr::str_split(.,";",simplify = TRUE)%>%
  .[,1]

# Step 3 - Gene and Gene List created 
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

# Step 6 - Run enricher and GSEA function - msigdbr data vs. local data
enrichr_msigDB <- enricher(geneID, TERM2GENE = m_t2g)

GSEA_msigDB <- GSEA(gene_num_vec,pvalueCutoff = 0,TERM2GENE = m_t2g,verbose = FALSE, by = 'DOSE')
View(GSEA_msigDB)

# ------ Verification of values in variables ----
head(enrichr_msigDB@result) 
head(GSEA_msigDB)

