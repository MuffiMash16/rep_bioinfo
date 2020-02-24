library(clusterProfiler)
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)

data_file <- readRDS("./data/data14Features.rds")
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

# --- KEGG Analysis -------
# Over Representation Test 
kk_enrich <- enrichKEGG(gene=gene,
                        organism = 'mmu',
                        pvalueCutoff = 0.05
)
kk_enrich_result <- kk_enrich@result

# --- KEGG gene set enrichment analysis -----
kk_gsea <- gseKEGG(geneList = gene_num_vec,
                   organism = 'mmu',
                   nPerm = 1000,
                   minGSSize = 120,
                   pvalueCutoff = 0.05,
                   verbose = FALSE)

kk_gsea_result <- kk_gsea@result


#---- MsigDb Analysis -----------
#Extract from MSigDb Analysis 
msigdbr_data <- msigdbr::msigdbr(species = "Mus musculus")%>%
  dplyr::select(gs_name, entrez_gene, gs_cat)

cat_vector <- c("H","C2","C5")

musculus_cat_data <- cat_vector %>% lapply(.,function(x){
  cat_based_dataset <- msigdbr_data %>% dplyr::filter(gs_cat == x)
  enrichr_cat_msigDB <- enricher(gene, TERM2GENE = cat_based_dataset)
  enrichr_cat_msigDB@result %>%
    mutate(Category = x) 
})

musculus_cat_df <- musculus_cat_data %>% do.call(rbind, .)

# Verify Results
H_data <- musculus_cat_df %>% filter(Category == 'H')
C2_data <- musculus_cat_df %>% filter(Category == 'C2')
C5_data <- musculus_cat_df %>% filter(Category == 'C5')




