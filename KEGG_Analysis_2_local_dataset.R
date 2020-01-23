
library(clusterProfiler)
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)

#Search by KEGG code or Scientific name
search_kegg_organism("mmu",by="kegg_code")

mmu <- search_kegg_organism('Mus musculus',by = 'scientific_name')
dim(mmu)
head(mmu)

# ----- Dataset -------------
data_file <- readRDS("../data/data14Features.rds")
analysis_data <- data_file@metadata$limma[, "Sirt_vs_Mock", ]%>%
  as.data.frame() %>%
  dplyr::select(effect,p) %>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)

analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>%
  stringr::str_split(.,";",simplify=TRUE) %>%
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

# --- KEGG over Representation Test ------
kk_enrich <- enrichKEGG(gene=gene,
                        organism = 'mmu',
                        pvalueCutoff = 0.05
)
head(kk_enrich)

# --- KEGG gene set enrichment analysis -----
kk_gsea <- gseKEGG(geneList = gene_num_vec,
                   organism = 'mmu',
                   nPerm = 1000,
                   minGSSize = 120,
                   pvalueCutoff = 0.05,
                   verbose = FALSE)
head(kk_gsea)

# --- KEGG Module - over representation ---
mkk_enrich <- enrichKEGG(gene= gene,
                         organism = 'mmu',
)
head(mkk_enrich)

# --- KEGG Module - gene set enrichment analaysis ----
mkk_gsea <- gseMKEGG(geneList = gene_num_vec, organism = 'mmu')

