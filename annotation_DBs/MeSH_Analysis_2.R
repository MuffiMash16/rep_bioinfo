# MeSH analysis using Local Dataset -----------------------------------------------------
# Data source from "gendoo", "gene2pubmed" and "RBBH" are all supported. User can select 
# interesting category to test. All 16 categories are supported. The analysis 
# supports >70 species listed in MeSHDb BiocView.

# Gene IDs were selected based on Index
# gsea () doesnt return results - tried out with different Categories
# dot plot - does not draw any diagrams

library(meshes)
library(MeSH.Mmu.eg.db)
library(clusterProfiler)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

# Dataset from DOSE
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
gene <- names(gene_num_vec)[1:100]
gene <- names(gene_num_vec)

# ----- enricher () ---------
en_x <- enrichMeSH(gene, MeSHDb = "MeSH.Mmu.eg.db", database='gendoo', category = 'C')
x_result <- en_x@result

# ----- gsea () ------------
gsea_y <- gseMeSH(gene_num_vec, 
                  MeSHDb = "MeSH.Mmu.eg.db", 
                  database = 'gene2pubmed', 
                  category = "V", 
                  by = "DOSE")

y_result <- gsea_y@result

#---- Visualization ----
#--- enrichplot (i.e.barplot, dotplot, cnetplot, emapplot and gseaplot) to visualize these enrichment results ---
goplot(en_x)
barplot(en_x, showCategory = 5)
dotplot(en_x)
gseaplot(y,y[1,1],title = y[1,2])