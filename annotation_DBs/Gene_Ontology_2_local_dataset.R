# GO Analysis on Local Dataset
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)
library(org.Mm.eg.db)

# Prepare Data
data_file <- readRDS("data/data14Features.rds")

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

# create Gene ID - filter for effect greater than 2
gene <- names(gene_num_vec)[abs(gene_num_vec)> 2]
head(gene)

# gene df is used in enrichr ()
gene.df <- bitr(gene, 
                fromType = "ENTREZID",
                toType = c("ENSEMBL","SYMBOL"),
                OrgDb = org.Mm.eg.db)

head(gene.df)


# --- GO Classification ---
ggo <- groupGO(gene = gene,
        OrgDb = org.Mm.eg.db,
        ont = "CC",
        level = 3,
        readable = TRUE)
head(ggo)

# --- GO Over Representation ---
e_go <- enrichGO(gene = gene,
                 universe = names(gene_num_vec),
                 OrgDb = org.Mm.eg.db,
                 ont = "CC" ,
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE)

head(e_go)

e_go_2 <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENSEMBL",
                   ont = "CC" ,
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

head(e_go_2)

# Map Gene ID to Gene Symbol using setreadable() or readble = TRUE
e_go_2_symbol <- setReadable(e_go_2, OrgDb = org.Mm.eg.db)

# ----- dropGO function, user can remove specific GO terms / GO level 
# --- from results obtained from both enrichGO and compareCluster.
drop_GO <- dropGO(e_go, level = 4, term = NULL)

# ---- testGO at specific levels - gofilter to restrict the result at 
# --- specific GO level ------
filter_GO <- gofilter(e_go, level = 5)

# ---- simplify method to reduce redundant GO terms from the outputs of ----
simplify_GO <- simplify(e_go,
                        cutoff =0.7,
                        by="p.adjust",
                        select_fun = min,
                        measure = "Wang",
                        semData = NULL
)

# --- GO Gene Set Enrichment ---
gsea_go <- gseGO(geneList = gene_num_vec,
                 OrgDb = org.Mm.eg.db,
                 ont = "CC",
                 nPerm = 1000,
                 minGSSize = 100 ,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE )


# --- GO Semantic Similarity ---
# use it to cluster genes/proteins into different clusters based on their 
# functional similarity and to measure the similarities among GO terms to 
# reduce the redundancy of GO enrichment results.
go_similarity <- GOSemSim::godata(OrgDb = org.Mm.eg.db,
                                  keytype = "ENSEMBL",
                                  ont = "CC",
                                  computeIC = TRUE)