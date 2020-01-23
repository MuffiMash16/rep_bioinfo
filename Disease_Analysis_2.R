# Disease Analysis on local Dataset 
# Code covers - DO + DGN - enrichr and gsea
# "DO" enrich - test throws an error - "No gene can be mapped"
# "DO" gsea -error - "no term enriched under specific pvalueCutoff..."
# "DGN" enrich & gsea - both throws errors

library(DOSE)
library(dplyr)

# ------------ Local Dataset ----------------------------------------
data_file <- readRDS("data14Features.rds")
colnames(rowData(data_file))

analysis_data <- metadata(data_file)$limma[,"Sirt_vs_Mock",] %>%
  as.data.frame()%>%
  dplyr::select(effect)%>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)

analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>% 
  stringr::str_split(.,";",simplify = TRUE)%>%
  .[,1]

geneList <- analysis_data %>% dplyr::filter(ENTREZ_ID != "" & effect != "")%>% 
  dplyr::select(ENTREZ_ID,effect)

# ---- Numeric vector - effect value
gene_num_vec <- geneList [,2]
#        - Character vector - ENTREZ IDs
names(gene_num_vec) <- as.character(geneList [,1])
#        - Ordered Vector
gene_num_vec <- sort(gene_num_vec,decreasing = TRUE, na.last = TRUE)

# Remove duplicate values
any(duplicated(names(gene_num_vec)))
gene_num_vec <- gene_num_vec[!duplicated(names(gene_num_vec))]

# create Gene ID - filter for effect greater than 1.5
gene <- names(gene_num_vec)[abs(gene_num_vec)> 1.5]
head(gene)



#-------*****-------- DO (Disease Ontology) --------*****-----------------------
# enrichDo function - is useful for identifying disease association of interesting genes
x <- enrichDO(gene          = gene,
              ont           = "DOLite",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "bonferroni",
              universe      = names(gene_num_vec),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
x_df <- x@result


# ---Check for Gene IDs displayed in the console- error message
which(apply(geneList, 1, function(r) any(r %in% c("8567", "727","5959","27241","7097","4627"))))
which(apply(geneList, 1, function(r) any(r %in% c("26896"))))


# ------------------gseDO function-------- 
y <- gseDO(gene_num_vec,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.02,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y)
y_df <- y@result



# -------***----- enrich DGN (gene-disease associations)-----*****-----------
eDGN <- enrichDGN(gene)

eDGN_result <- eDGN@result


# -----****------- gse DGN function --------**********-----------------
gse_dgn <- gseDGN(gene_num_vec,  
                  nPerm         = 100,
                  minGSSize     = 120,
                  pvalueCutoff  = 0.2,
                  pAdjustMethod = "BH",
                  verbose       = FALSE)

gse_dgn <- setReadable(gse_dgn,'org.Hs.eg.db')
head(gse_dgn)










