#From yulab book - 

library(DOSE)

#---- Dataset from DOSE --------------
data(geneList)

gene <- names(geneList)[abs(geneList) > 1.5]
head(gene)

#--------------- DO (Disease Ontology) -------------------------------
# enrichDo function - is useful for identifying disease association of interesting genes
x <- enrichDO(gene          = gene,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
x_df <- x@result

# gseDO function- 
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y)
y_df <- y@result

# ------------ NCG (National Cancer Gene) ------------------
gene_2 <- names(geneList)[abs(geneList) < 3]
# enrichNCG 
ncg <- enrichNCG(gene_2)
head(ncg)

# gseNCG 
gse_ncg <- gseNCG(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)

#Fecth gene symbols
gse_ncg <- setReadable(ncg, 'org.Hs.eg.db')

gse_ncg_result <- gse_ncg@result

# ------------ enrichDGN (gene-disease associations)-------
eDGN <- enrichDGN(gene)

eDGN_result <- eDGN@result


# ------------ DGNv (snp -gene-disease associations)------------

snp <- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392",
         "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396",
         "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646",
         "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221") 

dgnv <- enrichDGNv(snp)
dgnv_result <- dgnv@result

# ------------ gseDGN function -------------------------
gse_dgn <- gseDGN(geneList,  nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)

gse_dgn <- setReadable(gse_dgn,'org.Hs.eg.db')
head(gse_dgn)










