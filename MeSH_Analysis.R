# MeSH analysis using Local Dataset
# Data source from "gendoo", "gene2pubmed" and "RBBH" are all supported. User can selecte 
# interesting category to test. All 16 categories are supported. The analysis 
# supports >70 species listed in MeSHDb BiocView.
library(meshes)
library(MeSH.db)
library(MeSH.Hsa.eg.db)
library(enrichplot)

# Dataset from DOSE
data(geneList, package="DOSE")
de <- names(geneList)[1:100]

# --- enricher () ---------
# --- over-representation analysis, we use data source from gendoo and C (Diseases) category.
x <- enrichMeSH(de, 
                MeSHDb = 'MeSH.Hsa.eg.db', 
                database='gendoo', 
                category = 'C')
x_result <- x@result

#----- gsea() -------
#-- In the following example,  data source from gene2pubmed and 
# --test category G (Phenomena and Processes) using GSEA
y <- gseMeSH(geneList, 
             MeSHDb = "MeSH.Hsa.eg.db", 
             database = 'gene2pubmed', 
             category = "G",
             by = "DOSE")

y_result <- y@result

#---- Visualization ----
#--- enrichplot (i.e.barplot, dotplot, cnetplot, emapplot and gseaplot) to visualize these enrichment results ---
dotplot(x)
barplot(x)
gseaplot(y,y[1,1],title = y[1,2])


