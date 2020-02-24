# Example from yulab
library(clusterProfiler)

# Dataset
data("geneList",package = "DOSE")

gene <- names(geneList)[abs(geneList) > 2]

# gene df is used in enrichr ()
gene.df <- bitr(gene, fromType = "ENTREZID",
            toType = c("ENSEMBL","SYMBOL"),
            OrgDb = org.Hs.eg.db)
head(gene.df)

# ----- GO Classification - groupGO function - Classifies genes based on GO distribution at a specific level.
ggo <- groupGO(gene = gene,
               OrgDb    = org.Hs.eg.db,
               ont = "CC",
               level = 3,
               readable = TRUE
               )

head(ggo)

# ----- GO Over Representation Analysis-------
#------ 1- Over Representation Analysis
e_go <- enrichGO(gene = gene, 
                 universe = names(geneList),
                 OrgDb = org.Hs.eg.db, 
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05, 
                 readable = TRUE)
head(ego)

# ----- 2- Over Representation Analysis,specify 'keyType'
e_go_2 <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL",
                   ont= "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

# Map Gene ID to Gene Symbol using setreadable() or readble = TRUE
e_go_2 <- setReadable(e_go_2, OrgDb = org.Hs.eg.db)

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
# ----- *** ----- GO gsea -----*** -----

gsea_go <- gseGO(geneList = geneList,
                 OrgDb = org.Hs.eg.db,
                 ont = "CC",
                 nPerm = 1000,
                 minGSSize = 100 ,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE )


#-------- GO Semantic Similarity Analysis ( Not sure about the values passed) ---
# use it to cluster genes/proteins into different clusters based on their 
# functional similarity and to measure the similarities among GO terms to 
# reduce the redundancy of GO enrichment results.
go_similarity <- GOSemSim::godata(OrgDb = org.Hs.eg.db,
                 keytype = "ENSEMBL",
                 ont = "CC",
                 computeIC = TRUE)


