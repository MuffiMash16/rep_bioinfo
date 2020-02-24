library(clusterProfiler)
 
#Search by KEGG code or Scientific name
search_kegg_organism("ece",by="kegg_code")

ecoli <- search_kegg_organism('Escherichia coli',by = 'scientific_name')
dim(ecoli)
head(ecoli)

# --- KEGG over Representation Test ------
data(geneList,package = "DOSE")
gene <- names(geneList)[abs(geneList) > 2]

#change organism to mouse for local dataset
kk_enrich <- enrichKEGG(gene=gene,
           organism = 'hsa',
           pvalueCutoff = 0.05
           )
head(kk_enrich)

# --- KEGG gene set enrichment analysis -----
kk_gsea <- gseKEGG(geneList = geneList,
               organism = 'hsa',
               nPerm = 1000,
               minGSSize = 120,
               pvalueCutoff = 0.05,
               verbose = FALSE)
head(kk_gsea)

# --- KEGG Module - over representation ---
mkk_enrich <- enrichKEGG(gene= gene,
           organism = 'hsa',
           )

# --- KEGG Module - gene set enrichment analaysis ----
mkk_gsea <- gseMKEGG(geneList = geneList, organism = 'hsa')


