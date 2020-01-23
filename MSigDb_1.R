# Code from d book - Homo Sapiens
# msigdbr library 
library(msigdbr)
library(magrittr)
library(clusterProfiler)

msigdbr_show_species()
#Retrieve all gene sets related to a mouse from msigdbr
#m_df <- msigdbr("Mus musculus") 

#All genes related to Humans - Cat 6 
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)

head(m_t2g)

#use existing dataset from DOSE
data(geneList, package="DOSE")

gene <- names(geneList)[abs(geneList) > 2]

em <- enricher(gene, TERM2GENE=m_t2g)
head(em)

em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em2)


#--- code to check datsets ----
head(geneList)
head(gene)

result <- em@result %>% dplyr::filter(ID == "RPS14_DN.V1_DN")
result

result2 <- em2@result %>% dplyr::filter(ID== "VEGF_A_UP.V1_DN")
result2
