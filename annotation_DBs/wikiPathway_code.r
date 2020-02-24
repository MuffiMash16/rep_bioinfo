library(magrittr)
library(clusterProfiler)

data(geneList, package= "DOSE")
#gene <- geneList[abs(geneList) > 2]
#View(gene)
#Print the names attribute from the geneList dataset, where abs. value is > 2
gene <- names(geneList)[abs(geneList) > 2]
View(gene)

#Load file from WIKIPATHWAY
wikiPathwayData <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
gene_data<- read.gmt(wikiPathwayData)
head(gene_data)

#Seperate 1 column into multiple columns using seperate(col_name,new_cols, seperator symbol)
multi_col_gene_data <-  gene_data %>% tidyr::separate(ont,c("name","version","wpid","org"),"%")

#Select the required columns for Gene ID
gene_id <- multi_col_gene_data %>% dplyr::select(wpid,gene)

#Select the required columns for Gene Name
gene_name <- multi_col_gene_data %>% dplyr::select(wpid,name)

#enricher function - 3 arguments passed
enricher_output <- enricher(gene, TERM2GENE = gene_id, TERM2NAME = gene_name)
head(enricher_output)
View(enricher_output)

#GSEA function - geneList contains ID and names from preloaded data 
GSEA_output <-GSEA(geneList, TERM2GENE = gene_id ,TERM2NAME = gene_name,verbose = FALSE)
head(GSEA_output)

#Convert Gene ID to Gene Symbols using setReadable()
library(org.Hs.eg.db)

enricher_output_symbol <- setReadable(enricher_output,org.Hs.eg.db,keyType = "ENTREZID")
head(enricher_output_symbol)

GSEA_output_symbol <- setReadable(GSEA_output,org.Hs.eg.db,keyType = "ENTREZID")
head(GSEA_output_symbol)




