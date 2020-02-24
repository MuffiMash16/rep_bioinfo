library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tidytext)

#---- 1. Load File ---
data_file <- readRDS("~/r_projects/mashood_ora/data/data14Features.rds")
colnames(rowData(data_file))

# --- 2. Select columns from metadata -----
analysis_data <-metadata(data_file)$limma[,"Sirt_vs_Mock",]%>%
  as.data.frame()%>%
  dplyr::select(effect,p)%>%
  cbind(.,ENTREZ_ID = rowData(data_file)$ENTREZID)

#--- 3. Clean ENTREZ_ID ------
analysis_data$ENTREZ_ID <- analysis_data$ENTREZ_ID %>%
  stringr::str_split(.,";",simplify=TRUE)%>%
  .[,1]

#--- 4. Grouping data file based on P & Effect,store Group ID in a new column---- (ROW 31)
analysis_data <- mutate(analysis_data, Grp_ID = ifelse(p <= 0.05 & effect <= 0,1,
                                                       ifelse(p <= 0.05 & effect >= 0,2,3)))

#--- 5. Load Wikipathway info-------
wiki <- read.gmt("~/r_projects/mashood_ora/data/wikipathways-20191210-gmt-Mus_musculus.gmt")%>%
  tidyr::separate(ont,c("name","version","wpid","org"),"%")

wiki_geneID <- wiki %>% dplyr::select(wpid,gene)
wiki_geneName <- wiki%>% dplyr::select(wpid,name)

#----6.Create Gene Universe 
gene_universe <- analysis_data %>% dplyr::filter(!is.na(Grp_ID) & ENTREZ_ID != "" )


#----7.Use a Loop to find enriched value for each group + fetch Log vals + filter based on log vals
#------unique () selects unique values and cols not available = false
no_of_grps <- analysis_data$Grp_ID %>% unique() %>% .[!is.na(.)]

ora_list <- no_of_grps %>% lapply(., function(x){
  gene_oi <- analysis_data %>% dplyr::filter(Grp_ID == x & ENTREZ_ID != "")
  enrichr_result <- enricher(gene_oi$ENTREZ_ID, TERM2GENE = wiki_geneID,TERM2NAME = wiki_geneName,universe = gene_universe$ENTREZ_ID)
  enrichr_result@result %>%
    dplyr::mutate(p_log = -1*log(pvalue, 10))%>%
    dplyr::filter(p_log > 1.3) -> enrchresult_Grp
   enrchresult_Grp %>% mutate(Group_ID = x)
})

#--- 8. Bind data frames returned from all groups  and save it as a data frame
#   do.call(rbind, .) only accepts lists. hence we convert it into a DF in the next line
ora_list_df <- ora_list %>% do.call(rbind, .)

#--- 9. D/F ordered by p_log
plot <- ora_list_df%>%
  # 1. Arrange by - i.  facet group & ii. bar height
  arrange(Group_ID, p_log) %>%
  # 2. Add order column of row numbers
  mutate(order = row_number())
View(plot)


#---10. Plotting -----------------------------
plot %>%
  ggplot(aes(x=order,y= p_log)) +
  geom_bar(stat = "identity", fill="salmon2", alpha=0.8)+
  facet_wrap(~Group_ID)+
  theme_bw()+
  # Add Description to x axis
  scale_x_continuous(
    breaks = plot$order,labels = plot$Description, c(0,0)
  )+
  coord_flip()

View(plot)


#----- Supporting/verification commands ---------
View()
count()
rowData()
colnames()












  