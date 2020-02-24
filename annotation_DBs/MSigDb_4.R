#MSigDbr data is fetched from "CLUSTER PROFILER" instead of MSigDbr package
# Refer Chapter 7 from the Book

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

