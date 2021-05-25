library('devtools')
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

keytypes(org.Hs.eg.db)

MyGeneSet2 <- read.table("G:/R/data/nlp/03/GeneNames.csv",header=FALSE)
MyGeneSet2$V1 <- as.character(MyGeneSet2$V1) 
geneid <- MyGeneSet2$Geneid

MyGeneIDSet2 = bitr(MyGeneSet2$V1,
                    fromType="SYMBOL",
                    toType=c("ENSEMBL","ENTREZID", "GO"),
                    OrgDb="org.Hs.eg.db")
head(MyGeneIDSet2,2)

ego_ALL <- enrichGO(gene = MyGeneIDSet2$ENTREZID, 
                    universe = names(geneList),
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01, 
                    qvalueCutoff = 0.05,
                    readable = TRUE)

write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)

dotplot(ego_ALL,title="EnrichmentGO_ALL_dot")
barplot(ego_ALL, showCategory=20,title="EnrichmentGO_ALL")

ego_MF <- enrichGO(gene = MyGeneIDSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
plotGOgraph(ego_MF,firstSigNodes = 10, useInfo = "all", sigForAll = TRUE,
            useFullNames = TRUE)

ego_CC <- enrichGO(gene = MyGeneIDSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
plotGOgraph(ego_CC,firstSigNodes = 10, useInfo = "all", sigForAll = TRUE,
            useFullNames = TRUE)

ego_BP <- enrichGO(gene = MyGeneIDSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
plotGOgraph(ego_BP,firstSigNodes = 10, useInfo = "all", sigForAll = TRUE,
            useFullNames = TRUE)