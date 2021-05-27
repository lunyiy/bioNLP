library(ggplot)
library(org.Hs.eg.db)
library(pathview)
library(clusterProfiler)

HPO.EnrichToFile <- function(GeneList, genesToPhenotypeToDisease, SaveFileName = "HPO") {
  
  print("HPO enrichment analysis begining, please be patient")
  GeneList$V1<-as.character(GeneList$V1)
  GeneList2=bitr(GeneList$V1,fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Hs.eg.db')

  n <- length(GeneList2$ENTREZID)

  hpoId <- unique(genesToPhenotypeToDisease$HPO.Term.ID)
  mValue2 <- data.frame(mValue = character())
  for(valueId in hpoId)
  {
    allgeneInCategory <- unique(genesToPhenotypeToDisease[genesToPhenotypeToDisease$HPO.Term.ID == valueId,])
    M <- nrow(allgeneInCategory)
    mValue <- c(M)
    M_Frame <- data.frame(mValue)
    mValue2 <- rbind(mValue2,M_Frame)
  }
  sumData <- cbind(hpoId,mValue2)

  genesToPhenotypeToDiseaseTotalGene <- unique(genesToPhenotypeToDisease$entrez.gene.id)
  N <- length(genesToPhenotypeToDiseaseTotalGene)

  print("Completed 25%")
  sampleEntrez <- unique(GeneList2$ENTREZID)

  genes <- data.frame(genes = character())
  kvalue2 <- data.frame(kValue = character())
  for(valueId in hpoId)
  {
    allgeneInCategory <- unique(genesToPhenotypeToDisease[genesToPhenotypeToDisease$HPO.Term.ID== valueId,])
    totalHpoGene <- allgeneInCategory$entrez.gene.id
    intersectGeneEntrezs <-intersect(sampleEntrez,totalHpoGene)
    richGeneSymbol = NULL
    if (length(intersectGeneEntrezs) != 0){
      for (intersect_gene_entrez in intersectGeneEntrezs){
        intersect_genes <- GeneList2[GeneList2$ENTREZID == intersect_gene_entrez,]
        intersect_gene <- intersect_genes$SYMBOL
        richGeneSymbol<- paste(richGeneSymbol,intersect_gene,seq = ",")}
    }else
      richGeneSymbol <- " "
    gene <- c(richGeneSymbol)
    gene <- data.frame(gene)
    genes <- rbind(genes,gene)
    k = sum(sampleEntrez %in% totalHpoGene)
    kValue = c(k)
    kFrame <- data.frame(kValue)
    kvalue2 <- rbind(kvalue2,kFrame)
  }
  sumData <- cbind(sumData,kvalue2)
  sumData <- cbind(sumData,genes)

  pValue2 <- data.frame(pValue = character())
  for(valueId in hpoId){
    pgeneInCategory <- sumData[sumData$hpoId== valueId,]
    k <- pgeneInCategory$kValue
    M <- pgeneInCategory$mValue
    p <- phyper(k-1,M, N-M, n, lower.tail=FALSE)
    pValue <- c(p)
    pFrame <- data.frame(pValue)
    pValue2 <- rbind(pValue2,pFrame)
  }
  print("Completed 50%")
  hpoId <- data.frame(hpoId)
  hpoData <- cbind(hpoId,pValue2)

  hpoData <- unique(hpoData[hpoData$pValue<0.05,])

  hpoData <- hpoData[order(as.numeric(as.character(hpoData$pValue))),]
  fdr <- fdrtool(hpoData$pValue,statistic = "pvalue")
  qValue <- fdr$qval
  qValue <- data.frame(qValue)
  hpoData <- cbind(hpoData,qValue)
  
  pHpoId <- unique(hpoData$hpoId)
  sumHpo <- data.frame(description = character())
  for(valueId in pHpoId){
    PgeneInCategory <- unique(genesToPhenotypeToDisease[genesToPhenotypeToDisease$HPO.Term.ID== valueId,])
    description <- unique(PgeneInCategory$HPO.Term.Name)
    description <- data.frame(description)
    sumHpo <- rbind(sumHpo,description)
  }
  gene_number <- data.frame(Count = character())
  gene.name <- data.frame(gene.name = character())
  for(valueId in pHpoId){
    numGeneInCategory <- sumData[sumData$hpoId== valueId,]
    k <- numGeneInCategory$kValue
    Count <- c(k)
    Count <- data.frame(Count)
    Gene_names <- numGeneInCategory$gene
    Gene_name <- data.frame(Gene_names)
    gene_number <- rbind(gene_number,Count)
    gene.name <- rbind(gene.name,Gene_name)
  }
  sumHpo <- cbind(sumHpo,gene_number,hpoData,gene.name)
  sumHpo <- sumHpo[order(as.numeric(as.character(sumHpo$pValue))),]
  print("Completed 75%")

  hpoFileName = paste(c(SaveFileName, "csv"), collapse = '.')
  write.csv(sumHpo,hpoFileName,row.names = FALSE)
  
  print("Completed 100%, success!")
  cat(hpoFileName,"file has been generated in the working directory, please note that check!\n")
}

HPO.FileToPlot <- function(hpoFileName = "HPO.csv") {
  if(substring(hpoFileName,nchar(hpoFileName)-3) != '.csv'){
    hpoFileName <- paste(c(hpoFileName, 'csv'), collapse = '.')
  }
  data <- read.csv(hpoFileName, header = TRUE)

  ggplot(data = data,aes(x =Description,
                         y = Count,
                         fill = -log10(pValue))) +
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=data$description) +
  coord_flip() + labs(title = "EnrichmentHPO") + 
  theme(plot.title = element_text(size = 15,face = "bold"),
        axis.text = element_text(size = 8,face = "bold"),
        axis.title.x =element_text(size=14),
        axis.title.y=element_text(size=16),
        panel.background = element_rect(fill="white", colour='gray')) +
  scale_fill_gradient(low = 'blue', high = 'red')
}

setwd("F:/Users/Lun/Desktop/NLP/")
genes_to_phenotype <- read.table("genes_to_phenotype.txt",
                                 header = TRUE,
                                 sep="\t",
                                 stringsAsFactors = FALSE,
                                 quote = "")

genesToPhenotypeToDisease <- genes_to_phenotype[,c(1,2,3,4,9)]
#colnames(genesToPhenotypeToDisease) <- c("entrez-gene-id",	"gene-symbol",	"HPO-Term-ID",	"HPO-Term-Name",	"disease-ID")

#Gene
up_regulated_gene <- read.table('up-regulated-gene.txt', header = FALSE)
down_regulated_gene <- read.table('down-regulated-gene.txt', header = FALSE)

HPO.EnrichToFile(up_regulated_gene, genesToPhenotypeToDisease, 'HPO-up')
HPO.FileToPlot("HPO-up")
HPO.EnrichToFile(down_regulated_gene, genesToPhenotypeToDisease, 'HPO-down')
HPO.FileToPlot("HPO-down")
