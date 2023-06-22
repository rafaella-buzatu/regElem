library(dplyr)
library(openxlsx)


getGeneTopicScores <- function (regionData, linkedPeaks){
  
  peaks = c()
  geneNames = c()
  
  print ('Extracting genes from peaks...')
  #Define progress bar
  pb = txtProgressBar(min = 0, max = nrow(regionData), style = 3, width = 50) 
  
  #Extract peak and gene name combinations to keep
  for (row in 1:nrow(regionData)){
    #Extract peak position in linkedPeaks format
    peakPosition = paste (regionData [row, 'seqnames'], '-',
                          regionData [row, 'start'], '-',
                          regionData [row, 'end'], sep = "")
    #Extract genes assigned to peak of interest
    genes = linkedPeaks [which(linkedPeaks$peak == peakPosition), c( 'gene', 'pvalue')]
   
    if (nrow(genes) >0) {
      #Append positions in df and gene names to the corresponding lists
      for (i in 1:nrow(genes)){
        peaks = rep(append (peaks, row))
        }
      geneNames = append(geneNames, unlist(genes['gene']))
    }
    
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  #Get topic names
  topics = colnames(regionData[ , grepl( "Topic" , names( regionData))])
  #Create list to save vector scores
  TopicGeneScoresList <- vector("list", length( topics))
  #Name the elements
  names(TopicGeneScoresList) <- topics
  
  #Only look at the regions associated with genes
  regionDataSubset = regionData[peaks,]
  #Add gene nemes as variables
  regionDataSubset['geneName'] = geneNames
  
  print ('Creating topic score vectors...')
  
  for (topic in topics){
    topicScores = regionDataSubset [, c(topic, 'geneName')]
    #Remove columns where topic score is 0
    topicScores = topicScores[!(topicScores[, topic] == 0), ]
    #Order by topic score
    topicScores <- dplyr::arrange(topicScores, desc(topicScores[topic]))
    #Remove gene duplicates so as to keep only the highest topic score
    topicScoresNoDup = topicScores[!duplicated(topicScores['geneName']),]
    
    #Convert to vector
    vectorTopicScores = unlist(topicScoresNoDup[topic])
    #Add gene names to vector
    names(vectorTopicScores) = unlist(topicScoresNoDup['geneName'])
    
    #Add to list of vectors
    TopicGeneScoresList[[topic]] <- vectorTopicScores
    
  }  
  
  return (TopicGeneScoresList)
}


readGMT <- function(path, delim = "\t", store_metadata = FALSE) {
  gmt <- readLines(path)
  gmt <- strsplit(gmt, "\t")
  metadata <- unlist(lapply(gmt, FUN = function(x) x[2]))
  pathways <- lapply(gmt, FUN = function(x) toupper(x[3:length(x)]))
  names(pathways) <-  unlist(lapply(gmt, FUN = function(x) x[1]))
  names(metadata) <- names(pathways)
  if(!store_metadata) pathways else list(geneSets = pathways, metadata = metadata)
}

geneGOmtx = function(geneSets, pvals, min_gene = 5, max_gene = 500){
  # Optional min_gene & max_gene: Keep only GO terms with > 5 and < 500 genes
  mtx = lapply(1:length(geneSets),
               FUN = function(i, data, gene_sets) {
                 x <- data.frame(names(data) %in% as.character(unlist(gene_sets[i])))
                 names(x) <- names(gene_sets[i])
                 x
               },
               data = pvals,
               gene_sets = geneSets)
  
  mtx = dplyr::bind_cols(mtx)
  sums_matrix = colSums(mtx)
  keep_index = between(sums_matrix, min_gene, max_gene)
  keep_index = which(keep_index == TRUE)
  keep_matrix = mtx[,keep_index]
  rownames(keep_matrix) = names(pvals)
  keep_matrix
}

GSEA = function(gene_go_mtx, pvals){
  gsea = data.frame(matrix(ncol = 4, nrow = ncol(gene_go_mtx)))
  colnames(gsea) = c('GO_term','p', 'coef', 'p.adj')
  gsea[,1] = colnames(gene_go_mtx)
  for (i in 1:ncol(gene_go_mtx)) {
    res = lm(pvals ~ gene_go_mtx[,i]  )
    gsea[i,"p"] = coef(summary(res))[2,4] #res$prob[2]
    gsea[i,"coef"] = coef(summary(res))[2,1] #res$coefficients[2]
  }
  gsea$p.adj = p.adjust(gsea$p, method = "fdr")
  gsea[order(gsea$p.adj),]
}


runGSEAallTopics = function (topicScoreGenes, geneSets, pathToOutputsDir){
  
  #Get topics
  topics <- names(topicScoreGenes)
  
  wb <- createWorkbook()
  fileName = file.path (pathToOutputsDir, 'GSEAresults.xlsx')
 
  for (topic in topics){
  
  print (paste ('Running GSEA on', topic, '...'))
  #Extract corresponding topic scores
  topicScores = topicScoreGenes[[topic]]
  #Create input data 
  gomtx_nygc = geneGOmtx(geneSets, topicScores)
  #Run GSEA analysis
  gsea_nygc = GSEA(gomtx_nygc, topicScores)
  
  #Write results to workbook
  addWorksheet(wb, sheetName=topic)
  writeData(wb, sheet = topic, gsea_nygc)
  
  }
  
  saveWorkbook(wb, file=fileName, overwrite = TRUE)

}