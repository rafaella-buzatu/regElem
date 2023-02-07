library(BSgenome.Hsapiens.UCSC.hg38)

createRegionDataFrame <- function (regionScoresAllTopics, regionScoresPerTopic, 
                                   scoreType = 'threshold')
{
  #' Generates a dataframe combining the information from the input variables.
  #' 
  #' regionScoresAllTopics = dataframe containing information about each region,
  #'    such as chromosome, start and end location, the score for each topic, 
  #'    the number of number of counts mapped to that region accross the data set 
  #'    (nCounts),number of cells in which the region is accessible (nCells)
  #'
  #' regionScoresPerTopic = list of dataframes containing the scores for the each
  #'    topic's regions
  #'    
  #' scoreType = string specifying the format in which to keep the topic scores 
  #'             per region:
  #'   binary = binary above threshold,
  #'   threshold = numerical score above threshold, 
  #'   all = all numerical scores 
  
  
  #If we want to keep all scores, we simply rename the columns of the input df
  if (scoreType == 'all'){
    
    regionData = regionScoresAllTopics
    columns = c('seqnames','start','end', 'width','nCounts','nCells')
    for (i in 1:9){
      columns = append (columns,  paste(c('Topic', i), collapse = ""))
    }
    colnames(regionData) = columns
    
  } else {
  #Otehrwise, we create dataframe by subsetting the input df with all topics
  regionData <- subset(regionScoresAllTopics, select=c('seqnames','start','end',
                                                       'width','nCounts','nCells'))
  
  #Define progress bar
  pb = txtProgressBar(min = 0, max = length(regionScoresPerTopic), style = 3, width = 50) 
  
  #Iterate over all topics
  for (topic in  1:length(regionScoresPerTopic)){
    
    #Create new column with zeros
    columnName <- paste (c('Topic'),topic, sep = '')
    regionData[columnName] <- integer (nrow(regionData))
    
    #Iterate over all regions in the topic
    for (region in row.names(regionScoresPerTopic[[topic]])){
      
      if (scoreType == 'binary') {
      #Mark region with '1' if it belongs to given topic
      regionData[region, c(columnName)] = 1
      }
      
     else if (scoreType == 'threshold' ){
     #Add score for regions belonging to a certain topic
     regionData[region, c(columnName)] = regionScoresPerTopic[[topic]][region, 1]
     }
    }
    setTxtProgressBar(pb, topic)
  }
  
  close(pb)
  }
  #Return new df
  return (regionData)
}

addDNAsequences <- function (regionData){
  #' Adds the corresponding DNA sequences to each row of the input dataframe
  #' 
  #' regionData = A dataframe with information about each DNA region.
  
  regionData ['DNAseq'] <- NA
  #Define progress bar
  pb = txtProgressBar(min = 0, max = nrow(regionData), style = 3, width = 50) 
  
  for (row in 1:nrow(regionData)) {
  dnaSeq <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                            unlist(regionData[row, 'seqnames']), 
                                            unlist(regionData[row, 'start']), 
                                            unlist(regionData[row, 'end']) ))
  regionData[row, 'DNAseq'] = dnaSeq
  setTxtProgressBar(pb, row)
  }
  close(pb)
  return (regionData)
}

createCellTopicDataFrame <- function (cellTopicAssignments, UMAPscores, metadata) {

  #' Generates a dataframe combining the information about the topic assignment
  #' scores for each input cell with the UMAP scores for plotting.
  #'
  #' cellTopicAssignments = A dataframe containing the topic scores for each cell.
  #' 
  #' UMAPscores = A matrix containing the UMAP coordinates of the cell-topic assignments
  #' 
  #' metadata = A dataframe with metadata for each cell
  
  # Create a Vector with the column names
  columns = c('sampleName', 'UMAP1', 'UMAP2', 'cellType')
  for (i in 1:9){
    columns = append (columns, paste(c('Topic', i), collapse = ""))
  }
 
  #Create an empty dataframe
  cellData = data.frame(matrix(nrow = ncol(cellTopicAssignments), ncol = length(columns))) 
  # Assign column names
  colnames(cellData) = columns
  #Add cell names as column to new df
  cellData [, 'sampleName'] = colnames(cellTopicAssignments)
  
  #Iterate over all cells
  for (row in 1:nrow(cellData)) {
    #Add topic scores for each cell entry
    for (topic in 1:9) { 
      cellData[row, paste (c('Topic'),topic, sep = '')] = cellTopicAssignments [topic, row]
    }
    #Add UMAP coordinates 
    cellData[row, 'UMAP1'] = UMAPscores[cellData[row, 'sampleName'], 'UMAP1']
    cellData[row, 'UMAP2'] = UMAPscores[cellData[row, 'sampleName'], 'UMAP2']
    
    #Add cell type
    cellData[row, 'cellType'] = metadata['subclass'][metadata['sample_name'] == cellData[row, 'sampleName']]
  }
  return (cellData)
}