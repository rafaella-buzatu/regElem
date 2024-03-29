### CISTOPIC

library(BSgenome.Hsapiens.UCSC.hg38)

createRegionDataFrame <- function (regionScoresAllTopics, regionScoresPerTopic, 
                                   scoreType = 'all')
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
  #'   all = all numerical scores (DEFAULT)
  
  #Get number of topics
  topicNo = ncol(regionScoresAllTopics[ , grepl( "Topic" , names( regionScoresAllTopics ) )])
  
  #If we want to keep all scores, we simply rename the columns of the input df
  if (scoreType == 'all'){
    
    regionData = regionScoresAllTopics
    columns = c('seqnames','start','end', 'width','nCounts','nCells')
    for (i in 1:topicNo){
      columns = append (columns,  paste(c('Topic', i), collapse = ""))
    }
    colnames(regionData) = columns
    
  } else {
  #Otherwise, we create dataframe by subsetting the input df with all topics
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


getBaseWindow <- function (regionData, windowSize) {
  #' Extracts a  window from each region in the input dataframe, where
  #' the middle is the center of the initial region.
  
  pb = txtProgressBar(min = 0, max = nrow(regionData), style = 3, width = 50) 
  for (row in 1:nrow(regionData)) {
    
    if(regionData$width[row] %% 2 == 1) {middle = (regionData$width[row]-1)/2}
    if(regionData$width[row] %% 2 == 0) {middle = regionData$width[row]/2}
    
    newMid = regionData$start[row] + middle
    regionData$start[row] = newMid - (windowSize/2 -1)
    regionData$end[row] = newMid + (windowSize/2)
    regionData$width[row] = regionData$end [row] - regionData$start[row] + 1
    
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  return (regionData)
}

createCellTopicDataFrame <- function (cellTopicAssignments, dimReductionCoords , metadata) {

  #' Generates a dataframe combining the information about the topic assignment
  #' scores for each input cell with the UMAP scores for plotting.
  #'
  #' cellTopicAssignments = A dataframe containing the topic scores for each cell.
  #' 
  #' dimReductionCoords  = A matrix containing the UMAP coordinates of the cell-topic assignments
  #' 
  #' metadata = A dataframe with metadata for each cell
  
  #Get number of topics
  topicNo = nrow(cellTopicAssignments)
  
  #Extract tSNE and UMAP coordinates
  tSNEscores = dimReductionCoords$tSNE
  UMAPscores = dimReductionCoords$Umap
  
  # Create a Vector with the column names
  columns = c('sampleName', 'UMAP1', 'UMAP2','tSNE1', 'tSNE2', 'cellType')
  
  for (i in 1:topicNo){
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
    for (topic in 1:topicNo) { 
      cellData[row, paste (c('Topic'),topic, sep = '')] = cellTopicAssignments [topic, row]
    }
    #Add UMAP coordinates 
    cellData[row, 'UMAP1'] = UMAPscores[cellData[row, 'sampleName'], 'UMAP1']
    cellData[row, 'UMAP2'] = UMAPscores[cellData[row, 'sampleName'], 'UMAP2']
    
    #Add tSNE coordinates 
    cellData[row, 'tSNE1'] = tSNEscores[cellData[row, 'sampleName'], 'tSNE1']
    cellData[row, 'tSNE2'] = tSNEscores[cellData[row, 'sampleName'], 'tSNE2']
    
    #Add cell type
    cellData[row, 'cellType'] = metadata['subclass'][metadata['sample_name'] == cellData[row, 'sampleName']]
  }
  return (cellData)
}

getCellTypePerRegion  <- function (ATACseqData, metadata, countType) {
  #' Generates a dataframe with ATAC-seq read counts per cell Type for each DNA region.
  
  #Create empty dataframe with the corresponding columns to specify the location
  #of each region, followed by the number of reads in each type of cell
  columns = c ('seqnames', 'start', 'end', 'width')
  cellTypes = unique( unlist(metadata['level1']))
  columns = append(columns, cellTypes)
  cellTypePerRegion = data.frame(matrix(nrow = nrow(ATACseqData), ncol = length(columns))) 
  colnames(cellTypePerRegion) = columns
  
  #Fill with zeros
  cellTypePerRegion[is.na(cellTypePerRegion)] <- 0
  
  
  pb = txtProgressBar(min = 0, max = nrow(ATACseqData), style = 3, width = 50) 
  for (row in 1:nrow(ATACseqData)){
    
    #Add location information to corresponding columns
    location = row.names(ATACseqData)[row]
    cellTypePerRegion [row, 'seqnames'] = strsplit(location, ":")[[1]][1]
    cellTypePerRegion [row, 'start'] = as.numeric(strsplit(strsplit(location, ":")[[1]][2], '-')[[1]][1])
    cellTypePerRegion [row, 'end'] = as.numeric(strsplit(strsplit(location, ":")[[1]][2], '-')[[1]][2])
    cellTypePerRegion [row, 'width'] = cellTypePerRegion[row, 'end'] - cellTypePerRegion [ row,'start'] +1
  
    setTxtProgressBar(pb, row) 
  }
  close(pb)
 
  #Add peak counts per cell type
   for (cellType in cellTypes) {
    #Get index for columns of a given type
    cellsOfAType = metadata[metadata$level1 == cellType, 'sample_name']
    if (countType == 'reads'){
      #Count the reads
      cellTypeSums = rowSums(ATACseqData [, cellsOfAType]) }
    else if (countType == 'cells') {
      #Count the cells
      cellTypeSums = rowSums(ATACseqData [, cellsOfAType]>0) }
    
    cellTypePerRegion [, cellType] = cellTypeSums
    
  }
  
  #Remove chromosomes X and Y
  cellTypePerRegion <- subset (cellTypePerRegion, cellTypePerRegion$seqnames!= 'chrX' & cellTypePerRegion$seqnames!= 'chrY')
  
  return (cellTypePerRegion)
}

