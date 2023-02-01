library(BSgenome.Hsapiens.UCSC.hg38)

createRegionDataFrame <- function (regionScoresAllTopics, regionScoresPerTopic)
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
  
  
  #Create dataframe by subsetting the input matrix with all regions
  regionData <- subset(regionScoresAllTopics, select=c('seqnames','start','end','width','nCounts','nCells'))
  
  #Define progress bar
  pb = txtProgressBar(min = 0, max = length(regionScoresPerTopic), style = 3, width = 50) 
  
  
  #Iterate over all topics
  for (topic in  1:length(regionScoresPerTopic)){
    
    #Create new column with zeros
    stringTopic <- c('Topic')
    columnName <- paste (stringTopic,topic, sep = '')
    regionData[columnName] <- integer (nrow(regionData))
    
    #Iterate over all regions in the topic
    for (region in row.names(regionScoresPerTopic[[topic]])){
      #Mark region with '1' if it belongs to given topic
      regionData[region, c(columnName)] = 1
  
    }
    setTxtProgressBar(pb, topic)
  }
  
  close(pb)
  #Return new df
  return (regionData)
}

addDNAsequences <- function (regionData){
  #' Adds the corresponding DNA sequences to each row of the input dataframe
  
  regionData ['DNAseq'] <- NA
  #Define progress bar
  pb = txtProgressBar(min = 0, max = length(regionData), style = 3, width = 50) 
  
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