#Load packages
library (cisTopic)
library(writexl)
library(readxl)
source('utils/altered_functions.R')
source('utils/data.R')
source('utils/plots.R')

#Define path to store plots
pathToPlotsDir = 'plots/cisTopic/run3'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir = 'outputs/cisTopic/run3'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

#Load ATACseq data
ATACseqData = readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
#ATACseqData = ATACseqData[, 1:2000]

#Transform matrix into cisTopic object
cisTopicObject <- createcisTopicObject(ATACseqData, project.name='ATACseq_clustering')

#Read metadata from fileATAC
metadataOrig <- read.delim('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#metadataOrig = metadataOrig [1:2000, ]

#Give object and dataframe the same row names
metadata <- data.frame(metadataOrig, row.names = 1)
#Add metadata to cisTopic Object from the dataframe
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)

#Run models with chosen topic numbers
#Run1
#cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 3, 5, 7, 10, 20, 30, 40, 50), seed=123, nCores=1, addModels=FALSE)
#Run2
#cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:10), seed=123, nCores=1, addModels=FALSE)
#Run3
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(5, 9, 10:20), seed=123, nCores=1, addModels=FALSE)

#Select best model
cisTopicObject <- selectModelModified(cisTopicObject, pathFolder = pathToPlotsDir, type='derivative')



######   CELL TYPES
#Run Umap
cisTopicObject <- runUmap(cisTopicObject, target='cell')
#Run tSNE
cisTopicObject <- runtSNE(cisTopicObject, target='cell')

#Extract Umap and tSNE scores for plotting
dimReductionCoords = cisTopicObject@dr$cell
#Read from Rds
#dimReductionCoords = readRDS(file.path(pathToOutputsDir,'dimReductionCoords.Rds'))

#Extract scores for topic assignment per cell
cellTopicAssignments <- cisTopicObject@selected.model$document_expects 
#read from Rds
#cellTopicAssignments <- readRDS(file.path(pathToOutputsDir,'cellTopicAssignments.Rds'))
#Convert to dataframe with all info
cellData <- createCellTopicDataFrame (cellTopicAssignments, dimReductionCoords, metadataOrig)

#Save dataframe
write_xlsx (cellData, file.path(pathToOutputsDir,'cellData.xlsx'))
#Read df from excel
#cellData = read_xlsx(file.path(pathToOutputsDir,'cellData.xlsx'))

#Create plots
createUMAPsPerTopic(cellData, pathToPlotsDir)
createtSNEsPerTopic(cellData, pathToPlotsDir)


###### REGULATORY REGIONS

#Get scores showing likelihood of region to belong to a topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

#Extract dataframe with scores for all regions for all topics
regionScoresAllTopics = cisTopicObject@region.data
#Read from Rds
#regionScoresAllTopics = readRDS(file.path(pathToOutputsDir,'regionScoresAllTopics.Rds'))

#Binarize topics; get most representative regions per topic
cisTopicObject <- binarizecisTopicsModified(cisTopicObject, pathFolder = pathToPlotsDir, thrP=0.975)

#Get regions for each topic and corresponding scores
regionScoresPerTopic = cisTopicObject@binarized.cisTopics
#Read from Rds
#regionScoresPerTopic = readRDS(file.path(pathToOutputsDir,'regionScoresPerTopic.Rds'))

#Combine information in final dataframe
regionData <- createRegionDataFrame (regionScoresAllTopics, regionScoresPerTopic )
#Subset
#regionData <- regionData [1:2000,]
#Save dataframe
write_xlsx (regionData, file.path(pathToOutputsDir,'regionData.xlsx'))
#Read df from excel
#regionData = read_xlsx(file.path(pathToOutputsDir,'regionData.xlsx'))

# RESULTS VISUALISATION

#Run tSNE 
cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)

#ANNOTATION
#Annotate GO terms to topics
library(org.Hs.eg.db)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

#Plot heatmap
pdf(file.path(pathToPlotsDir, 'GOsignaturesPerTopic.pdf'))
par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

#Save the cisTopicObject
saveRDS(cisTopicObject, file= file.path(pathToOutputsDir,'cisTopicObject.Rds'))
