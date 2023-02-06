#Load packages
library (cisTopic)
library(writexl)
library(readxl)
source('utils/altered_functions.R')
source('utils/data.R')
source('utils/plots.R')


#Define path to store plots
pathToPlotsDir = 'plots'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir = 'outputs'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputDir))
}

#Load ATACseq data
ATACseqData = readRDS('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
#ATACseqData = ATACseqData[, 1:2000]

#Transform matrix into cisTopic object
cisTopicObject <- createcisTopicObject(ATACseqData, project.name='ATACseq_clustering')

#Read metadata from fileATAC
metadata <- read.delim('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#metadata = metadata [1:2000, ]

#Give object and dataframe the same row names
metadata <- data.frame(metadata, row.names = 1)
#Add metadata to cisTopic Object
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)

#Run models with chosen topic numbers
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=1, addModels=FALSE)

#Select best model
cisTopicObject <- selectModelModified(cisTopicObject, pathFolder = pathToPlotsDir, type='derivative')

### CELL TYPES
#Run Umap
cisTopicObject <- runUmap(cisTopicObject, target='cell')

#Create Umap Plot
pdf(file.path(pathToPlotsDir, 'Umap.pdf'))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('subClass'))
dev.off ()

#Extract Umap scores for plotting
#UMAPscores = cisTopicObject@dr$cell
#Extract scores for topic assignment per cell
cellTopicAssignments <- cisTopicObject@selected.model$document_expects 
#read from Rds
cellTopicAssignments <- readRDS(file.path(pathToOutputsDir,'cellTopicAssignments.Rds'))
#Convert to dataframe with all info
cellData <- createCellTopicDataFrame (cellTopicAssignments, UMAPscores)

#Create plot



### REGULATORY REGIONS
#Get scores showing likelihood of region to belong to a topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)


#GATHER INFORMATION
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
regionDataTopics <- createRegionDataFrame (regionScoresAllTopics, regionScoresPerTopic)
#Add DNA sequences for each region
regionData <- addDNAsequences (regionDataTopics)
#Save dataframe
write_xlsx (regionData, file.path(pathToOutputsDir,'regionData.xlsx'))
#Read df from excel
regionData = read_xlsx(file.path(pathToOutputsDir,'regionData.xlsx'))

# RESULTS VISUALISATION
#Export region sets to bed files
getBedFiles(cisTopicObject, path=file.path(pathToOutputsDir,'cisTopicsBed'))

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
