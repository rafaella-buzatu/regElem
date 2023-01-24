#Load packages
library (cisTopic)
source('altered_functions.R')

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
#ATACseq_data = ATACseq_data[, 1:2000]

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
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('subClass', 'level1'))
dev.off ()

pdf(file.path(pathToPlotsDir, 'heatmap.pdf'))
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('subclass'))
dev.off()

### REGULATORY REGIONS
#Get scores showing likelihood of region to belong to a topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

#Generate BigWig files to observe scored regions on the genome
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
getBigwigFiles(cisTopicObject, path= file.path(pathToOutputsDir, 'cisTopicBigWig'), seqlengths=seqlengths(txdb))

#Save the cisTopicObject
saveRDS(cisTopicObject, file= file.path(pathToOutputsDir,'cisTopicObject.Rds'))