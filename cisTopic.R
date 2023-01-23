library (cisTopic)
source('altered_functions.R')

#Load ATACseq data
ATACseq_data = readRDS('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')

subset = ATACseq_data[1:2000,1:2000]

#Transform matrix into cisTopic object
cisTopicObject <- createcisTopicObject(subset, project.name='ATACseq_clustering')

#Add metadata
#metadata <- read.delim('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')

#metadata = metadata [1:2000, ]

#cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)

#Run models with chosen topic numbers
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=1, addModels=FALSE)

#Select best model
cisTopicObject <- selectModelModified(cisTopicObject, type='derivative')

### CELL TYPES
#Run Umap
cisTopicObject <- runUmap(cisTopicObject, target='cell')

#Create Umap Plot
pdf('Umap.pdf')
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL)#, colorBy=c('subClass', 'level1'))
dev.off ()

pdf('heatmap.pdf')
cellTopicHeatmap(cisTopicObject, method='Probability')#, colorBy=c('subclass'))
dev.off()

### REGULATORY REGIONS
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)


#Save the cisTopicObject
saveRDS(cisTopicObject, file='cisTopicObject.Rds')