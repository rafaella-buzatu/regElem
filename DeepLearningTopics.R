library(keras)
library(writexl)
library(readxl)
source('utils/CNNtopicsFunctions.R')
library(dplyr)

#Define path to store plots
pathToPlotsDir <- 'plots/cellTypeFromTopic'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir <- 'outputs/cellTypeFromTopic'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

###    GET INPUT
#Read input file
#ATACseqData <- readRDS('data/ATACsubsetMedian.Rds')
#ATACseqData <- readRDS('data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
#metadata <- read.delim('data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#cellData <- read_xlsx('outputs/cisTopic/run3/cellData.xlsx')

#Create dataframe with read counts per type of cell for each region
#cellTypesFromTopics <- getCellTypeFromTopics ( metadata, cellData )
#write_xlsx (cellTypesFromTopics, file.path(pathToOutputsDir,'cellTypesFromTopics.xlsx'))
cellTypesFromTopics= read_xlsx('outputs/cellTypeFromTopic/cellTypesFromTopics.xlsx')


###   PROCESS INPUT FOR CNN
##Subset Glutamatergic cells because they are over-represented
set.seed(10)
#Get indices of glutamatergic cells
indicesGlutamatergic <- which(cellTypesFromTopics$cellType == 'Glutamatergic')
#Get max between number of NN and GB cells
keep = max (length(which(cellTypesFromTopics['cellType'] == 'GABAergic')),
            length(which(cellTypesFromTopics['cellType'] == 'Non-Neuronal')))
cellTypesFromTopics = cellTypesFromTopics[-sample(indicesGlutamatergic,
                                                  length(indicesGlutamatergic) - keep), ]


#Drop sample name column
cellTypesFromTopics <- cellTypesFromTopics[ , !(names(cellTypesFromTopics) %in% c('sampleName'))]
#Convert to tensors and split in training, test and validation sets
inputData = splitTrainTest (cellTypesFromTopics, testPercentage = 0.2, seed = 10)

#Extract train set
xTrain = inputData$train$x
yTrain = inputData$train$y
#Extract test set
xTest = inputData$test$x
yTest = inputData$test$y

#Get labels
labels = inputData$train$labels

#Create model
inputShape <- dim(xTrain)[2]
nClasses <- dim(yTrain)[2]
model <- createModelMultiLabelClassification (inputShape, nClasses)

#Train model
nnModel = trainModel(xTrain, yTrain,
                      model,
                      batchSize = 128, 
                      epochs =50, 
                      patience = 10,
                      valSplit = 0.2,
                      pathToPlotsDir = pathToPlotsDir)

#Load model
#nnModel <- load_model_hdf5(file.path(pathToOutputsDir,"nnModel.hdf5"))

#Evaluate model on test set
nnModel %>% evaluate(xTest, yTest)
#      loss   accuracy 
#0.03237376 0.99664336 

#Get predictions
yPred <- getPredictions (nnModel, xTest)
write_xlsx (yPred, file.path(pathToOutputsDir,'predictions.xlsx'))
yTest <-data.frame(argmax(yTest))
#Save model
save_model_hdf5(nnModel, file.path(pathToOutputsDir,"nnModel.hdf5"))
print (summary(nnModel))

#For Regression
plotConfusionMatrix (yPred, yTest, labels, pathToPlotsDir)



####  Get average scores of topics per cell type

cellTypesFromTopics= read_xlsx('outputs/cellTypeFromTopic/cellTypesFromTopics.xlsx')
#Extract names of cell types and topics
cellTypes = unlist(unique(cellTypesFromTopics['cellType']))
topics= colnames(cellTypesFromTopics [ , grepl( "Topic" , names( cellTypesFromTopics ) )])

#Create df to store averages
averages = data.frame(matrix(nrow = length(topics), ncol = length(cellTypes)))
colnames(averages) = unlist(cellTypes)
row.names(averages) = topics

#Fill in df
for (cellType in cellTypes){
  for (topic in topics){
    average = colMeans(cellTypesFromTopics[which(cellTypesFromTopics['cellType'] == cellType), topic])
    averages[topic, cellType] = average
  }
}

#Extract all topic scores per region
regionData <- read_xlsx('outputs/cisTopic/run3/regionData.xlsx')

#TOPIC DIFFERENCE PER CELL TYPE (from plots and averages):
#GB -> {7, 13}
#GL -> {4, 6, 10}

#GB Regions
#Extract highest 10 regions for each relevant topic
Topic7<- regionData[order(regionData$Topic7, decreasing = TRUE ), ][1:10, c('seqnames', 'start', 'end', 'Topic7')]
colnames(Topic7)[4] <- "TopicScore"
Topic13<- regionData[order(regionData$Topic13, decreasing = TRUE ), ][1:10, c('seqnames', 'start', 'end', 'Topic13')]
colnames(Topic13)[4] <- "TopicScore"
#Merge in one dataframe
topicGB <- rbind (Topic7, Topic13)
#Remove duplicated
topicGB <- topicGB[!duplicated(topicGB[ , c('seqnames', 'start', 'end')]),]

#GL Regions
#GB Regions
Topic4<- regionData[order(regionData$Topic4, decreasing = TRUE ), ][1:10, c('seqnames', 'start', 'end', 'Topic4')]
colnames(Topic4)[4] <- "TopicScore"
Topic6<- regionData[order(regionData$Topic6, decreasing = TRUE ), ][1:10, c('seqnames', 'start', 'end', 'Topic6')]
colnames(Topic6)[4] <- "TopicScore"
Topic10<- regionData[order(regionData$Topic10, decreasing = TRUE ), ][1:10, c('seqnames', 'start', 'end', 'Topic10')]
colnames(Topic10)[4] <- "TopicScore"
topicGL <- rbind (Topic4, Topic6, Topic10)
topicGL <- topicGL[!duplicated(topicGL[ , c('seqnames', 'start', 'end')]),]



