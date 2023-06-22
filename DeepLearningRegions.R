library(keras)
library(writexl)
library(readxl)
source('utils/plots.R')
source('utils/CNNregionsFunctions.R')
source('utils/data.R')

#Define path to store plots
pathToPlotsDir = 'plots/topicPrediction'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir = 'outputs/topicPrediction'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

#Read input data
regionData = read_xlsx(file.path('outputs', 'cisTopic','run3', 'regionData.xlsx'))

#Remove chromosomes X and Y
#regionData <- subset (regionData, regionData$seqnames!= 'chrX' & regionData!= 'chrY')

#Subset each region to 500 bases
regionData<- get500baseWindow (regionData)
#Add DNA sequences
regionData <- addDNAsequences(regionData)

#Save dataframe
write_xlsx (regionData, file.path(pathToOutputsDir,'regionDataCNN.xlsx'))
#Read df from excel
#regionData = read_xlsx(file.path(pathToOutputsDir,'regionDataCNN.xlsx'))
#Subset
#regionData = regionData[1:20000,]

#Get indices for test and train sets
chromosomesTest = c('chr9', 'chr8', 'chr13', 'chr14')
indicesSplit = getTrainTestIndicesFromChr (chromosomesTest, regionData)
trainIndex = indicesSplit$trainIndex
testIndex = indicesSplit$testIndex

#Convert to tensors and split in training, test and validation sets
inputData = getInputCNN (regionData, 
                         testPercentage = 20,
                         trainIndex = trainIndex,
                         testIndex = testIndex)

#Extract train set
xTrain = inputData$train$x
yTrain = inputData$train$y
#Extract test set
xTest = inputData$test$x
yTest = inputData$test$y

#Create model
inputShape <- c(dim(xTrain)[2], dim(xTrain)[3])
nClasses <- dim(yTrain)[2]
cnnModel <- createModel(inputShape, nClasses)

#Train model
cnnModel = trainModel(xTrain, yTrain,
                      cnnModel,
                      batchSize = 128, 
                      epochs = 30, 
                      patience = 10,
                      valSplit = 0.2,
                      pathToPlotsDir = pathToPlotsDir)

#Load model
#cnnModel <- load_model__hdf5(file.path(pathToOutputsDir,"cnnModel.hdf5"))

#Evaluate model on test set
cnnModel %>% evaluate(xTest, yTest)

#Get predictions
yPred <- getPredictions (cnnModel, xTest, pathToOutputsDir)

#Save model
save_model_hdf5(cnnModel, file.path(pathToOutputsDir,"cnnModel.hdf5"))
print (summary(cnnModel))

#Plot the predictions vs true values
yPred <- read_xlsx(file.path(pathToOutputsDir,'TopicPredictions.xlsx'))
yTest<- regionData [testIndex, ]

plotPredictedvsTrue (yPred, yTest, pathToPlotsDir)