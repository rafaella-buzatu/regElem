library(keras)
library(writexl)
library(readxl)
source('utils/plots.R')
source('utils/CNNfunctions.R')
source('utils/data.R')

#Define path to store plots
pathToPlotsDir = 'plots'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir = 'outputs'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

#Read input data
regionData = read_xlsx(file.path(pathToOutputsDir,'regionData.xlsx'))

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

#Convert to tensors and split in training, test and validation sets
inputData = getInputCNN (regionData, testPercentage = 20)

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



