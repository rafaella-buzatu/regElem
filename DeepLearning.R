library(keras)
library(writexl)
library(readxl)
source('utils/plots.R')
source('utils/CNNfunctions.R')

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
regionData = regionData[1:2000,]
#Convert to tensors and split in training, test and validation sets
inputData = getInputCNN (regionData, testPercentage = 20)

#Extract train set
xTrain = inputData$train$x
yTrain = inputData$train$y
#Extract test set
xTest = inputData$test$x
yTest = inputData$test$y

#Create model
inputShape <- c(dim(xTrain)[2], dim(xTrain)[3], 1)
nClasses <- dim(yTrain)[2]
cnnModel <- createModel(inputShape, nClasses)

#Train model
batchSize <- 128
epochs <- 50

cnnHistory <- cnnModel %>% fit(
  xTrain, yTrain,
  batch_size = batchSize,
  epochs = epochs,
  validation_split = 0.2
)

#Plot evolution of loss function and performance metrics
plotCNNhistory(cnnHistory, pathToPlotsDir)

#Evaluate model on test set
cnnModel %>% evaluate(xTest, yTest)
#cnnModel$acc

