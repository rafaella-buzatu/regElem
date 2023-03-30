library(keras)
library(writexl)
library(readxl)
source('utils/plots.R')
source('utils/CNNcellTypesFunctions.R')
source('utils/data.R')

#Define path to store plots
pathToPlotsDir <- 'plots/cellTypePrediction/subsetInputMed_readCountOutput'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir <- 'outputs/cellTypePrediction/subsetInputMed_readCountOutput'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

###    GET INPUT
#Read input file
#ATACseqData <- readRDS('data/ATACsubsetMedian.Rds')
ATACseqData <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
metadata <- read.delim('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')

#Create dataframe with read counts per type of cell for each region
cellTypesPerRegion <- getCellTypePerRegion (ATACseqData, metadata, countType = 'reads')
write_xlsx (cellTypesPerRegion, file.path(pathToOutputsDir,'cellTypesPeRegion.xlsx'))
#cellTypesPerRegion= read_xlsx('outputs/cellTypePrediction/readCountOutput/cellTypesPeRegion.xlsx')


###   PROCESS INPUT FOR CNN

#Correct the read count by cell number and transform to log
cellTypesPerRegion <- transformReadCountsToLog(cellTypesPerRegion, metadata)
#Transform read counts to percentages
#cellTypesPerRegion <- transformReadCountsToPercentages(cellTypesPerRegion, metadata)
#Subset each region to 500 bases
cellTypesPerRegion<- getBaseWindow (cellTypesPerRegion, windowSize = 500)
#Add DNA sequences
cellTypesPerRegion <- addDNAsequences(cellTypesPerRegion)
#Binarize
#cellTypesPerRegion <- binarizeCellTypeScore (cellTypesPerRegion, thresholdNumeric= 5)

#Save dataframe
write_xlsx (cellTypesPerRegion, file.path(pathToOutputsDir,'cellTypesPeRegionCNN.xlsx'))
#cellTypesPerRegion= read_xlsx(file.path(pathToOutputsDir,'cellTypesPeRegionCNN.xlsx'))

#Get indices for test and train sets
chromosomesTest = c('chr9', 'chr8', 'chr13', 'chr14')
indicesSplit = getTrainTestIndicesFromChr (chromosomesTest, cellTypesPerRegion)
trainIndex = indicesSplit$trainIndex
testIndex = indicesSplit$testIndex

#Convert to tensors and split in training, test and validation sets
inputData = getInputCNN (cellTypesPerRegion, 
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
cnnModel <- createModelRegression (inputShape, nClasses)
#cnnModel <- createModelMultiLabelClassification (inputShape, nClasses)

#Train model
cnnModel = trainModel(xTrain, yTrain,
                      cnnModel,
                      batchSize = 128, 
                      epochs = 100, 
                      patience = 20,
                      valSplit = 0.2,
                      pathToPlotsDir = pathToPlotsDir)

#Load model
#cnnModel <- load_model_hdf5(file.path(pathToOutputsDir,"cnnModel.hdf5"))

#Evaluate model on test set
cnnModel %>% evaluate(xTest, yTest)

#Get predictions
listCellTypes <- colnames(cellTypesPerRegion[, 5:(ncol(cellTypesPerRegion)-1)])
yPred <- getPredictions (cnnModel, xTest, listCellTypes,pathToOutputsDir, binary = FALSE)

#Save model
save_model_hdf5(cnnModel, file.path(pathToOutputsDir,"cnnModel.hdf5"))
print (summary(cnnModel))

#Plot the predictions vs true values
yPred <- read_xlsx(file.path(pathToOutputsDir,'CellTypePredictions.xlsx'))
yTest<- cellTypesPerRegion [testIndex, ]

#For Regression
plotPredictedvsTrue (yPred, yTest,  pathToPlotsDir)

