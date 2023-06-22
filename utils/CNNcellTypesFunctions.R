### DL functions
library (keras)
library(tensorA)
library(tensorflow)
library(ptw)
library(writexl)

source('utils/data.R')


transformReadCountsToLog <- function (cellTypesPerRegion, metadata){
  
  ### Get number of cells per type
  #Get names of cell types
  cellTypes =  unlist(unique( metadata['level1']))
  
  
  #Create empty list to store counts
  numberCellsPerType = vector(mode = "list", length = (length(cellTypes)))
  #Add counts to list
  for (type in 1: length(cellTypes)){
    numberCellsPerType [[type]] = nrow(subset(metadata, level1 == cellTypes[[type]]))
  }
  names(numberCellsPerType) = cellTypes
  
  #Get mean of cell numbers
  meanCells <- mean(unlist(numberCellsPerType))
  
  numberCellsPerTypeWeight = vector(mode = "list", length = (length(cellTypes)))
  #Get weight of each cell type
  for (i in 1: length(numberCellsPerType)){
    numberCellsPerTypeWeight [[i]] = numberCellsPerType[[i]]/meanCells
  }
  names(numberCellsPerTypeWeight) = cellTypes
  
  #Initialize progress bar
  pb = txtProgressBar(min = 0, max = nrow(cellTypesPerRegion), style = 3, width = 50) 
  #Iterate over every row in the input dataframe
  for (row in 1:nrow (cellTypesPerRegion)) {
    #Normalize coutns using weight
    for (type in cellTypes){
      cellTypesPerRegion[row, type] = log(cellTypesPerRegion[row, type]/unlist(numberCellsPerTypeWeight[type]) +1)
    }
    #Update progress bar
    setTxtProgressBar(pb, row) 
  }
  close (pb)
  
  return (cellTypesPerRegion)
}


transformReadCountsToPercentages <- function (cellTypesPerRegion, metadata){
  
  #Get names of cell types
  cellTypes =  unlist(unique( metadata['level1']))
  
  #Initialize progress bar
  pb = txtProgressBar(min = 0, max = nrow(cellTypesPerRegion), style = 3, width = 50) 
  #Extract inputs and outputs
  #Iterate over every row in the input dataframe
  for (row in 1:nrow (cellTypesPerRegion)) {
    totalReads = sum (cellTypesPerRegion [row, cellTypes])
    for (type in cellTypes){
      cellTypesPerRegion[row, type] = cellTypesPerRegion[row, type] /totalReads 
    }
    #Update progress bar
    setTxtProgressBar(pb, row) 
  }
  close (pb)
  
  return (cellTypesPerRegion)
}



binarizeCellTypeScore <- function (cellTypesPerRegion, thresholdMetric = NULL, thresholdNumeric = NULL) {
  
  for (i in 5:7) {
    if (!is.null(thresholdMetric)){
      threshold = summary (unlist(cellTypesPerRegion [, i]))[thresholdMetric]
    }
    else if (!is.null(thresholdNumeric)) {
      threshold = thresholdNumeric
    }
    cellTypesPerRegion [which (cellTypesPerRegion [, i] <= threshold), i] = 0
    cellTypesPerRegion [which (cellTypesPerRegion [, i] >   threshold), i] = 1
  }
  
  return (cellTypesPerRegion)
}

oneHotEncode <- function(dnaSeq) {
  #'One hot encodes the input dna sequence.
  
  # Construct matrix
  #Define nucleotides
  nucleotides = c("A", "C", "G", "T")
  #Get length of dnaSequence
  seqLength <- nchar(dnaSeq)
  #Split DNA sequence
  seqSplit <- unlist(strsplit(x = dnaSeq, split = ""))
  
  #Create matrix of zeros of size sequence length 
  oneHotMatrix <- matrix(data = rep(0), ncol = seqLength, nrow = 4)
  #rename columns and rows with nucleotides
  rownames(oneHotMatrix) <- nucleotides
  colnames(oneHotMatrix) <- seqSplit
  
  # Encode
  for (i in nucleotides) {
    oneHotMatrix[rownames(oneHotMatrix) == i, colnames(oneHotMatrix) == i] <- 1
  }
  
  oneHotMatrix <- t(oneHotMatrix)
  row.names(oneHotMatrix) <- c(1:dim(oneHotMatrix)[1])
  
  return(oneHotMatrix)
  
}

getTrainTestIndicesFromChr <-function (chromosomesTest, cellTypesPerRegion){
  
  testIndex = which (cellTypesPerRegion$seqnames %in% chromosomesTest)
  trainIndex = which (!(cellTypesPerRegion$seqnames %in% chromosomesTest))
  
  
  indicesSplit = list (trainIndex = trainIndex, testIndex = testIndex)
  
  return (indicesSplit)
}

splitTrainTest <- function (listOneHotMatrices, listLabels, testPercentage,
                            trainIndex = NULL, testIndex = NULL){
  #' Split the input list of matrices into training and test sets
  
  if (is.null(trainIndex) & is.null(testIndex)){
  ### SPLIT INDICES
  #Get indices of input dataset
  allIndex = c((1:length(listOneHotMatrices)))
  
  #SPLIT TEST-TRAIN
  #Get number of records of test set
  noTest = as.integer(testPercentage/100*length(listOneHotMatrices))
  #Randomly sample indices of test set
  set.seed(15)
  testIndex = sample(allIndex, noTest )
  #Get indices of training + val set
  trainIndex = allIndex[!(allIndex %in% testIndex)]
  }
 
  #EXTRACT ACTUAL VALUES
  
  #Extract training set X in a matrix
  xTrain <-  array(numeric(),c(length(listOneHotMatrices[trainIndex]), dim(listOneHotMatrices[[1]])[1], dim(listOneHotMatrices[[1]])[2]))
  listTraining <- listOneHotMatrices[trainIndex]
  for (i in 1: length(listTraining)){
    xTrain [i,  , ] = listTraining [[i]]
  }
  #Extract training set labels in a matrix   
  yTrain <- array(numeric(), c(length(listLabels [trainIndex]), length(listLabels[[1]])))
  listTrainingLabels <- listLabels[trainIndex]
  for (i in 1: length(listTrainingLabels)){
    yTrain [i, ] = listTrainingLabels [[i]]
  }
  #Save X and Y in list
  trainSet = list (x = as.tensor(xTrain), y = as.tensor(yTrain))
  
  
  #Extract test set X in a matrix
  xTest <-  array(numeric(),c(length(listOneHotMatrices[testIndex]), dim(listOneHotMatrices[[1]])[1], dim(listOneHotMatrices[[1]])[2]))
  listTest <- listOneHotMatrices[testIndex]
  for (i in 1: length(listTest)){
    xTest [i,  , ] = listTest [[i]]
  }
  #Extract training set labels in a matrix   
  yTest <- array(numeric(), c(length(listLabels [testIndex]), length(listLabels[[1]])))
  listTestLabels <- listLabels[testIndex]
  for (i in 1: length(listTestLabels)){
    yTest [i, ] = listTestLabels [[i]]
  }
  #Save X and Y in list
  testSet = list (x = as.tensor(xTest), y = as.tensor (yTest))
  
  #Save test and training set in a list
  inputData = list (train = trainSet, test = testSet)
  
  return (inputData)
  
}


getInputCNN <- function(cellTypesPerRegion, testPercentage, trainIndex = NULL,
                        testIndex = NULL){
  #' Formats the input dataframe into input data for the CNN
  
  #Define empty lists to save the inputs
  listOneHotMatrices <- vector(mode = "list", length = (nrow(cellTypesPerRegion)))
  listOutputs <- vector (mode = "list", length = nrow(cellTypesPerRegion))
  
  #Get names of columns containing cell type counts 
  cellTypeColNames = colnames(cellTypesPerRegion[, 5:(ncol(cellTypesPerRegion)-1)])
  
  print ('Encoding DNA sequences...')
  #Initialize progress bar
  pb = txtProgressBar(min = 0, max = nrow(cellTypesPerRegion), style = 3, width = 50) 
  #Extract inputs and outputs
  for (i in (1:nrow(cellTypesPerRegion))){
    #Convert DNA sequences of regions to one hot encoded matrices
    encodedMatrix = oneHotEncode (unlist(cellTypesPerRegion [i, 'DNAseq']))
    #Add matrix to list
    listOneHotMatrices [[i]] = encodedMatrix
    #Extract normalized reads per type
    listOutputs[[i]] = unlist(cellTypesPerRegion[i,  cellTypeColNames])
    
    #Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  print ('Splitting into training and test set...')
  #Split into training, test sets
  inputDataset = splitTrainTest(listOneHotMatrices, listOutputs, testPercentage,
                                trainIndex, testIndex)
  
  return (inputDataset)
}
  
createModelRegression <- function (inputShape, nClasses){
  #' Creates and compiles the convolutional model architecture
  
  cnnModel <- keras_model_sequential() %>%
    layer_conv_1d(filters=64, kernel_size=12, input_shape = inputShape, activation="relu") %>% 
    layer_conv_1d(filters=64, kernel_size=12, activation="relu") %>% 
    layer_max_pooling_1d(pool_size=4) %>%
    layer_dropout(0.25) %>%
    layer_conv_1d(filters=32, kernel_size=2, activation="relu") %>%
    layer_conv_1d(filters=32, kernel_size=12, activation="relu") %>%
    layer_max_pooling_1d(pool_size=4) %>%
    layer_dropout(0.25) %>%
    layer_flatten() %>%
    layer_dense(units = nClasses, activation = 'relu')
  
  cnnModel %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adadelta()
  )
  
  return (cnnModel)
}


createModelMultiLabelClassification <- function (inputShape, nClasses){
  #' Creates and compiles the convolutional model architecture
  
  cnnModel <- keras_model_sequential() %>%
    layer_conv_1d(filters=64, kernel_size=12, input_shape = inputShape, activation="relu") %>% 
    layer_conv_1d(filters=64, kernel_size=12, activation="relu") %>% 
    layer_max_pooling_1d(pool_size=4) %>%
    layer_dropout(0.25) %>%
    layer_conv_1d(filters=32, kernel_size=2, activation="relu") %>%
    layer_conv_1d(filters=32, kernel_size=12, activation="relu") %>%
    layer_max_pooling_1d(pool_size=4) %>%
    layer_dropout(0.25) %>%
    layer_flatten() %>%
    layer_dense(200,activation="relu")%>%
    layer_dropout(0.25)%>%
    layer_dense(50,activation="relu")%>%
    layer_dropout(0.25)%>%
    layer_dense(units = nClasses, activation = 'sigmoid')
  
  cnnModel %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_adadelta())

  
  return (cnnModel)
}

trainModel <-function (xTrain, yTrain, cnnModel, batchSize, epochs, patience,
                       valSplit, pathToPlotsDir){
  
  #Initialize earlyStopping
  callbacks <- list(callback_early_stopping(monitor = "val_loss", patience = 20,
                                            restore_best_weights = TRUE, mode = "auto"))
  
  #Train model on training set
  cnnHistory <- cnnModel %>% fit(
    xTrain, yTrain,
    batch_size = batchSize,
    epochs = epochs,
    callbacks = callbacks,
    validation_split = valSplit
  )
  
  #Plot evolution of loss function and performance metrics
  pdf(file = file.path(pathToPlotsDir, 'lossMetric.pdf'),
      width = 7,
      height =10)
  p<- plot(cnnHistory )
  print (p)
  dev.off()
  
  #Save model
  save_model_tf(cnnModel, file.path(pathToOutputsDir,"cnnModel.hdf5"))
  
  return (cnnModel)
}

getPredictions <- function (cnnModel, xTest, listCellTypes, pathToOutputDir, binary = FALSE){
  #Extract predictions
  yPred <- as.data.frame(predict(cnnModel, xTest))
  
  if (binary == TRUE) {
    for (i in 1:ncol(yPred)){
      yPred [which (yPred[, i] > 0.5), i] = 1
      yPred [which (yPred[, i] < 0.5), i] = 0
    }
  }
  
  colnames (yPred) <- listCellTypes
  
  write_xlsx (yPred, file.path(pathToOutputsDir,'CellTypePredictions.xlsx'))
  
  return (yPred)
}

