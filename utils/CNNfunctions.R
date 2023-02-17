### DL functions
library (keras)
library(tensorA)
library(tensorflow)
library(ptw)
library(writexl)

source('utils/data.R')


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


splitTrainTest <- function (listOneHotMatrices, listLabels, testPercentage){
  #' Split the input list of matrices into training and test sets
  
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


getInputCNN <- function(regionData, testPercentage){
  #' Formats the input dataframe into input data for the CNN
  
  #Define empty lists to save the inputs
  listOneHotMatrices <- vector(mode = "list", length = (nrow(regionData)))
  listLabels <- vector (mode = "list", length = nrow(regionData))
  
  #Get names of columns containing topic probabilities
  topicColNames = colnames(regionData[ , grepl( "Topic" , names( regionData))])
  
  print ('Encoding DNA sequences...')
  #Initialize progress bar
  pb = txtProgressBar(min = 0, max = nrow(regionData), style = 3, width = 50) 
  #Extract inputs and outputs
  for (i in (1:nrow(regionData))){
    #Convert DNA sequences of regions to one hot encoded matrices
    encodedMatrix = oneHotEncode (unlist(regionData [i, 'DNAseq']))
    #Add matrix to list
    listOneHotMatrices [[i]] = encodedMatrix
    #Extract values of topic probabilities
    listLabels[[i]] = unlist(regionData[i, topicColNames])
    
    #Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  print ('Splitting into training and test set...')
  #Split into training, test sets
  inputDataset = splitTrainTest(listOneHotMatrices, listLabels, testPercentage)
  
  return (inputDataset)
}
  
createModel <- function (inputShape, nClasses){
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
    layer_dense(500,activation="relu")%>%
    layer_dropout(0.25)%>%
    layer_dense(units = nClasses, activation = 'relu')
  
  cnnModel %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adadelta(),
    metrics = c("mean_squared_error")
  )
  
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
  plot(cnnHistory )
  dev.off()
  
  #Save model
  save_model_tf(cnnModel, file.path(pathToOutputsDir,"cnnModel.hdf5"))
  
  return (cnnModel)
}

getPredictions <- function (cnnModel, xTest, pathToOutputDir){
  #Extract predictions
  yPred <- as.data.frame(predict(cnnModel, xTest))
  
  #Get no of topics
  topicNo = dim(yPred)[2]
  #Rename columns
  columns = c()
  for (i in 1:topicNo){
    columns = append (columns,  paste(c('Topic', i), collapse = ""))
  }
  colnames (yPred) <- columns
  
  write_xlsx (yPred, file.path(pathToOutputsDir,'TopicPredictions.xlsx'))
  
  return (yPred)
}