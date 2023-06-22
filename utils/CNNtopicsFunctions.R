### DL functions
library (keras)
library(tensorA)
library(tensorflow)
library(ptw)
library(writexl)
library(caret)
library (ramify)

source('utils/data.R')

getCellTypeFromTopics<- function ( metadata, cellData ){
  
  
  topics= colnames(cellData [ , grepl( "Topic" , names( cellData  ) )])
  colsToKeep = append(topics, 'sampleName',after=0)
  
  cellDataTopics = cellData[c(colsToKeep)]
  cellDataTopics['cellType'] <- NA
  
  for (i in 1:nrow(cellDataTopics)){
    sample = unlist(cellDataTopics[i, 'sampleName'])
    cellType = metadata[metadata['sample_name'] == sample, 'level1']
    cellDataTopics[i, 'cellType'] = cellType
  }
  
  return (cellDataTopics)
  
  
}


splitTrainTest <- function (cellTypesFromTopics, testPercentage, seed = 10){
  #' Split the input dataframe into training and test sets
  
  set.seed(seed)
  
  sample <- sample(c(TRUE, FALSE), nrow(cellTypesFromTopics), replace=TRUE, prob=c(1-testPercentage,testPercentage))
  train  <- cellTypesFromTopics[sample, ]
  test   <- cellTypesFromTopics[!sample, ]
  
  #EXTRACT ACTUAL VALUES
  
  xTrain = train[, !(names(train) %in% c('cellType'))]
  yTrain = train[, 'cellType']
  labels = unique(yTrain)[[1]]
  #One-hot-encode yTrain
  yTrain = oneHotEncodeLabels(yTrain, labels)
  #Save X and Y in list
  trainSet = list (x = as.tensor(data.matrix(xTrain)), y = as.tensor(data.matrix(yTrain)), labels = labels)
  
  xTest = test[,  !(names(train) %in% c('cellType'))]
  yTest= test[, 'cellType']
  #One-hot-encode yTest
  yTest = oneHotEncodeLabels(yTest, labels)
  #Save X and Y in list
  testSet = list(x = as.tensor(data.matrix(xTest)), y = as.tensor(data.matrix(yTest)), labels = labels)
  
  #Save test and training set in a list
  inputData = list (train = trainSet, test = testSet)
  
  return (inputData)
  
}

oneHotEncodeLabels <- function (dataset, labels){
  
  #Create matrix of zeros of size sequence length 
  oneHotMatrix <- matrix(data = rep(0), ncol = length(labels), nrow = nrow(dataset))
  
  for (i in 1:nrow(dataset)){
    oneHotMatrix[i, which (labels == dataset[i, ][[1]])] = 1
  }
  
  return (oneHotMatrix)
}

createModelMultiLabelClassification <- function (inputShape, nClasses){
  #' Creates and compiles the convolutional model architecture
  
  cnnModel <- keras_model_sequential() %>%
    layer_dense(7,activation="relu",input_shape = c(inputShape))%>%
    layer_dropout(0.25)%>%
    layer_dense(units = nClasses, activation = 'softmax')
  
  cnnModel %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = 'adam',
    metrics = 'accuracy')
  
  
  return (cnnModel)
}

trainModel <-function (xTrain, yTrain, cnnModel, batchSize, epochs, patience,
                       valSplit, pathToPlotsDir){
  
  #Initialize earlyStopping
  callbacks <- list(callback_early_stopping(monitor = "val_loss", patience = patience,
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

getPredictions <- function (nnModel, xTest){
  #Extract predictions
  yPred <-data.frame(argmax(predict(nnModel, xTest), rows = TRUE))
  colnames(yPred) = c('cellType')
  
  return (yPred)
}

plotConfusionMatrix <-function (yPred, yTrue, labels, pathToPlotsDir){
  cm <- confusionMatrix(factor(unlist(yPred)), factor(unlist(yTrue)), dnn = c("Prediction", "True"))
  pdf(file = file.path(pathToPlotsDir, 'confusionMatrix.pdf'),
      width = 7,
      height =6)
  plt <- as.data.frame(cm$table)
  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  
  p <- ggplot(plt, aes(Prediction,True, fill= Freq)) +
    geom_tile() + geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="#009194") +
    labs(x = "True",y = "Prediction") +
    scale_x_discrete(labels= labels) +
    scale_y_discrete(labels= rev(labels))
  
  print (p)
  
  dev.off()
}