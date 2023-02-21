library(keras)
library(writexl)
library(readxl)
source('utils/plots.R')
source('utils/CNNcellTypesFunctions.R')
source('utils/data.R')

#Define path to store plots
pathToPlotsDir <- 'plots/cellTypePrediction'
if (!dir.exists(file.path (pathToPlotsDir))){
  dir.create(file.path (pathToPlotsDir))
}

#Define path to store outputs
pathToOutputsDir <- 'outputs/cellTypePrediction'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

###    GET INPUT
#Read input files
ATACseqData <- readRDS('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
metadata <- read.delim('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')

#Create dataframe with read counts per type of cell for each region
cellTypesPerRegion <- getCellTypePerRegion (ATACseqData, metadata)

###   PROCESS INPUT FOR CNN

#Correct the read count by cell number and transform to log
cellTypesPerRegion <- transformReadCountsToLog(cellTypesPerRegion, metadata)
#Subset each region to 500 bases
cellTypesPerRegion<- get500baseWindow (cellTypesPerRegion)
#Add DNA sequences
cellTypesPerRegion <- addDNAsequences(cellTypesPerRegion)

#Save dataframe
write_xlsx (regionData, file.path(pathToOutputsDir,'cellTypesPeRegionCNN.xlsx'))

#Get indices for test and train sets
chromosomesTest = c('chr9', 'chr8', 'chr13', 'chr14')
indicesSplit = getTrainTestIndicesFromChr (chromosomesTest, cellTypesPerRegion)
trainIndex = indicesSplit$trainIndex
testIndex = indicesSplit$testIndex

