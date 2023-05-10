library(matrixStats)
library(readxl)
library(writexl)
source('utils/data.R')

### Subset based on variance between reads 

#Read input data
ATACpeaks <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
#Empty vector to store variance per region
rowSD <- c()
pb = txtProgressBar(min = 0, max = nrow (ATACpeaks)/3456, style = 3, width = 50) 
#Calculate variance per region and fill in vector
for (i in seq(1, nrow(ATACpeaks), 3457)){
  rowSD <- c(rowSD, rowSds(as.matrix(ATACpeaks[i:(i+3456),])))
  setTxtProgressBar(pb, i)
}
close(pb)

#Subset based on variance 
ATACsubset3rdQ <- ATACpeaks[which(rowSD> summary(rowSD)['3rd Qu.']), ]
saveRDS(ATACsubset3rdQ, file= 'data/ATACsubset3rdQ.Rds')

ATACsubset1stQ <- ATACpeaks[which(rowSD> summary(rowSD)['1st Qu.']), ]
saveRDS(ATACsubset1stQ, file= 'data/ATACsubset1stQ.Rds')

ATACsubsetMedian<- ATACpeaks[which(rowSD> summary(rowSD)['Median']), ]
saveRDS(ATACsubsetMedian, file= 'data/ATACsubsetMedian.Rds')



###  Subset based on read count per cell type

#Read input data
cellTypesPerRegion= read_xlsx('outputs/cellTypePrediction/readCountOutput/cellTypesPeRegion.xlsx')
metadata <- read.delim('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')

#Correct the read count by cell number and transform to log
cellTypesPerRegion <- transformReadCountsToLog(cellTypesPerRegion, metadata)

#Empty vector to store variance per region
rowSD <- rowSds(as.matrix(cellTypesPerRegion[,5:7]))


cellTypesPerRegion1stQ <- cellTypesPerRegion[which(rowSD> summary(rowSD)['1st Qu.']), ]
write_xlsx (cellTypesPerRegion1stQ, 'outputs/cellTypePrediction/readCountOutput/cellTypesPeRegion1stQ.xlsx')
