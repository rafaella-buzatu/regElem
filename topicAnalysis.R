library(Seurat)
library(Signac)
library(writexl)
library(readxl)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
source('utils/GSEA.R')

###     LINK PEAKS TO GENES

#Read input data
ATACpeaks <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
RNAcounts <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_RNA_Counts.RDS')
metadata <- read.delim('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#Make cell names the row names
metadata <- data.frame(metadata, row.names = 1)

#Create Seurat Object 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotation) <- 'UCSC'

#Add ATACpeaks as chromatin assay 
chromAssay <- CreateChromatinAssay( counts = ATACpeaks,
                                    sep = c(":", "-"),
                                    genome = 'hg38',
                                    annotation = annotation)
SeuratObject <- CreateSeuratObject(counts = chromAssay, assay = "ATACseq", meta.data = metadata)
#Add RNA assay
SeuratObject[["RNA"]] <- CreateAssayObject(counts = RNAcounts)

#Link Peaks to genes
SeuratObject[['ATACseq']] <- RegionStats( object = SeuratObject[['ATACseq']],
                             genome = BSgenome.Hsapiens.UCSC.hg38)

LinkedPeaks <- LinkPeaks(SeuratObject,
                         "ATACseq",
                         "RNA",
                         peak.slot = 'counts',
                         expression.slot = 'data',
                         method = "pearson",
                         gene.coords = NULL,
                         distance = 5e+05,
                         min.distance = NULL,
                         min.cells = 10,
                         genes.use = NULL,
                         n_sample = 200,
                         pvalue_cutoff = 0.05,
                         score_cutoff = 0.05,
                         gene.id = FALSE,
                         verbose = TRUE)
                         
# Extract information into a dataframe
dataFrameLinkedPeaks <- as.data.frame(Links(LinkedPeaks)@elementMetadata@listData)
#Save dataframe
write_xlsx (dataFrameLinkedPeaks, 'outputs/linkedPeaks.xlsx')


###      PATHWAY ANALYSIS (GSEA)

#Define path for outputs
pathToOutputsDir = 'outputs/cisTopic/run3'
if (!dir.exists(file.path (pathToOutputsDir))){
  dir.create(file.path (pathToOutputsDir))
}

#Read region-topic assignments
regionData = read_xlsx(file.path(pathToOutputsDir,'regionData.xlsx'))
#Read linked peaks
linkedPeaks = read_xlsx('outputs/linkedPeaks.xlsx')

#Create vectors of Topic Scores 
topicScoreGenes = getGeneTopicScores(regionData, linkedPeaks)

#Read gene set assignments for each pathway
geneSets = readGMT("data/c5.go.v7.4.symbols.gmt") ### GO terms

#Run GSE analysis on all topics and save results in an excel file
runGSEAallTopics (topicScoreGenes, geneSets, pathToOutputsDir)