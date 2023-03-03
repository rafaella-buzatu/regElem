library(Seurat)
library(Signac)
library(writexl)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

ATACpeaks <- readRDS('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
RNAcounts <- readRDS('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_RNA_Counts.RDS')
metadata <- read.delim('ATAC-seq_data/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#Make cell names the row names
metadata <- data.frame(metadata, row.names = 1)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotation) <- 'UCSC'
chromAssay <- CreateChromatinAssay( counts = ATACpeaks,
                                    sep = c(":", "-"),
                                    genome = 'hg38',
                                    annotation = annotation)

SeuratObject <- CreateSeuratObject(counts = chromAssay, assay = "ATACseq", meta.data = metadata)

SeuratObject[["RNA"]] <- CreateAssayObject(counts = RNAcounts)


SeuratObject[['ATACseq']] <- RegionStats( object = SeuratObject[['ATACseq']],
                             genome = BSgenome.Hsapiens.UCSC.hg38)

LinkedPeaks <- LinkPeaks(
  SeuratObject,
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
  verbose = TRUE
   )

dataFrameLinkedPeaks <- as.data.frame(Links(LinkedPeaks)@elementMetadata@listData)
write_xlsx (dataFrameLinkedPeaks, 'outputs/linkedPeaks.xlsx')