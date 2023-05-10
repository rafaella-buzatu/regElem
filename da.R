library(Seurat)
library(Signac)
library(writexl)
library(readxl)
library (EnsDb.Hsapiens.v86)

#Read input data
ATACpeaks <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS')
#RNAcounts <- readRDS('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_RNA_Counts.RDS')
metadata <- read.delim('data/SNAREseqData/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
#Make cell names the row names
metadata <- data.frame(metadata, row.names = 1)

metadata['NonNeuronal'] = 0
metadata[metadata['level1'] == 'Non-Neuronal', 'NonNeuronal'] = 'Yes'
metadata[metadata['level1'] != 'Non-Neuronal', 'NonNeuronal'] = 'No'

metadata['Glutamatergic'] = 0
metadata[metadata['level1'] == 'Glutamatergic', 'Glutamatergic'] = 'Yes'
metadata[metadata['level1'] != 'Glutamatergic', 'Glutamatergic'] = 'No'

metadata['GABAergic'] = 0
metadata[metadata['level1'] == 'GABAergic', 'GABAergic'] = 'Yes'
metadata[metadata['level1'] != 'GABAergic', 'GABAergic'] = 'No'

#Add ATACpeaks as chromatin assay 
chromAssay <- CreateChromatinAssay( counts = ATACpeaks,
                                    sep = c(":", "-"),
                                    genome = 'hg38')
SeuratObject <- CreateSeuratObject(counts = chromAssay, assay = "ATACseq", meta.data = metadata)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(SeuratObject) <- annotations

DefaultAssay(SeuratObject) <- 'ATACseq'


SeuratObject <- SetIdent(SeuratObject, value = 'NonNeuronal')

da_peaks_NN <- FindMarkers(
  object = SeuratObject,
  ident.1 = 'Yes',
  ident.2 = 'No',
  test.use = 'LR',
  min.pct = 0.05,
  logfc.threshold = 0
)

SeuratObject <- SetIdent(SeuratObject, value = "Glutamatergic")

da_peaks_GL <- FindMarkers(
  object = SeuratObject,
  ident.1 = 'Yes',
  ident.2 = 'No',
  test.use = 'LR',
  min.pct = 0.05,
  logfc.threshold = 0
)


SeuratObject <- SetIdent(SeuratObject, value = "GABAergic")

da_peaks_GB <- FindMarkers(
  object = SeuratObject,
  ident.1 = 'Yes',
  ident.2 = 'No',
  test.use = 'LR',
  min.pct = 0.05,
  logfc.threshold = 0
)

write_xlsx (da_peaks_NN, file.path('outputs/da/da_peaks_NN.xlsx'))

