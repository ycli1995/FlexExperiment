
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Roxygen2 calls ###############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.doc_links <- function(nm) {
  pkg <- .pkg_map[nm]
  paste0("\\code{\\link[", pkg, "]{", nm, "}}")
}

.dot_param <- "Arguments passed to other metheds."
.val_param <- "An object of a class specified in the S4 method signature."
.vb_param <- "Print progress."

.pkg_map <- c(

  "SelfHits" = "S4Vectors:Hits-class",

  "drop0" = "Matrix",
  "CsparseMatrix" = "Matrix:CsparseMatrix-class",

  "SVT_SparseMatrix" = "SparseArray",

  "DataFrame" = "S4Vectors",

  "GenomicRanges" = "GenomicRanges:GenomicRanges-class",
  "GRangesList" = "GenomicRanges:GRangesList-class",
  "GRanges" = "GenomicRanges:GRanges-class",

  "SummarizedExperiment" = "RangedSummarizedExperiment-class",
  "assay" = "SummarizedExperiment:SummarizedExperiment-class",
  "assays" = "SummarizedExperiment:SummarizedExperiment-class",
  "colData" = "SummarizedExperiment:SummarizedExperiment-class",
  "rowData" = "SummarizedExperiment:SummarizedExperiment-class",

  "SCE-internals" = "SingleCellExperiment:SCE-internals",
  "int_colData" = "SingleCellExperiment:SCE-internals",
  "int_colData<-" = "SingleCellExperiment:SCE-internals",
  "int_elementMetadata" = "SingleCellExperiment:SCE-internals",
  "int_elementMetadata<-" = "SingleCellExperiment:SCE-internals",
  "reducedDim" = "SingleCellExperiment",
  "reducedDims" = "SingleCellExperiment",
  "reducedDimNames" = "SingleCellExperiment",
  "colPair" = "SingleCellExperiment",
  "colPairs" = "SingleCellExperiment",
  "colPairNames" = "SingleCellExperiment",
  "rowPair" = "SingleCellExperiment",
  "rowPairs" = "SingleCellExperiment",
  "rowPairNames" = "SingleCellExperiment",
  "altExp" = "SingleCellExperiment",
  "altExps" = "SingleCellExperiment",
  "altExpNames" = "SingleCellExperiment",


  "MultiAssayExperiment" = "MultiAssayExperiment",
  "getWithColData" = "MultiAssayExperiment",
  "experiments" = "MultiAssayExperiment:MultiAssayExperiment-methods",
  "subsetBy" = "MultiAssayExperiment",

  "Fragment" = "Signac",
  "Annotation" = "Signac",
  "ChromatinAssay" = "Signac:ChromatinAssay-class",

  "CsparseMatrix" = "Matrix:CsparseMatrix-class",
  "RsparseMatrix" = "Matrix:RsparseMatrix-class",
  "TsparseMatrix" = "Matrix:TsparseMatrix-class"
)
