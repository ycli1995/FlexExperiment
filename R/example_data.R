#' @importFrom magrittr %>%
#' @importFrom SingleCellExperiment LinearEmbeddingMatrix SingleCellExperiment
#'
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
NULL

#' Load the example datasets
#'
#' Four subsetted versions of 10x Genomics' PBMC multi-omics datasets
#'
#' @param sample One of the following sample names:
#' \itemize{
#' \item "pbmc3k_unsorted"
#' \item "pbmc3k_sorted"
#' \item "pbmc10k_unsorted"
#' \item "pbmc10k_sorted"
#' }
#'
#' @name load-examples
#' @export
exampleFE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    data("pbmc_examples", envir = environment())
    obj <- .load_example_RNA(pbmc_examples[[sample]], to = "FE")
    return(obj)
  }
  data("atac_examples", envir = environment())
  obj <- .load_example_RNA(atac_examples[[sample]], to = "FE")
  return(obj)
}

#' @rdname load-examples
#' @export
exampleSCFE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    data("pbmc_examples", envir = environment())
    obj <- .load_example_RNA(pbmc_examples[[sample]], to = "SCFE")
    return(obj)
  }
  data("atac_examples", envir = environment())
  obj <- .load_example_RNA(atac_examples[[sample]], to = "SCFE")
  return(obj)
}

#' @rdname load-examples
#' @export
exampleSCE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    data("pbmc_examples", envir = environment())
    obj <- .load_example_RNA(pbmc_examples[[sample]], to = "SCE")
    return(obj)
  }
  data("atac_examples", envir = environment())
  obj <- .load_example_RNA(atac_examples[[sample]], to = "SCE")
  return(obj)
}


load_example <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    data("pbmc_examples", envir = environment())
    obj <- .load_example_RNA(pbmc_examples[[sample]], to = "SCFE")
    return(obj)
  }
  data("atac_examples", envir = environment())
  obj <- .load_example_RNA(atac_examples[[sample]], to = "SCFE")
  return(obj)
}

#' @rdname load-examples
#' @export
load_example_SCE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    data("pbmc_examples", envir = environment())
    obj <- .load_example_RNA(pbmc_examples[[sample]], to = "SCE")
    return(obj)
  }
  data("atac_examples", envir = environment())
  obj <- .load_example_RNA(atac_examples[[sample]], to = "SCE")
  return(obj)
}

.all_examples <- c(
  "pbmc3k_unsorted",
  "pbmc3k_sorted",
  "pbmc10k_unsorted",
  "pbmc10k_sorted"
)

.load_LEM <- function(input) {
  lem <- LinearEmbeddingMatrix(
    sampleFactors = input$embeddings,
    featureLoadings = input$loadings,
    factorData = DataFrame(
      stdev = input$stdev,
      row.names = colnames(input$embeddings)
    )
  )
  rownames(lem) <- rownames(input$embeddings)
  lem
}

.load_example_RNA <- function(example_data, to = c("SCFE", "FE", "SCE")) {
  to <- match.arg(to)

  pca <- .load_LEM(example_data$pca)
  tsne <- example_data$tsne
  umap <- example_data$umap
  if (to == "FE") {
    obj <- FlexExperiment(
      assays = list(
        counts = example_data$counts,
        logcounts = example_data$logcounts,
        scale.data = example_data$scale.data
      ),
      rowData = example_data$rowData,
      colData = example_data$colData
    )
    return(obj)
  }

  knn <- nnToPairedHits(example_data$knn$idx, example_data$knn$dist)
  snn <- matToPairedHits(example_data$snn)
  if (to == "SCFE") {
    obj <- SingleCellFlexExperiment(
      assays = list(
        counts = example_data$counts,
        logcounts = example_data$logcounts,
        scale.data = example_data$scale.data
      ),
      rowData = example_data$rowData,
      colData = example_data$colData,
      reducedDims = list(PCA = pca, TSNE = tsne, UMAP = umap),
      colPairs = list(KNN = knn, SNN = example_data$snn)
    )
    return(obj)
  }
  knn <- as(knn, "SelfHits")
  snn <- as(snn, "SelfHits")
  obj <- SingleCellExperiment(
    assays = list(
      counts = example_data$counts,
      logcounts = example_data$logcounts
    ),
    rowData = DataFrame(example_data$rowData),
    colData = DataFrame(example_data$colData),
    reducedDims = list(PCA = pca, TSNE = tsne, UMAP = umap),
    colPairs = list(kNN = knn, SNN = snn)
  )
  obj
}

.load_example_ATAC <- function(example_data, to = c("SCFE", "FE", "SCE")) {
  to <- match.arg(to)

  lsi <- .load_LEM(example_data$lsi)
  tsne <- example$tsne
  umap <- example_data$umap
  if (to == "FE") {
    obj <- FlexExperiment(
      assays = list(
        counts = example_data$counts,
        logcounts = example_data$logcounts,
        scale.data = example_data$scale.data
      ),
      rowRanges = example_data$rowRanges,
      rowData = example_data$rowData,
      colData = example_data$colData
    )
    return(obj)
  }



  knn <- nnToPairedHits(example_data$knn$idx, example_data$knn$dist)
  snn <- matToPairedHits(example_data$snn)
  if (to == "SCFE") {
    obj <- SingleCellFlexExperiment(
      assays = list(
        counts = example_data$counts,
        logcounts = example_data$logcounts,
        scale.data = example_data$scale.data
      ),
      rowRanges = example_data$rowRanges,
      rowData = example_data$rowData,
      colData = example_data$colData,
      reducedDims = list(LSI = lsi, TSNE = tsne, UMAP = umap),
      colPairs = list(KNN = knn, SNN = example_data$snn)
    )
    return(obj)
  }
  knn <- as(knn, "SelfHits")
  snn <- as(snn, "SelfHits")
  obj <- SingleCellExperiment(
    assays = list(
      counts = example_data$counts,
      logcounts = example_data$logcounts
    ),
    rowRanges = example_data$rowRanges,
    colData = DataFrame(example_data$colData),
    reducedDims = list(PCA = pca, TSNE = tsne, UMAP = umap),
    colPairs = list(kNN = knn, SNN = snn)
  )
  rowData(obj) <- DataFrame(example_data$rowData)
  obj
}
