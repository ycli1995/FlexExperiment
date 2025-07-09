#' @importFrom SingleCellExperiment LinearEmbeddingMatrix SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
NULL

#' A small example version of the PBMC dataset
"pbmc_examples"

#' A small example version of the PBMC ATAC dataset
"atac_examples"

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
    return(.load_example_RNA(
      get(data("pbmc_examples", envir = environment()))[[sample]],
      to = "FE"
    ))
  }
  .load_example_ATAC(
    get(data("atac_examples", envir = environment()))[[sample]],
    to = "FE"
  )
}

#' @rdname load-examples
#' @export
exampleSCFE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    return(.load_example_RNA(
      get(data("pbmc_examples", envir = environment()))[[sample]],
      to = "SCFE"
    ))
  }
  .load_example_ATAC(
    get(data("atac_examples", envir = environment()))[[sample]],
    to = "SCFE"
  )
}

#' @rdname load-examples
#' @export
exampleSE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    return(.load_example_RNA(
      get(data("pbmc_examples", envir = environment()))[[sample]],
      to = "SE"
    ))
  }
  .load_example_ATAC(
    get(data("atac_examples", envir = environment()))[[sample]],
    to = "SE"
  )
}

#' @rdname load-examples
#' @export
exampleSCE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    return(.load_example_RNA(
      get(data("pbmc_examples", envir = environment()))[[sample]],
      to = "SCE"
    ))
  }
  .load_example_ATAC(
    get(data("atac_examples", envir = environment()))[[sample]],
    to = "SCE"
  )
}


load_example <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    return(.load_example_RNA(
      get("pbmc_examples", envir = environment()),
      to = "SCFE"
    ))
  }
  .load_example_ATAC(
    get("atac_examples", envir = environment()),
    to = "SCFE"
  )
}

#' @rdname load-examples
#' @export
load_example_SCE <- function(sample = NULL, type = c("RNA", "ATAC")) {
  type <- match.arg(type)
  sample <- sample %||% .all_examples[1]
  sample <- match.arg(sample, .all_examples)
  if (type == "RNA") {
    return(.load_example_RNA(
      get("pbmc_examples", envir = environment()),
      to = "SCE"
    ))
  }
  .load_example_ATAC(
    get("atac_examples", envir = environment()),
    to = "SCE"
  )
}

.all_examples <- c(
  "pbmc3k_unsorted",
  "pbmc3k_sorted",
  "pbmc10k_unsorted",
  "pbmc10k_sorted"
)

.load_LEM <- function(input) {
  lem <- LinearEmbeddingMatrix(
    sampleFactors = unname(input$embeddings),
    featureLoadings = unname(input$loadings),
    factorData = DataFrame(
      stdev = input$stdev,
      row.names = colnames(input$embeddings)
    )
  )
  rownames(lem) <- rownames(input$embeddings)
  lem
}

.load_example_RNA <- function(
    example_data,
    to = c("SCFE", "FE", "SCE", "SE")
) {
  to <- match.arg(to)

  pca <- .load_LEM(example_data$pca)
  tsne <- example_data$tsne
  umap <- example_data$umap
  assays <- list(
    counts = example_data$counts,
    logcounts = example_data$logcounts,
    scale.data = example_data$scale.data
  )
  if (to == "FE") {
    return(FlexExperiment(
      assays = assays,
      rowData = example_data$rowData,
      colData = example_data$colData
    ))
  }
  if (to == "SE") {
    return(SummarizedExperiment(
      assays = assays[c("counts", "logcounts")],
      rowData = DataFrame(example_data$rowData),
      colData = DataFrame(example_data$colData)
    ))
  }

  knn <- nnToPairedHits(example_data$knn$idx, example_data$knn$dist)
  snn <- matToPairedHits(example_data$snn)
  if (to == "SCFE") {
    return(SingleCellFlexExperiment(
      assays = assays,
      rowData = example_data$rowData,
      colData = example_data$colData,
      reducedDims = list(PCA = pca, TSNE = tsne, UMAP = umap),
      colPairs = list(KNN = knn, SNN = example_data$snn)
    ))
  }
  knn <- as(knn, "SelfHits")
  snn <- as(snn, "SelfHits")
  SingleCellExperiment(
    assays = assays[c("counts", "logcounts")],
    rowData = DataFrame(example_data$rowData),
    colData = DataFrame(example_data$colData),
    reducedDims = list(PCA = pca, TSNE = tsne, UMAP = umap),
    colPairs = list(kNN = knn, SNN = snn)
  )
}

.load_example_ATAC <- function(
    example_data,
    to = c("SCFE", "FE", "SCE", "SE")
) {
  to <- match.arg(to)

  lsi <- .load_LEM(example_data$lsi)
  tsne <- example_data$tsne
  umap <- example_data$umap
  assays <- list(
    counts = example_data$counts,
    logcounts = example_data$logcounts,
    scale.data = example_data$scale.data
  )
  if (to == "FE") {
    return(FlexExperiment(
      assays = assays,
      rowRanges = example_data$rowRanges,
      rowData = example_data$rowData,
      colData = example_data$colData
    ))
  }
  if (to == "SE") {
    obj <- SummarizedExperiment(
      assays = assays[c("counts", "logcounts")],
      rowRanges = example_data$rowRanges,
      colData = example_data$colData
    )
    rowData(obj) <- DataFrame(example_data$rowData)
    return(obj)
  }

  knn <- nnToPairedHits(example_data$knn$idx, example_data$knn$dist)
  snn <- matToPairedHits(example_data$snn)
  if (to == "SCFE") {
    return(SingleCellFlexExperiment(
      assays = assays,
      rowRanges = example_data$rowRanges,
      rowData = example_data$rowData,
      colData = example_data$colData,
      reducedDims = list(LSI = lsi, TSNE = tsne, UMAP = umap),
      colPairs = list(KNN = knn, SNN = example_data$snn)
    ))
  }
  knn <- as(knn, "SelfHits")
  snn <- as(snn, "SelfHits")
  obj <- SingleCellExperiment(
    assays = assays[c("counts", "logcounts")],
    rowRanges = example_data$rowRanges,
    colData = DataFrame(example_data$colData),
    reducedDims = list(LSI = lsi, TSNE = tsne, UMAP = umap),
    colPairs = list(kNN = knn, SNN = snn)
  )
  rowData(obj) <- DataFrame(example_data$rowData)
  obj
}
