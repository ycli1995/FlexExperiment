# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

col_merge_dgcmatrix_cpp <- function(mat_list, mat_rownames, all_rownames) {
    .Call('_FlexExperiment_col_merge_dgcmatrix_cpp', PACKAGE = 'FlexExperiment', mat_list, mat_rownames, all_rownames)
}

row_merge_dgcmatrix_cpp <- function(mat_list, mat_colnames, all_colnames) {
    .Call('_FlexExperiment_row_merge_dgcmatrix_cpp', PACKAGE = 'FlexExperiment', mat_list, mat_colnames, all_colnames)
}

rcpp_knn_to_sparse_weights <- function(nn_idx, nn_dist, self_loops) {
    .Call('_FlexExperiment_rcpp_knn_to_sparse_weights', PACKAGE = 'FlexExperiment', nn_idx, nn_dist, self_loops)
}

rcpp_knn_to_sparse_idx <- function(nn_idx, self_loops) {
    .Call('_FlexExperiment_rcpp_knn_to_sparse_idx', PACKAGE = 'FlexExperiment', nn_idx, self_loops)
}

