#include <Rcpp.h>
#include <RcppCommon.h>
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export(rng = false)]]
Rcpp::List rcpp_knn_to_sparse_weights(
    Rcpp::IntegerMatrix & nn_idx,
    Rcpp::NumericMatrix & nn_dist,
    bool self_loops
) {
  int n_obs = nn_idx.nrow();
  int n_neighbors = nn_idx.ncol();
  int max_len = n_obs * n_neighbors;
  Rcpp::IntegerVector row_idx(max_len);
  Rcpp::IntegerVector col_idx(n_obs * n_neighbors);
  Rcpp::NumericVector vals(max_len);

  int idx = 0;
  for (int i = 0; i < n_obs; ++i) {
    for (int j = 0; j < n_neighbors; ++j) {
      if (nn_idx(i, j) == 0) continue;
      if ((nn_idx(i, j) == i + 1) & !self_loops) continue;
      col_idx[idx] = i;
      row_idx[idx] = nn_idx(i, j) - 1;
      vals[idx] = nn_dist(i, j);
      idx += 1;
    }
  }
  if (idx < max_len) {
    row_idx = Rcpp::head(row_idx, idx);
    col_idx = Rcpp::head(col_idx, idx);
    vals = Rcpp::head(vals, idx);
  }
  Rcpp::List res = Rcpp::List::create(
    Rcpp::_["i"] = row_idx,
    Rcpp::_["j"] = col_idx,
    Rcpp::_["x"] = vals
  );
  return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List rcpp_knn_to_sparse_idx(
    Rcpp::IntegerMatrix & nn_idx,
    bool self_loops
) {
  int n_obs = nn_idx.nrow();
  int n_neighbors = nn_idx.ncol();
  int max_len = n_obs * n_neighbors;
  Rcpp::IntegerVector row_idx(max_len);
  Rcpp::IntegerVector col_idx(max_len);

  int idx = 0;
  for (int i = 0; i < n_obs; ++i) {
    for (int j = 0; j < n_neighbors; ++j) {
      if (nn_idx(i, j) == 0) continue;
      if ((nn_idx(i, j) == i + 1) & !self_loops) continue;
      col_idx[idx] = i;
      row_idx[idx] = nn_idx(i, j) - 1;
      idx += 1;
    }
  }
  if (idx < max_len) {
    row_idx = Rcpp::head(row_idx, idx);
    col_idx = Rcpp::head(col_idx, idx);
  }
  Rcpp::List res = Rcpp::List::create(
    Rcpp::_["i"] = row_idx,
    Rcpp::_["j"] = col_idx,
    Rcpp::_["x"] = 1
  );
  return res;
}
