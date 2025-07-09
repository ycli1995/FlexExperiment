#include <RcppEigen.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> col_merge_dgcmatrix_cpp(
    std::vector<Eigen::SparseMatrix<double>> mat_list,
    std::vector<std::vector<std::string>> mat_rownames,
    std::vector<std::string> all_rownames
) {
  // row names map
  std::unordered_map<std::string, int> name_map;
  for (int i = 0; i < all_rownames.size(); i++) {
    name_map[all_rownames[i]] = i;
  }
  // row indices map for each matrix
  std::vector<std::vector<int>> mat_map;
  mat_map.reserve(mat_rownames.size());
  std::vector<int> offsets(mat_list.size(), 0);
  int num_cols = 0;
  int num_nZero = 0;
  for (int i = 0; i < mat_rownames.size(); i++) {
    std::vector<int> map(mat_rownames[i].size(), 0);
    for (int j = 0; j < mat_rownames[i].size(); j++) {
      map[j] = name_map[mat_rownames[i][j]];
    }
    mat_map.emplace_back(map);
    offsets[i] = num_cols;
    num_cols += mat_list[i].cols();
    num_nZero += mat_list[i].nonZeros();
  }

  int num_rows = all_rownames.size();
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_cols);
  combined_mat.reserve(num_nZero);
  for (int i = 0; i < mat_list.size(); ++i) {
    for (int col = 0; col < mat_list[i].cols(); ++col) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat_list[i], col); it; ++it) {
        combined_mat.insert(mat_map[i][it.index()], col + offsets[i]) = it.value();
      }
    }
  }
  combined_mat.makeCompressed();
  return combined_mat;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> row_merge_dgcmatrix_cpp(
    std::vector<Eigen::SparseMatrix<double>> mat_list,
    std::vector<std::vector<std::string>> mat_colnames,
    std::vector<std::string> all_colnames
) {
  // column names map
  std::unordered_map<std::string, int> name_map;
  for (int i = 0; i < all_colnames.size(); i++) {
    name_map[all_colnames[i]] = i;
  }
  // column indices map for each matrix
  std::vector<Eigen::SparseVector<int>> mat_map;
  mat_map.reserve(mat_colnames.size());
  std::vector<int> offsets(mat_list.size(), 0);
  int num_rows = 0;
  int num_nZero = 0;
  for (int i = 0; i < mat_colnames.size(); i++) {
    Eigen::SparseVector<int> map(all_colnames.size());
    for (int j = 0; j < mat_colnames[i].size(); j++) {
      std::string key = mat_colnames[i][j];
      if (name_map.count(key)) {
        map.insert(name_map[key], 0) = j + 1;
      }
    }
    mat_map.emplace_back(map);
    offsets[i] = num_rows;
    num_rows += mat_list[i].rows();
    num_nZero += mat_list[i].nonZeros();
  }

  int num_cols = all_colnames.size();
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_cols);
  combined_mat.reserve(num_nZero);
  for (int col = 0; col < num_cols; ++col) {
    for (int i = 0; i < mat_list.size(); ++i) {
      int c = mat_map[i].coeff(col, 0);
      if (c == 0) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat_list[i], c - 1); it; ++it) {
        combined_mat.insert(it.index() + offsets[i], col) = it.value();
      }
    }
  }
  combined_mat.makeCompressed();
  return combined_mat;
}
