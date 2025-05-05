#include <RcppEigen.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> col_merge_dgCMatrices_rcpp(
    List mat_list,
    List mat_colnames,
    std::vector<std::string> all_colnames
) {
  // Convert Rcpp lists to c++ vectors
  std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> mat_vec;
  mat_vec.reserve(mat_list.size());

  std::vector<std::vector<std::string>> rownames_vec;
  rownames_vec.reserve(mat_colnames.size());

  std::vector<std::unordered_map<std::string, int>> map_vec;
  map_vec.reserve(mat_list.size());

  int num_rows = 0;
  int num_nZero = 0;
  // offsets keep track of which row to add in to
  std::vector<int> offsets;

  for (size_t i = 0; i < mat_list.size(); i++) {
    mat_vec.emplace_back(
      Rcpp::as<Eigen::SparseMatrix<double, Eigen::ColMajor>>(mat_list.at(i))
    );
    rownames_vec.emplace_back(mat_colnames[i]);
    // Set up hash maps for rowname based lookup
    std::unordered_map<std::string, int> mat_map;
    for (size_t j = 0; j < rownames_vec[i].size(); j++) {
      mat_map[rownames_vec[i][j]] = j;
    }
    map_vec.emplace_back(mat_map);
    offsets.push_back(num_rows);
    num_rows += mat_vec[i].rows();
    num_nZero += mat_vec[i].nonZeros();
  }
  // set up tripletList for new matrix creation
  std::vector<Eigen::Triplet<double>> tripletList;
  int num_cols = all_colnames.size();
  tripletList.reserve(num_nZero);
  // loop over all rows and add nonzero entries to tripletList
  for(int i = 0; i < num_cols; i++) {
    std::string key = all_colnames[i];
    for(int j = 0; j < mat_vec.size(); j++) {
      if (map_vec[j].count(key)) {
        for(Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it1(mat_vec[j], map_vec[j][key]); it1; ++it1){
          tripletList.emplace_back(it1.row() + offsets[j], i, it1.value());
        }
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_cols);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}
