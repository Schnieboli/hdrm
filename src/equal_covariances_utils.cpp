#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export()]]
double A1_eq_cpp(arma::mat &mat){
  int d = mat.n_rows;
  int n = mat.n_cols;
  double tmp = 0.0;
  double out = 0.0;
  arma::vec col1(d), col2(d);

  for(int l2 = 0; l2 < (n-1); ++l2){
    col1 = mat.col(l2);
    for(int l1 = l2+1; l1 < n; ++l1){
      col2 = mat.col(l1);
      for(int j = 0; j < d; ++j){
        tmp = col2(j) - col1(j);
        out += pow(tmp, 2);
      }
    }
  }
  return out;
}

// TO-DO: A2_eq_cpp
