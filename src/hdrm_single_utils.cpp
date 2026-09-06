#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double B0_cpp(arma::mat& mat){
  int d = mat.n_rows;
  int N = mat.n_cols;
  double out = 0;
  double c = 0;
  arma::vec col(d);

  for(int i = 0; i < N; ++i){
    c = 0;
    col = mat.col(i);
    for(int j = 0; j < d; ++j){
      c += pow(col(j), 2);
    }
    out += c;
  }
  return out;
}


// [[Rcpp::export]]
double B2_cpp(arma::mat& mat) {
  int N = mat.n_cols;

  double out = 0.0;
  double tmp = 0.0;

  for (int k = 0; k < N; ++k) {
    for (int l = 0; l < N; ++l) {
      if(k != l){
        tmp = arma::dot(mat.col(k), mat.col(l));
        out += pow(tmp, 2);
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
double B3_cpp(arma::mat& mat) {
  int d = mat.n_rows;
  int N = mat.n_cols;
  double out = 0.0;
  arma::vec vec_k(d), vec_l(d), vec_r(d);

  for (int k = 0; k < N - 2; ++k) {
    for (int l = k + 1; l < N - 1; ++l) {
      for (int r = l + 1; r < N; ++r) {
        
        vec_k = mat.col(k);
        vec_l = mat.col(l);
        vec_r = mat.col(r);
        
        out += arma::dot(vec_k, vec_l) * 
          arma::dot(vec_l, vec_r) * 
          arma::dot(vec_r, vec_k);
      }
    }
  }
  return out;
}
