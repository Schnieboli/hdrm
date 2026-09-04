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
  return out/N;
}


// [[Rcpp::export]]
double B2_cpp(arma::mat& mat) {
  int d = mat.n_rows;
  int N = mat.n_cols;

  double out = 0.0;
  double temp = 0.0;

  double c1 = 0.0;
  double c2 = 0.0;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      temp = 0.0;
      if(i != j){
        for(int k = 0; k < d; ++k){
          c1 = mat(k,i);
          c2 = mat(k,j);
          temp += c1 * c2;
        }
        out += pow(temp, 2);
      }
    }
  }
  return out / (N*(N-1));
}


// [[Rcpp::export]]
double B3_cpp(arma::mat& mat) {
  int d = mat.n_rows;
  int N = mat.n_cols;
  double out = 0.0;
  double dot_product12 = 0.0;
  double dot_product23 = 0.0;
  double dot_product31 = 0.0;
  arma::vec vec1(d);
  arma::vec vec2(d);
  arma::vec vec3(d);

  for (int i = 0; i < N - 2; ++i) {
    for (int j = i + 1; j < N - 1; ++j) {
      for (int k = j + 1; k < N; ++k) {
        dot_product12 = 0.0;
        dot_product23 = 0.0;
        dot_product31 = 0.0;

        vec1 = mat.col(i);
        vec2 = mat.col(j);
        vec3 = mat.col(k);

        for (int r = 0; r < d; ++r) {

          dot_product12 += vec1(r) * vec2(r);
          dot_product23 += vec2(r) * vec3(r);
          dot_product31 += vec3(r) * vec1(r);
        }

        out += dot_product12 * dot_product23 * dot_product31;
      }
    }
  }
  return out;
}
