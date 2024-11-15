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

// [[Rcpp::export()]]
double A2_eq_cpp(arma::mat &mat){
  int d = mat.n_rows;
  int n = mat.n_cols;
  double tmp = 0.0;
  double out = 0.0;
  arma::vec col1(d), col2(d), col3(d), col4(d);

  for (int l2 = 0; l2 < (n-1); ++l2) {
    col2 = mat.col(l2);

    for(int l1 = l2+1; l1 < n; ++l1){
      col1 = mat.col(l1);

      for(int k2 = 1; k2 < (n-1); ++k2){
        col4 = mat.col(k2);

        if(k2 != l1 && k2!= l2){

          for(int k1 = k2 + 1; k1 < n; ++k1){

            if(k1 != l1 && k1 != l2){
              col3 = mat.col(k1);
              tmp = 0;

              for(int j = 0; j < d; ++j){
                tmp += (col1(j) - col2(j)) * (col3(j) - col4(j));
              }

              out += pow(tmp, 2);
            }

          }
        }

      }
    }
  }
  return out;
}
