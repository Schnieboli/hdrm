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
    for(int l1 = l2 + 1; l1 < n; ++l1){
      col1 = mat.col(l1);
      for(int k2 = 0; k2 < (n-1); ++k2){
        col4 = mat.col(k2);
        // only call loop if condition is true
        if(k2 != l1 && k2!= l2){
          for(int k1 = k2 + 1; k1 < n; ++k1){
            // only do calculations if condition is true
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


// [[Rcpp::export()]]
double C1_eq_cpp(arma::mat &mat){
  int n = mat.n_cols;
  int d = mat.n_rows;
  double tmp1 = 0;
  double tmp2 = 0;
  double tmp3 = 0;
  double out = 0.0;
  arma::vec v1(d), v2(d), v3(d), v4(d), v5(d), v6(d);

  for(int l1 = 0; l1 < n; ++l1){
    v1 = mat.col(l1);
    for(int l2 = 0; l2 < n; ++l2){
      v2 = mat.col(l2);
      for(int l3 = 0; l3 < n; ++l3){
        v3 = mat.col(l3);
        for(int l4 = 0; l4 < n; ++l4){
          v4 = mat.col(l4);
          for(int l5 = 0; l5 < n; ++l5){
            v5 = mat.col(l5);
            for(int l6 = 0; l6 < n; ++l6){
              if(l1 != l2 && l1 != l3 && l1 != l4 && l1 != l5 && l1 != l6 && l2 != l3 && l2 != l4 && l2 != l5 && l2 != l6 && l3 != l4 && l3 != l5 && l3 != l6 && l4 != l5 && l4 != l6 && l5 != l6){
                v6 = mat.col(l6);
                tmp1 = 0;
                tmp2 = 0;
                tmp3 = 0;
                for(int j = 0; j < d; ++j){
                  tmp1 += (v1(j) - v2(j)) * (v3(j) - v4(j));
                  tmp2 += (v3(j) - v4(j)) * (v5(j) - v6(j));
                  tmp3 += (v5(j) - v6(j)) * (v1(j) - v2(j));
                }
                out+= tmp1 * tmp2 * tmp3;
              }
            }
          }
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export()]]
double C1_star_eq_cpp(arma::mat &mat, int B){
  int n = mat.n_cols;
  int d = mat.n_rows;
  double tmp1 = 0;
  double tmp2 = 0;
  double tmp3 = 0;
  double out = 0.0;
  arma::vec v1(d), v2(d), v3(d), v4(d), v5(d), v6(d);
  arma::uvec ind(6);

  for(int b = 0; b < B; ++b){
    ind = arma::randperm(n, 6);
    v1 = mat.col(ind(0));
    v2 = mat.col(ind(1));
    v3 = mat.col(ind(2));
    v4 = mat.col(ind(3));
    v5 = mat.col(ind(4));
    v6 = mat.col(ind(5));
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    for(int j = 0; j < d; ++j){
      tmp1 += (v1(j) - v2(j)) * (v3(j) - v4(j));
      tmp2 += (v3(j) - v4(j)) * (v5(j) - v6(j));
      tmp3 += (v5(j) - v6(j)) * (v1(j) - v2(j));
    }
    out += tmp1 * tmp2 * tmp3;
  }
  return out;
}
