#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export()]]
double A1_eq_cpp(arma::mat &mat){
  int d = mat.n_rows;
  int n = mat.n_cols;
  double out = 0.0;
  arma::vec col_l1(d);
  
  for(int l1 = 0; l1 < (n-1); ++l1){
    col_l1 = mat.col(l1);
    for(int l2 = l1+1; l2 < n; ++l2){
      out += arma::accu(arma::square(col_l1 - mat.col(l2)));
    }
  }
  return out;
}

// [[Rcpp::export()]]
double A2_eq_cpp(arma::mat &mat){
  int d = mat.n_rows;
  int n = mat.n_cols;
  double out = 0.0;
  arma::vec col_l1(d), col_l2(d), col_k1(d), col_k2(d);
  
  for (int l2 = 0; l2 < (n-1); ++l2) {
    col_l2 = mat.col(l2);
    for(int l1 = l2 + 1; l1 < n; ++l1){
      col_l1 = mat.col(l1);
      for(int k2 = 0; k2 < (n-1); ++k2){
        col_k2 = mat.col(k2);
        // only call loop if condition is true
        if(k2 != l1 && k2!= l2){
          for(int k1 = k2 + 1; k1 < n; ++k1){
            // only do calculations if condition is true
            if(k1 != l1 && k1 != l2){
              out += pow(arma::dot(col_l2 - col_l1, 
                                   col_k2 - mat.col(k1)),
                                   2);
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
  double out = 0.0;
  arma::vec col_l1(d), diff_12(d), col_l3(d), diff_34(d), col_l5(d), diff_56(d);
  double prod_1234 = 0.0;
  
  for(int l1 = 0; l1 < n; ++l1){
    col_l1 = mat.col(l1);
    for(int l2 = 0; l2 < n; ++l2){
      if(l1 == l2) continue;
      diff_12 = col_l1 - mat.col(l2);
      for(int l3 = 0; l3 < n; ++l3){
        if((l1 == l3) | (l2 == l3)) continue;
        col_l3 = mat.col(l3);
        for(int l4 = 0; l4 < n; ++l4){
          if((l1 == l4) | (l2 == l4) | (l3 == l4)) continue;
          diff_34 = col_l3 - mat.col(l4);
          prod_1234 = arma::dot(diff_12, diff_34);
          for(int l5 = 0; l5 < n; ++l5){
            if((l1 == l5) | (l2 == l5) | (l3 == l5) | (l4 == l5)) continue;
            col_l5 = mat.col(l5);
            for(int l6 = 0; l6 < n; ++l6){
              if((l1 == l6) | (l2 == l6) | (l3 == l6) | (l4 == l6) | (l5 == l6))
                continue;
              diff_56 = col_l5 - mat.col(l6);
              out += prod_1234 * 
                arma::dot(diff_34, diff_56) * 
                arma::dot(diff_56, diff_12);
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
  double out = 0.0;
  arma::vec diff_12(d), diff_34(d), diff_56(d);
  arma::uvec ind(6);
  
  for(int b = 0; b < B; ++b){
    ind = arma::randperm(n).head(6);
    diff_12 = mat.col(ind(0)) - mat.col(ind(1));
    diff_34 = mat.col(ind(2)) - mat.col(ind(3));
    diff_56 = mat.col(ind(4)) - mat.col(ind(5));
    
    out += arma::dot(diff_12, diff_34) * 
      arma::dot(diff_34, diff_56) * 
      arma::dot(diff_56, diff_12);
  }
  return out;
}
