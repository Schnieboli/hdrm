#include <RcppArmadillo.h>
using namespace Rcpp;


// exact versions -----------------------------------------------------------
// [[Rcpp::export]]
double A1_i_cpp(arma::mat& mat){
  int N = mat.n_cols;
  int d = mat.n_rows;
  arma::vec vec_l1(d), vec_l2(d), diff (d);
  
  double out = 0.0;
  
  for(int l2 = 0; l2 < N-1; ++l2){
    vec_l2 = mat.col(l2);
    for(int l1 = l2+1; l1 < N; ++l1){
      vec_l1 = mat.col(l1);
      diff = vec_l1 - vec_l2;
      out += arma::dot(diff, diff);
    }
  }
  return out / (N * (N - 1));
}


//[[Rcpp::export]]
double A2_ir_cpp(arma::mat& mat1, arma::mat& mat2){
  int n1 = mat1.n_cols, n2 = mat2.n_cols, d = mat1.n_rows;
  double out = 0.0;
  arma::vec col_l2(d), diff_l(d), col_k2(d), diff_k(d);
  
  for(int l2 = 0; l2 < n1-1; ++l2){
    col_l2 = mat1.col(l2);
    for(int l1 = l2+1; l1 < n1; ++l1){
      diff_l = mat1.col(l1) - col_l2;
      for(int k2 = 0; k2 < n2-1; ++k2){
        col_k2= mat2.col(k2);
        for(int k1 = k2+1; k1 < n2; ++k1){
          diff_k = mat2.col(k1) - col_k2;
          out += pow(arma::dot(diff_l, diff_k), 2);
        }
      }
    }
  }
  return out / (4 * R::choose(n1, 2) * R::choose(n2, 2));
}


// [[Rcpp::export]]
double A3_i_cpp(arma::mat& mat){
  int n = mat.n_cols;
  double out;
  double Part1 = 0.0, Part2 = 0.0, Part3 = 0.0, Part4 = 0.0, Part5 = 0.0, Part7 = 0.0;
  double Part6 = n * n * arma::accu(arma::square(arma::mean(mat, 1)));
  double a12 = 0.0, a22 = 0.0, a13 = 0.0, a23 = 0.0;
  
  for(int l2 = 0; l2 < n; ++l2){
    a22 = arma::dot(mat.col(l2), mat.col(l2));
    Part7 += a22;
    for(int l1 = 0; l1 < n; ++l1){
      a12 = arma::dot(mat.col(l1), mat.col(l2));
      Part1 += a12*a12 * (l1!=l2);
      for(int l3 = 0; l3 < n; ++l3){
        a23 = arma::dot(mat.col(l2), mat.col(l3));
        a13 = arma::dot(mat.col(l1), mat.col(l3));
        
        Part2 += a12 * a13 * (l1 != l2) * (l2!=l3) * (l1!=l3);
        Part3 += a13 * (a23 + a12) * (l2!=l3) * (l1!=l3);
        Part4 += a13 * a22 * (l1 != l2) *(l2!=l3) * (l1!=l3);
        Part5 += a12 * a23 * (l1!=l2);
      }
    }
  }
  
  Part1 *= (n-2) * (n-3);
  Part2 *= (2*n) - 5;
  
  out = (Part1 - Part2 - Part3 - Part4 - Part5 + (Part6 * (Part6 - Part7)));
  out /= n * (n-1) * (n-2) * (n-3);
  
  return(out);
}

// subsampling versions -------------------------------------------------------

// [[Rcpp::export]]
double A1star_i_cpp(const arma::mat& mat, int& B){
  int n = mat.n_cols;
  arma::uvec ind(2);
  double out = 0.0;
  
  for(int b = 0; b < B; ++b){
    ind = arma::randperm(n).head(2);
    out += arma::dot(mat.col(ind(0)) - mat.col(ind(1)),
                     mat.col(ind(0)) - mat.col(ind(1)));
  }
  return(out/(2*B));
}


// [[Rcpp::export]]
double A2star_ir_cpp(const arma::mat& mat1, arma::mat& mat2, int B){
  int n1 = mat1.n_cols;
  int n2 = mat2.n_cols;
  arma::uvec ind1(2), ind2(2);
  double out = 0.0, tmp;
  
  for(int b = 0; b < B; ++b){
    ind1 = arma::randperm(n1).head(2);
    ind2 = arma::randperm(n2).head(2);
    
    tmp = arma::dot(mat1.col(ind1(0)) - mat1.col(ind1(1)),
                    mat2.col(ind2(0)) - mat2.col(ind2(1)));
    out += pow(tmp, 2);
  }
  return(out/(4*B));
}


// [[Rcpp::export]]
double A3star_i_cpp(const arma::mat& mat, int& B){
  int n = mat.n_cols;
  arma::uvec ind(4);
  double out = 0.0, tmp;
  
  for(int b = 0; b < B; ++b){
    ind = arma::randperm(n).head(4);
    tmp = arma::dot(mat.col(ind(0)) - mat.col(ind(1)),
                    mat.col(ind(2)) - mat.col(ind(3)));
    out += pow(tmp, 2);
  }
  return(out/(4*B));
}


// [[Rcpp::export]]
double C5star_cpp(const Rcpp::List& X_list, const int B) {
  int a = X_list.size();
  
  std::vector <arma::mat> mats(a);
  for(int i = 0; i < a; ++i){
    mats[i] = Rcpp::as<arma::mat>(X_list[i]);
  }
  int d = mats[0].n_rows;
  
  double out = 0.0;
  arma::vec Z12(d), Z34(d), Z56(d);
  arma::mat sigma(d, 6 * a);
  arma::uvec ind(6);
  
  for (int b = 0; b < B; ++b) {
    Z12.zeros();
    Z34.zeros();
    Z56.zeros();
    
    for (int i = 0; i < a; ++i) {
      arma::mat Xi = mats[i];
      int n_i = Xi.n_cols;
      ind = arma::randperm(n_i).head(6);
      
      for (int j = 0; j < 6; ++j) {
        sigma.col(6 * i + j) = Xi.col(ind(j));
      }
      
      Z12 += sigma.col(0 + 6 * i) - sigma.col(1 + 6 * i);
      Z34 += sigma.col(2 + 6 * i) - sigma.col(3 + 6 * i);
      Z56 += sigma.col(4 + 6 * i) - sigma.col(5 + 6 * i);
    }

    double dot_12_34 = arma::dot(Z12, Z34);
    double dot_34_56 = arma::dot(Z34, Z56);
    double dot_56_12 = arma::dot(Z56, Z12);
    
    out += dot_12_34 * dot_34_56 * dot_56_12;
  }
  return out / (8 * B);
}
