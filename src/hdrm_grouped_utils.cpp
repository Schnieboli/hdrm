#include <RcppArmadillo.h>
using namespace Rcpp;


// exact versions -----------------------------------------------------------
// [[Rcpp::export]]
double A1_cpp(arma::mat& mat){
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

// [[Rcpp::export]]
double A3_cpp(arma::mat& mat){
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
double A1star_cpp(const arma::mat& mat, int& B){
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
double A2star_cpp(const arma::mat& mat1, arma::mat& mat2, int& B){
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
double A3star_cpp(const arma::mat& mat, int& B){
  int n = mat.n_cols;
  int d = mat.n_rows;
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
double C5star_cpp_internal(arma::mat& X, arma::vec& group, const int& B, arma::uvec& n){ // Matrix X ist schon mit TM multipliziert und schon mit sqrt(N/n) multipliziert
  // ausserdem muss X nach Gruppen sortiert sein!!!
  int a = unique(group).index_max() + 1;
  int m = X.n_rows;
  double cout = 0.0;
  arma::vec Z12(m), Z34(m), Z56(m);
  arma::mat sigma(m, 6*a);
  arma::uvec indizes(6);
  int ind = 0;
  
  for(int b = 0; b < B; ++b){
    Z12.zeros();
    Z34.zeros();
    Z56.zeros();
    int shift = 0;
    for(int i = 0; i < a; ++i){
      indizes = arma::randperm(n(i)).head(6); // einfach so lassen!!!
      for(int j = 0; j < 6; ++j){
        ind = shift + indizes(j);
        sigma.col(6*i + j) = X.col(ind);
      }
      shift += n(i); // damit immer die richtige Gruppe ausgewaehlt wird...
      Z12 += sigma.col(0 + 6*i) - sigma.col(1 + 6*i);
      Z34 += sigma.col(2 + 6*i) - sigma.col(3 + 6*i);
      Z56 += sigma.col(4 + 6*i) - sigma.col(5 + 6*i);
    }
    cout += arma::accu(Z12.t() * Z34) * arma::accu(Z34.t() * Z56) * arma::accu(Z56.t() * Z12);
  }
  return cout/(8*B);
}
