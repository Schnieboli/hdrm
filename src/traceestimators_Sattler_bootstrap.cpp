#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double A1star_cpp(const arma::mat& X, int& B){
  int n = X.n_cols;
  int d = X.n_rows;
  arma::vec v1(d), v2(d);
  arma::uvec ind(2);
  double out = 0.0;
  double temp = 0.0;

  for(int b = 0; b < B; ++b){
    temp = 0.0;
    ind = arma::conv_to<arma::uvec>::from(arma::randperm(n).head(2));
    v1 = X.col(ind(0));
    v2 = X.col(ind(1));
    for(int j = 0; j < d; ++j){
      temp += pow((v1(j) - v2(j)), 2);
    }
    out += temp;
  }
return(out/(2*B));
}


// [[Rcpp::export]]
double A2star_cpp(const arma::mat& X, arma::mat& Y, int& B){
  int nX = X.n_cols;
  int nY = Y.n_cols;
  int d = X.n_rows;
  arma::vec v1(d), v2(d), v3(d), v4(d);
  arma::uvec indX(2), indY(2);
  double out = 0.0;
  double temp = 0.0;

  for(int b = 0; b < B; ++b){
    temp = 0.0;
    indX = arma::conv_to<arma::uvec>::from(arma::randperm(nX).head(2));
    indY = arma::conv_to<arma::uvec>::from(arma::randperm(nY).head(2));
    v1 = X.col(indX(0));
    v2 = X.col(indX(1));
    v3 = Y.col(indY(0));
    v4 = Y.col(indY(1));
    for(int j = 0; j < d; ++j){
      temp += pow((v1(j) - v2(j)) * (v3(j) - v4(j)), 2);
    }
    out += temp;
  }
  return(out/(4*B));
}



// [[Rcpp::export]]
double A3star_cpp(const arma::mat& X, int& B){
  int n = X.n_cols;
  int d = X.n_rows;
  arma::vec v1(d), v2(d), v3(d), v4(d);
  arma::uvec ind(2);
  double out = 0.0;
  double temp = 0.0;

  for(int b = 0; b < B; ++b){
    temp = 0.0;
    ind = arma::conv_to<arma::uvec>::from(arma::randperm(n).head(4));
    v1 = X.col(ind(0));
    v2 = X.col(ind(1));
    v3 = X.col(ind(2));
    v4 = X.col(ind(3));
    for(int j = 0; j < d; ++j){
      temp += pow((v1(j) - v2(j)) * (v3(j) - v4(j)), 2);
    }
    out += temp;
  }
  return(out/(4*B));
}

// [[Rcpp::export]]
double C5star_cpp_internal(arma::mat& X, arma::vec& group, const int B, arma::uvec& n){ // Matrix X ist schon mit TM multipliziert und schon mit sqrt(N/n) multipliziert
  // außerdem muss X nach Gruppen sortiert sein!!!
  int a = unique(group).index_max() + 1;
  int d = X.n_rows/a;
  double cout = 0.0;
  arma::vec Z12(d*a), Z34(d*a), Z56(d*a);
  arma::mat sigma(d*a, 6*a);
  arma::vec indizes(6);
  int ind = 0;

  for(int b = 0; b < B; ++b){
    arma::vec Z12(d*a), Z34(d*a), Z56(d*a);
    int shift = 0;
    for(int i = 0; i < a; ++i){
      indizes = shift + arma::conv_to<arma::vec>::from(arma::randperm(n(i)).head(6)); // einfach so lassen!!!
      shift += n(i); // damit immer die richtige Gruppe ausgewählt wird...
      for(int j = 0; j < 6; ++j){
        ind = indizes(j);
        sigma.col(6*i + j) = X.col(ind);
      }
      Z12 += sigma.col(0 + 6*i) - sigma.col(1 + 6*i);
      Z34 += sigma.col(2 + 6*i) - sigma.col(3 + 6*i);
      Z56 += sigma.col(4 + 6*i) - sigma.col(5 + 6*i);
    }
    cout += arma::accu(Z12.t() * Z34) * arma::accu(Z34.t() * Z56) * arma::accu(Z56.t() * Z12);
  }
  return cout/(8*B);
}

