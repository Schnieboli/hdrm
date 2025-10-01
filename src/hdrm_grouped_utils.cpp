#include <RcppArmadillo.h>
using namespace Rcpp;


// exact versions -----------------------------------------------------------
// [[Rcpp::export]]
double A1_cpp(arma::mat& mat){
  int N = mat.n_cols;
  int d = mat.n_rows;
  arma::colvec col1(d), col2(d);

  double out = 0.0;
  double diff = 0.0;

  // Fall 1: Wenn die Matrix mehr als eine Zeile hat, wie im Originalcode
  if(d > 1){
    for(int i = 0; i < N - 1; ++i){
      col1 = mat.col(i);
      for(int j = i + 1; j < N; ++j){
        col2 = mat.col(j);
        for(int k = 0; k < d; ++k){
          diff = col1(k) - col2(k);
          out += diff * diff;
        }
      }
    }
    return out / (N * (N - 1));
  }

  // Fall 2: Wenn die Matrix nur eine Zeile hat
  else {
    for(int i = 0; i < N - 1; ++i){
      col1 = mat.col(i);
      for(int j = i + 1; j < N; ++j){
        col2 = mat.col(j);
        diff = col1(0) - col2(0); // Nur die Werte der ersten (und einzigen) Zeile vergleichen
        out += diff * diff;
      }
    }
    return out / (N * (N - 1));
  }
}


// [[Rcpp::export]]
double A3_cpp(arma::mat& mat, double Part6){ // Part6 = sum(rowMeans(mat)^2) in R

  int n = mat.n_cols;
  int d = mat.n_rows;
  double out;
  double Part1 = 0.0, Part2 = 0.0, Part3 = 0.0, Part4 = 0.0, Part5 = 0.0, Part7 = 0.0;
  double a12 = 0.0, a22 = 0.0, a13 = 0.0, a23 = 0.0;


  for(int l2 = 0; l2 < n; ++l2){
    a22 = 0;
    for(int i = 0; i< d; ++i){
      a22 += pow(mat(i,l2), 2);
    }
    Part7 += a22;

    for(int l1 = 0; l1 < n; ++l1){
      a12 = 0;
      for(int i = 0; i < d; ++i){
        a12 += mat(i,l1) * mat(i,l2);
      }
      Part1 += a12*a12 * (l1!=l2);
      for(int l3 = 0; l3 < n; ++l3){
        a23 = 0;
        a13 = 0;
        for(int i = 0; i < d; ++i){
          a23 += mat(i,l2) * mat(i,l3);
          a13 += mat(i,l1) * mat(i,l3);
        }
        Part5 += a12 * a23 * (l1!=l2);
        Part2 += a12 * a13 * (l1 != l2) * (l2!=l3) * (l1!=l3);
        Part3 += a13 * (a23 + a12) * (l2!=l3) * (l1!=l3);
        Part4 += a13 * a22 * (l1 != l2) *(l2!=l3) * (l1!=l3);
      }
    }
  }

  Part1 *= (n-2) * (n-3);
  Part2 *= (2*n) - 5;
  Part6 *= n*n;

  out = (Part1 - Part2 - Part3 - Part4 - Part5 + (Part6 * (Part6 - Part7)));
  out /= n * (n-1) * (n-2) * (n-3);

  return(out);
}

// bootstrap versions -------------------------------------------------------

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
    ind = arma::randperm(n, 2);
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
    indX = arma::randperm(nX, 2);
    indY = arma::randperm(nY, 2);
    v1 = X.col(indX(0));
    v2 = X.col(indX(1));
    v3 = Y.col(indY(0));
    v4 = Y.col(indY(1));
    for(int j = 0; j < d; ++j){
      temp += (v1(j) - v2(j)) * (v3(j) - v4(j));
    }
    out += pow(temp, 2);
  }
  return(out/(4*B));
}



// [[Rcpp::export]]
double A3star_cpp(const arma::mat& X, int& B){
  int n = X.n_cols;
  int d = X.n_rows;
  arma::colvec v1(d), v2(d), v3(d), v4(d);
  arma::uvec ind(2);
  double out = 0.0;
  double temp = 0.0;

  for(int b = 0; b < B; ++b){
    temp = 0.0;
    ind = arma::randperm(n, 4);
    v1 = X.col(ind(0));
    v2 = X.col(ind(1));
    v3 = X.col(ind(2));
    v4 = X.col(ind(3));
    for(int j = 0; j < d; ++j){
      temp += (v1(j) - v2(j)) * (v3(j) - v4(j));
    }
    out += pow(temp, 2);
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
      indizes = arma::randperm(n(i), 6); // einfach so lassen!!!
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
