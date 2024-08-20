#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;



// // [[Rcpp::export]]
// arma::mat A2_cpp(arma::mat& mat_i, arma::mat& mat_r, arma::mat& P_i, arma::mat& P_r) {
//
//   arma::mat M_i = (((P_i * mat_i.t()) * mat_r) * P_r.t());
//   arma::mat out = M_i % M_i;
//
//   return out;
// }


// // [[Rcpp::export]]
// double C5star_cpp(arma::mat& X, arma::vec& group, const int B, arma::mat& TM, arma::uvec& n){
//
//   int d = X.n_rows;
//   int a = unique(group).index_max() + 1;
//
//   double cout = 0.0;
//   arma::vec Z12(d*a), Z34(d*a), Z56(d*a);
//   arma::mat sigma(d, 6*a);
//   arma::vec indizes(6);
//   int ind = 0;
//
//   for(int b = 0; b < B; ++b){
//
//     for(int i = 0; i < a; ++i){
//       indizes = arma::conv_to<arma::vec>::from(arma::randperm(n(i)).head(6)); // einfach so lassen!!!
//       for(int j = 0; j < 6; ++j){
//         ind = indizes(j);
//         sigma.col(6*i + j) = X.col(ind);
//       }
//       for(int i = 0; i < a; ++i){
//         for(int j = 0; j < d; ++j){
//
//           Z12(j + d*i) = sigma(j, 0 + 6*i) - sigma(j, 1 + 6*i);
//           Z34(j + d*i) = sigma(j, 2 + 6*i) - sigma(j, 3 + 6*i);
//           Z56(j + d*i) = sigma(j, 4 + 6*i) - sigma(j, 5 + 6*i);
//
//
//         }
//       }
//       cout += arma::accu(Z12.t() * (TM * Z34)) * arma::accu(Z34.t() * (TM * Z56)) * arma::accu(Z56.t() * (TM * Z12));
//     }
//   }
//
//   return cout/(8*B);
// }
//
//
// arma::mat A2_cpp(arma::mat& mat_i, arma::mat& mat_r, arma::mat& P_i, arma::mat& P_r) {
//
//   arma::mat M_i = (((P_i * mat_i.t()) * mat_r) * P_r.t());
//   arma::mat out = M_i % M_i;
//
//   return out;
// }


// [[Rcpp::export]]
double C5star_cpp_neu(arma::mat& X, arma::vec& group, const int B, arma::uvec& n){ // Matrix X ist schon mit TM multipliziert und schon mit sqrt(N/n) multipliziert
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

    for(int i = 0; i < a; ++i){
      indizes = arma::conv_to<arma::vec>::from(arma::randperm(n(i)).head(6)); // einfach so lassen!!!
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
