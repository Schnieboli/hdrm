#include <Rcpp.h>
using namespace Rcpp;
//
// // [[Rcpp::export]]
// double A1_cpp(const NumericMatrix& mat){
//   int nrows = mat.nrow();
//   int ncols = mat.ncol();
//
//   double out = 0.0;
//   double diff = 0.0;
//   double val1 = 0.0;
//   double val2 = 0.0;
//
//
//   for(int i = 0; i < nrows - 1; ++i){
//     for(int j = i + 1; j < nrows; ++j){
//
//       for(int k = 0; k < ncols; ++k){
//
//         val1 = mat(i,k);
//         val2 = mat(j,k);
//         diff = val1 - val2;
//
//         out += diff*diff;
//       }
//
//     }
//   }
// return out;
// }
//
//
// // [[Rcpp::export]]
// double A2_cpp(const NumericMatrix& mat_i, const NumericMatrix& mat_r){
//   int n1 = mat_i.nrow();
//   int n2 = mat_r.nrow();
//   int d = mat_i.ncol();
//
//
//   double temp = 0.0;
//   double out = 0.0;
//   double p1;
//   double p2;
//
//   for(int i = 0; i < n1; ++i){
//     for(int j = 0; j < n2; ++j){
//       temp = 0;
//       for(int k = 0; k < d; ++k){
//         p1 = mat_i(i,k);
//         p2 = mat_r(j,k);
//         temp += p1 * p2;
//
//       }
//       out += pow(temp, 2);
//     }
//   }
// return out;
// }

// [[Rcpp::export]]
double A3_cpp(const NumericMatrix& mat, double Part6){ // gleiches Ergebnis wie A3_R -> etwa 12x schneller

  int n = mat.ncol(), d = mat.nrow();
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
      Part1 += a12*a12 *(l1!=l2);
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
  Part2 *= (2*n - 5);
  Part6 *= pow(n,2);

  out = (Part1 - Part2 - Part3 - Part4 - Part5 + (Part6 * (Part6 - Part7))) / (n * (n-1) * (n-2)* (n-3));

  return(out);

}








