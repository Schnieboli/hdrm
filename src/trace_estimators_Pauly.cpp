#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double B3_cpp(const NumericMatrix& mat) {
    int d = mat.nrow();
    int N = mat.ncol();
    double total_sum = 0.0;
    double dot_product12 = 0.0;
    double dot_product23 = 0.0;
    double dot_product31 = 0.0;
    double val1 = 0.0;
    double val2 = 0.0;
    double val3 = 0.0;

    for (int i = 0; i < N - 2; ++i) {
      for (int j = i + 1; j < N - 1; ++j) {
        for (int k = j + 1; k < N; ++k) {
          dot_product12 = 0.0;
          dot_product23 = 0.0;
          dot_product31 = 0.0;

          for (int r = 0; r < d; ++r) {
            val1 = mat(r, i);
            val2 = mat(r, j);
            val3 = mat(r, k);

            dot_product12 += val1 * val2;
            dot_product23 += val2 * val3;
            dot_product31 += val3 * val1;
          }

          total_sum += dot_product12 * dot_product23 * dot_product31;
        }
      }
    }
    return total_sum;
  }

// [[Rcpp::export]]
double B2_cpp(const NumericMatrix& mat) {
  int d = mat.nrow();
  int N = mat.ncol();

  double total_sum = 0.0;
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
        total_sum += pow(temp, 2);
      }
    }
  }
  return total_sum / (N*(N-1));
}
