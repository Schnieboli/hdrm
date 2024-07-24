#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double B3_cpp(const NumericMatrix& mat) {
    int nrows = mat.nrow();
    int ncols = mat.ncol();
    double total_sum = 0.0;
    double dot_product12 = 0.0;
    double dot_product23 = 0.0;
    double dot_product31 = 0.0;
    double val1 = 0.0;
    double val2 = 0.0;
    double val3 = 0.0;

    for (int i = 0; i < nrows - 2; ++i) {
      for (int j = i + 1; j < nrows - 1; ++j) {
        for (int k = j + 1; k < nrows; ++k) {
          dot_product12 = 0.0;
          dot_product23 = 0.0;
          dot_product31 = 0.0;

          for (int col = 0; col < ncols; ++col) {
            val1 = mat(i, col);
            val2 = mat(j, col);
            val3 = mat(k, col);

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
  int nrows = mat.nrow();
  int ncols = mat.ncol();

  double total_sum = 0.0;
  double temp = 0.0;

  double c1 = 0.0;
  double c2 = 0.0;

  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < nrows; ++j) {
      temp = 0.0;
      if(i != j){

        for(int k = 0; k < ncols; ++k){

          c1 = mat(i,k);
          c2 = mat(j,k);

          temp += c1 * c2;

        }
        total_sum += pow(temp, 2);
      }
    }
  }
  return total_sum;
}
