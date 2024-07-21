#library(Rcpp)


#' @export
Rcpp::cppFunction(' // etwa 10x schneller
double B3_cpp(NumericMatrix mat) {
  int nrows = mat.nrow();
  double total_sum = 0.0;

  for (int i = 0; i < nrows - 2; ++i) {
    for (int j = i + 1; j < nrows - 1; ++j) {
      for (int k = j + 1; k < nrows; ++k) {
        NumericVector row1 = mat(i, _);
        NumericVector row2 = mat(j, _);
        NumericVector row3 = mat(k, _);

        double dot_product12 = std::inner_product(row1.begin(), row1.end(), row2.begin(), 0.0);
        double dot_product23 = std::inner_product(row2.begin(), row2.end(), row3.begin(), 0.0);
        double dot_product31 = std::inner_product(row3.begin(), row3.end(), row1.begin(), 0.0);

        total_sum += dot_product12 * dot_product23 * dot_product31;
      }
    }
  }
  return total_sum;
}
')

#' @export
Rcpp::cppFunction(' // etwa 2x schneller
double B2_cpp(NumericMatrix mat) {
  int nrows = mat.nrow();
  double total_sum = 0.0;
  double temp = 0.0;

  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < nrows; ++j) {

        if(i != j){
        NumericVector row1 = mat(i, _);
        NumericVector row2 = mat(j, _);
        temp = std::inner_product(row1.begin(), row1.end(), row2.begin(), 0.0);
        total_sum += pow(temp, 2);
        }
      }
    }
  return total_sum;
}
')

@export
Rcpp::cppFunction(' // etwa weniger als 2x schneller
double B0_cpp(NumericMatrix mat) {
  int nrows = mat.nrow();
  double total_sum = 0.0;

  for (int i = 0; i < nrows; ++i) {
        NumericVector row1 = mat(i, _);
        total_sum += std::inner_product(row1.begin(), row1.end(), row1.begin(), 0.0);
    }
  return total_sum;
}
')


# library(microbenchmark)
# N <- 20
# d <- 30
# X <- matrix(rnorm(N*d), N, d)
# microbenchmark(x <- B0_cpp(X), times = 600)
# microbenchmark(x <- B0(X, N), times = 600)
