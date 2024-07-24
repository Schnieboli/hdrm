library(Rcpp)
#'@export
cppFunction(
'double B3_cpp(const NumericMatrix& mat) {
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
')

# Test the function in R




Rcpp::cppFunction(' // etwa 2x schneller
double B2_cpp(const NumericMatrix mat) {
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
')

#' #' @export
#' Rcpp::cppFunction(' // etwa weniger als 2x schneller
#' double B0_cpp(NumericMatrix mat) {
#'   int nrows = mat.nrow();
#'   double total_sum = 0.0;
#'
#'   for (int i = 0; i < nrows; ++i) {
#'         NumericVector row1 = mat(i, _);
#'         total_sum += std::inner_product(row1.begin(), row1.end(), row1.begin(), 0.0);
#'     }
#'   return total_sum;
#' }
#' ')


# library(microbenchmark)
# N <- 80
# d <- 30
# X <- matrix(rnorm(N*d), d, N)
# microbenchmark(x <- B0_cpp(t(X)), # Achtung: B0_cpp ist O(N), aber B0 mMn O(1) -> B0 ist bei N = 200, d = 30 schneller!
#                x <- B0(X),
#                times = 500,
#                unit = "relative"
#                )
# # Unit: relative
# #           expr      min       lq     mean   median       uq       max neval cld
# # x <- B0_cpp(X) 1.000000 1.000000 1.000000 1.000000 1.000000 1.0000000   500  a
# #     x <- B0(X) 2.353846 1.853933 1.735317 1.802083 1.675926 0.4571244   500   b
# microbenchmark(x <- B2_cpp(t(X)), # Achtung: B2_ccp ist O(N^2), B2 ist O(N) -> also sollte B2 irgendwann besser sein as B2_cpp -> für N = 200, d = 30 ist das auch schon der Fall!#
#                x <- B2(X),
#                times = 100,
#                unit = "relative"
#                )
# # Unit: relative
# #           expr      min       lq     mean   median      uq        max neval cld
# # x <- B2_cpp(X) 1.000000 1.000000 1.000000 1.000000 1.00000 1.00000000   500   a
# #     x <- B2(X) 1.857963 1.647732 1.151706 1.584923 1.58509 0.01822618   500   a
# microbenchmark(x <- B3_cpp(t(X)), # hier müsste der Unterschied erhalten bleiben, da beide Funktionen O(N^3) sind
#                x <- B3(X),
#                times = 100,
#                unit = "relative"
#                )
# # Unit: relative
# #           expr      min       lq     mean   median      uq      max neval cld
# # x <- B3_cpp(X)  1.00000  1.00000  1.00000  1.00000  1.0000   1.0000   500  a
# #     x <- B3(X) 15.12839 14.66063 15.68241 14.65206 14.0193 145.9329   500   b
