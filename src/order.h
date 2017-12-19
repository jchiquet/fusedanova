#include <Rcpp.h>
#include <vector>

Rcpp::IntegerVector get_order(const int N, const Rcpp::IntegerMatrix& merge, const Rcpp::IntegerVector& node_size) ;