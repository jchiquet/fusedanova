#ifndef _hclust_format_H
#define _hclust_format_H

#include <Rcpp.h>
#include <vector>

Rcpp::IntegerVector hc_order(const int_fast32_t, const Rcpp::IntegerMatrix&, const Rcpp::IntegerVector&) ;
Rcpp::IntegerVector hc_merge(const int_fast32_t, const int_fast32_t, const int_fast32_t) ;
  
#endif
  