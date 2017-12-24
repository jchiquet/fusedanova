#ifndef _hclust_format_H
#define _hclust_format_H

#include <Rcpp.h>
#include <vector>
#include "fusedanova.h"

Rcpp::IntegerVector hc_order(const Rcpp::IntegerMatrix&, const Rcpp::IntegerVector&) ;
Rcpp::IntegerMatrix hc_merge(const std::vector<node>) ;
  
#endif
