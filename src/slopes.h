#include <Rcpp.h>
#include <string>
#include <vector>

# define square(x) ((x)*(x))

Rcpp::NumericVector get_slopes( Rcpp::NumericVector &xm    ,
                                Rcpp::IntegerVector &ngroup,
                                Rcpp::NumericVector &xv    , 
                                std::string weights        ,
                                double gamma               ,
                                Rcpp::NumericMatrix &W     ) ;