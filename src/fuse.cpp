#include <Rcpp.h>
#include "nosplit.h"

using namespace Rcpp ;

//' @export
// [[Rcpp::export]]
DataFrame fuse(NumericVector x, NumericVector slopes, NumericVector ngroup){

	Group *G = maketree(&x[0], x.length(), &slopes[0], &ngroup[0]);

  int nrow = 2*G->len -1; // nb infos = K init gpes + K-1 fusions 
  
  NumericVector beta(nrow), lambda(nrow), slope(nrow);
  IntegerVector idown(nrow), iup(nrow);
  
  int row=0;
  
  add_results(G, &beta[0], &lambda[0], &slope[0], &idown[0], &iup[0], &row);
  
  DataFrame L = DataFrame::create(Named("beta"  , beta),
                                  Named("lambda", lambda),
                                  Named("slope" , slope),
                                  Named("idown" , idown),
                                  Named("iup"   , iup));
	
	delete_tree(G); // delete the old tree for memory

	return L;
}


// we process each dimension individually using this function
// RcppExport SEXP noSplitcv(SEXP R_x,SEXP R_xv,SEXP R_ngroup, SEXP R_xtest,SEXP R_ngrouptest ,SEXP R_args) {
// 
// 	NumericVector x(R_x);
// 	NumericVector xv(R_xv);
// 	NumericVector xtest(R_xtest);
// 	NumericVector ngroup(R_ngroup);
// 	NumericVector ngrouptest(R_ngrouptest);
// 	List args(R_args);
//   
// 	std::string weights = Rcpp::as<std::string>(args["weights"]); 
// 	double gamma = Rcpp::as<double>(args["gamma"]); 
// 	NumericMatrix W = args["W"];
// 	NumericVector lambdalist = args["lambdalist"];
// 
// 	NumericVector error(lambdalist.length());
// 
// 	vector<double> sl = calculateSlope(x,ngroup,xv,weights,gamma,W,x.length());
//  
// 	Group *G = maketree(&x[0], x.length(), &sl[0],&ngroup[0]);
// 
// 	error_cv(G,&lambdalist[0],lambdalist.length(),&xtest[0], &ngrouptest[0],&error[0]);
// 
// 	delete_tree(G); 
// 
// 	return(error);
// }
// 
