// [[Rcpp::plugins(cpp11)]]
#include "FusionTree.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame fusedanova_cpp(NumericVector beta0, NumericVector slope0, IntegerVector size0) {

  // VARIABLES DECLARATION

  FusionTree myTree (beta0, slope0, size0) ;

  // n-1 FUSIONS _MUST_ OCCUR
  for (int k = myTree.K; k < (2*myTree.K - 1);  k++) {
    // std::cout << "Fusion #" << k-n+1 << std::endl;
    
    // Find the next best active fusion
    while (!myTree.CandidateFusions.top().is_active() & !myTree.CandidateFusions.empty()) {
      myTree.CandidateFusions.pop() ;
    }
      
    if (myTree.CandidateFusions.empty()) {
      Rcpp::Rcout << "ouch!! no more active fusions: you obviously chose a too large gamma" << std::endl;
      myTree.nodes.resize(k);
      break;
    }
    Fusion fusion = myTree.CandidateFusions.top();
    myTree.CandidateFusions.pop() ;

    // Merge the two fusiong nodes
    myTree.fuse(k, fusion) ;

    // Update neighbors and check if some new rules must be added to the queue
    myTree.update();
  } // end fusion loop

  // outputing the path 
  myTree.export_tree() ;

  return myTree.tree ;
}
