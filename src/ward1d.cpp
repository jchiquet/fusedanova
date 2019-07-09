// [[Rcpp::plugins(cpp11)]]
#include "FusionTree_ward1d.h"

#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame ward1d_cpp(NumericVector sum_0, NumericVector sum2_0, IntegerVector size0) {

  // VARIABLES DECLARATION

  FusionTree_ward1d myTree (sum_0, sum2_0, size0) ;

  // n-1 FUSIONS _MUST_ OCCUR
  for (int k = myTree.K; k < (2*myTree.K - 1);  k++) {
    // std::cout << "Fusion #" << k-n+1 << std::endl;
    
    // Find the next best active fusion
    while (!myTree.CandidateFusions.top().is_active() & !myTree.CandidateFusions.empty()) {
      myTree.CandidateFusions.pop() ;
    }
      
    if (myTree.CandidateFusions.empty()) {
      std::cout << "ouch!! no more active fusions..." << std::endl;
      myTree.nodes.resize(k);
      break;
    }
    Fusion_ward1d fusion = myTree.CandidateFusions.top();
      
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
