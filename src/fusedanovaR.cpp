// [[Rcpp::plugins(cpp11)]]
#include "fusedanova.h"
#include "hclust_format.h"

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List fusedanova_cpp(NumericVector beta0, NumericVector slope0, IntegerVector size0) {

  // Problem Dimension
  int n = size0.size() ;
  
  // VARIABLES DECLARATION

  // FusionTree myTree (beta0, slope0, size0) ;
  
  // the vector of the successive nodes in the fusion tree
  vector<node> nodes ;  nodes.reserve(2 * n - 1) ;
  
  // a heap to handle to fusion events
  priority_queue <Fusion, vector<Fusion>, UpcomingFusions> CandidateFusions ; 
  
  // A matrix to encode in the hclust format the final fusion tree
  IntegerMatrix merge (n - 1, 2) ; 

  // Variable to ouput the fusion tree to R
  NumericVector beta(n-1), lambda(n-1) ;
  IntegerVector idown(n-1), iup(n-1), isplit(n-1), sizes(n-1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  // std::cout << "starting initialization" << std::endl;
  for (int k=0; k < n; k++) {
    node node_(n, k, beta0[k], slope0[k], size0[k]) ;
    nodes.push_back(node_) ;
  }

  // FIRST SET OF RULES BETWEEN THE SUCCESSIVE N-2 PAIRS OF NODES AT THE BOTTOM OF THE TREE
  for (int k=0; k < (n-1); k++) {
    Fusion candidate = Fusion(&nodes[k], &nodes[k+1]);
    if (candidate.get_lambda() >= 0)
      CandidateFusions.push(candidate);
  }

  // n-1 FUSIONS _MUST_ OCCUR
  for (int k = n; k < (2*n-1);  k++) {
    std::cout << "Fusion #" << k-n+1 << std::endl;
    
    // FIND THE NEXT ACTIVE RULE
    while (!CandidateFusions.top().is_active()) {
      CandidateFusions.pop() ;
      R_CheckUserInterrupt() ;
    }

    // while (!CandidateFusions.top().is_active() & !CandidateFusions.empty()) {
    //   CandidateFusions.pop() ;
    //   R_CheckUserInterrupt() ;
    // }
    // if (CandidateFusions.empty()) {
    //   std::cout << "ouch: no more active fusions" << std::endl;
    //   break;
    // }

    Fusion candidate = CandidateFusions.top();
    CandidateFusions.pop() ;

    // MERGE THE TWO FUSING NODES
    node node_ = *candidate.node1 + *candidate.node2;
    node_.label = k;
    candidate.node1->active = false;
    candidate.node2->active = false;
    nodes.push_back(node_) ;

    // Update neighbors and check if some new rules must be added to the queue
    
    // down neighbor...
    if (nodes.back().has_down()) {
      nodes[nodes.back().down ].up   =  nodes.size() - 1 ;
      if (nodes[nodes.back().down].active) {
        Fusion candidate_ = Fusion(&nodes[nodes.back().down], &nodes.back());
        if (candidate_.get_lambda() >= 0)
          CandidateFusions.push(candidate_);
      }
    }
    
    // up neighbor...
    if (nodes.back().has_up()) {
      nodes[nodes.back().up   ].down =  nodes.size() - 1 ;
      if (nodes[nodes.back().up].active) {
        Fusion candidate_ = Fusion(&nodes.back(), &nodes[nodes.back().up]);
        if (candidate_.get_lambda() >= 0)
          CandidateFusions.push(candidate_);
      }
    }

    // outputing in hclust merge format
    merge(k-n,_) = hc_merge(n, candidate.label1(), candidate.label2()) ;
      
    // outputing the path 
    beta   (k-n) = nodes.back().beta       ;
    lambda (k-n) = nodes.back().lambda     ;
    idown  (k-n) = nodes.back().idown + 1  ;
    iup    (k-n) = nodes.back().iup + 1    ;
    isplit (k-n) = nodes.back().isplit + 1 ;
    sizes  (k-n) = nodes.back().nsize      ;

  } 
    
  DataFrame path = DataFrame::create(
        Named("beta"  ) = beta  ,
        Named("lambda") = lambda,
        Named("down"  ) = idown ,
        Named("high"  ) = iup   ,
        Named("split" ) = isplit
  );

  print(sizes);
  IntegerVector order = hc_order(n, merge, sizes);

  return List::create( Named("path") = path, Named("merge") = merge, Named("order") = order);
  
}
