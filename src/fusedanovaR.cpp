// [[Rcpp::plugins(cpp11)]]
#include "fusedanova.h"
#include "hclust_format.h"

using namespace Rcpp;
using namespace std;

DataFrame format_path(const vector<node> nodes) {
  
  // Variable to ouput the fusion tree to R
  int K = nodes[0].range ;
  int nfusion = nodes.size()-K ;
  NumericVector beta(nfusion), lambda(nfusion), slope(nfusion)  ;
  IntegerVector idown(nfusion), iup(nfusion), isplit(nfusion), sizes(nfusion) ;
  
  for (int k = K; k < nodes.size(); k++) {
    beta   (k-K) = nodes[k].beta       ;
    lambda (k-K) = nodes[k].lambda     ;
    idown  (k-K) = nodes[k].idown  + 1 ;
    iup    (k-K) = nodes[k].iup    + 1 ;
    isplit (k-K) = nodes[k].isplit + 1 ;
    sizes  (k-K) = nodes[k].size       ;
    slope  (k-K) = nodes[k].slope      ;
  }
  
  return(DataFrame::create(
      Named("beta"  ) = beta  ,
      Named("lambda") = lambda,
      Named("down"  ) = idown ,
      Named("up"    ) = iup   ,
      Named("split" ) = isplit,
      Named("sizes" ) = sizes,
      Named("slopes" ) = slope
  ));
}

//' @export
// [[Rcpp::export]]
List fusedanova_cpp(NumericVector beta0, NumericVector slope0, IntegerVector size0) {

  // Problem Dimension
  int n = size0.size() ;
  
  // VARIABLES DECLARATION

  // FusionTree myTree (beta0, slope0, size0) ; // TODO
  
  // the vector of the successive nodes in the fusion tree
  vector<node> nodes ;  nodes.reserve(2 * n - 1) ;
  
  // a heap to handle to fusion events
  priority_queue <Fusion, vector<Fusion>, UpcomingFusions> CandidateFusions ; 
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
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
    // std::cout << "Fusion #" << k-n+1 << std::endl;
    
    while (!CandidateFusions.top().is_active() & !CandidateFusions.empty()) {
      std::cout << CandidateFusions.top().get_lambda() << std::endl;
      CandidateFusions.pop() ;
      R_CheckUserInterrupt() ;
    }
    if (CandidateFusions.empty()) {
      std::cout << "ouch: no more active fusions: you obviously chose a too large gamma" << std::endl;
      nodes.resize(k);
      break;
    }

    Fusion fusion = CandidateFusions.top();
    std::cout << "fusion at" << fusion.get_lambda() << std::endl;
    CandidateFusions.pop() ;

    // MERGE THE TWO FUSING NODES
    node node_ = *fusion.node1 + *fusion.node2;
    node_.label = k;
    fusion.node1->active = false;
    fusion.node2->active = false;
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
  } // end fusion loop

  // outputing the path 
  DataFrame path = format_path(nodes) ;

  // outputing in hclust merge format
  IntegerMatrix merge = hc_merge(nodes) ;
  
  print(merge);
  
  // recovering order for plotting dendrogram  
  IntegerVector order = hc_order(merge, path["sizes"]);

  print(order) ;
  // Send back everything
  return List::create( Named("path") = path, Named("merge") = merge, Named("order") = order);
  
}
