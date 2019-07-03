#include "FusionTree_ward1d.h"

using namespace Rcpp ;

FusionTree_ward1d::FusionTree_ward1d(const Rcpp::NumericVector sum_0, const Rcpp::NumericVector sum2_0, const Rcpp::IntegerVector size0) {
  
  // Problem Dimension
  K = size0.size() ;
  nodes.reserve(2 * K - 1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  for (int k=0; k < K; k++) {
    node_ward1d node_(K, k, sum_0[k], sum2_0[k], size0[k]) ;
    nodes.push_back(node_) ;
  }
  
  // FIRST SET OF RULES BETWEEN THE SUCCESSIVE N-2 PAIRS OF NODES AT THE BOTTOM OF THE TREE
  for (int k=0; k < (K-1); k++) {
    Fusion_ward1d candidate = Fusion_ward1d(&nodes[k], &nodes[k+1]);
    if (candidate.get_lambda() >= 0)
      CandidateFusions.push(candidate);
  }
  
};

void FusionTree_ward1d::fuse(int label, Fusion_ward1d fusion) {
  // MERGE THE TWO FUSING NODES
  node_ward1d node_ = *fusion.node1 + *fusion.node2;
  node_.label = label;
  fusion.node1->active = false;
  fusion.node2->active = false;
  nodes.push_back(node_) ;
};

void FusionTree_ward1d::update() {
  // down neighbor...
  if (nodes.back().has_down()) {
    nodes[nodes.back().down ].up   =  nodes.size() - 1 ;
    if (nodes[nodes.back().down].active) {
      Fusion_ward1d candidate_ = Fusion_ward1d(&nodes[nodes.back().down], &nodes.back());
      if (candidate_.get_lambda() >= 0)
        CandidateFusions.push(candidate_);
    }
  }
  
  // up neighbor...
  if (nodes.back().has_up()) {
    nodes[nodes.back().up   ].down =  nodes.size() - 1 ;
    if (nodes[nodes.back().up].active) {
      Fusion_ward1d candidate_ = Fusion_ward1d(&nodes.back(), &nodes[nodes.back().up]);
      if (candidate_.get_lambda() >= 0)
        CandidateFusions.push(candidate_);
    }
  }
};

void FusionTree_ward1d::export_tree() {
  
  // function to ouput the fusion tree to R
  int nfusion = nodes.size() - K ;
  Rcpp::NumericVector beta(nfusion), lambda(nfusion), slope(nfusion)  ;
  Rcpp::IntegerVector idown(nfusion), iup(nfusion), isplit(nfusion), sizes(nfusion) ;
  Rcpp::IntegerVector label(nfusion), child1(nfusion), child2(nfusion), down(nfusion), up(nfusion) ;
  
  for (int k = K; k < nodes.size(); k++) {
    beta   (k-K) = nodes[k].xplus       ;
    lambda (k-K) = nodes[k].lambda      ;
    idown  (k-K) = nodes[k].idown   + 1 ;
    iup    (k-K) = nodes[k].iup     + 1 ;
    isplit (k-K) = nodes[k].isplit  + 1 ;
    label  (k-K) = nodes[k].label   + 1 ;
    child1(k-K)  = nodes[k].parent1 + 1 ;
    child2(k-K)  = nodes[k].parent2 + 1 ;
    sizes  (k-K) = nodes[k].size        ;
  }

  tree = DataFrame::create(
    Rcpp::Named("beta"  )       = beta   ,
    Rcpp::Named("lambda")       = lambda ,
    Rcpp::Named("index_down"  ) = idown  ,
    Rcpp::Named("index_split" ) = isplit ,
    Rcpp::Named("index_up"    ) = iup    ,
    Rcpp::Named("child1")       = child1,
    Rcpp::Named("child2")       = child2,
    Rcpp::Named("label" )       = label  ,
    Rcpp::Named("sizes" )       = sizes  ,
    Rcpp::Named("slopes")       = slope
  );
};
