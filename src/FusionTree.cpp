#include "FusionTree.h"

using namespace Rcpp ;

FusionTree::FusionTree(const Rcpp::NumericVector beta0, const Rcpp::NumericVector slope0, const Rcpp::IntegerVector size0) {
  
  // Problem Dimension
  K = size0.size() ;
  nodes.reserve(2 * K - 1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  for (int k=0; k < K; k++) {
    node node_(K, k, beta0[k], slope0[k], size0[k]) ;
    nodes.push_back(node_) ;
  }
  
  // FIRST SET OF RULES BETWEEN THE SUCCESSIVE N-2 PAIRS OF NODES AT THE BOTTOM OF THE TREE
  for (int k=0; k < (K-1); k++) {
    Fusion candidate = Fusion(&nodes[k], &nodes[k+1]);
    if (candidate.get_lambda() >= 0)
      CandidateFusions.push(candidate);
  }
  
};

void FusionTree::fuse(int label, Fusion fusion) {
  // MERGE THE TWO FUSING NODES
  node node_ = *fusion.node1 + *fusion.node2;
  node_.label = label;
  fusion.node1->active = false;
  fusion.node2->active = false;
  nodes.push_back(node_) ;
};

void FusionTree::update() {
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
};

void FusionTree::export_path() {
  
  // function to ouput the fusion tree to R
  int nfusion = nodes.size() - K ;
  Rcpp::NumericVector beta(nfusion), lambda(nfusion), slope(nfusion)  ;
  Rcpp::IntegerVector idown(nfusion), iup(nfusion), isplit(nfusion), sizes(nfusion) ;
  
  for (int k = K; k < nodes.size(); k++) {
    beta   (k-K) = nodes[k].beta       ;
    lambda (k-K) = nodes[k].lambda     ;
    idown  (k-K) = nodes[k].idown  + 1 ;
    iup    (k-K) = nodes[k].iup    + 1 ;
    isplit (k-K) = nodes[k].isplit + 1 ;
    sizes  (k-K) = nodes[k].size       ;
    slope  (k-K) = nodes[k].slope      ;
  }
  
  path = DataFrame::create(
    Rcpp::Named("beta"  ) = beta  ,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("down"  ) = idown ,
    Rcpp::Named("up"    ) = iup   ,
    Rcpp::Named("split" ) = isplit,
    Rcpp::Named("sizes" ) = sizes ,
    Rcpp::Named("slopes") = slope
  );
};

void FusionTree::export_merge() {
  
  // function to ouput the dendrogram of fusion into the hclust format 
  merge = IntegerMatrix(nodes.size() - K, 2);
  
  for (int k = K; k < nodes.size(); k++) {
    
    int_fast32_t merge1, merge2;
    
    if (nodes[k].parent1 >= K) {
      merge1 = (nodes[k].parent1 + 1) - K;
    } else {
      merge1 = -(nodes[k].parent1 + 1);
    }
    if (nodes[k].parent2 >= K) {
      merge2 = (nodes[k].parent2 + 1) - K;
    } else {
      merge2 = -(nodes[k].parent2 + 1);
    }
    if (merge1 < merge2) {
      merge(k-K, 0) = merge1;
      merge(k-K, 1) = merge2;
    }
    else {
      merge(k-K, 0) = merge2;
      merge(k-K, 1) = merge1;
    }
  }
};


struct pos_node {
  int_fast32_t pos;
  int_fast32_t node;
};

void FusionTree::export_order() {
  // function to recover a correct order of the nodes for the dendrogram (required for hclust object) 
  
  /* Parameters:
    N         : number of data points
  merge     : (N-1)×2 array which specifies the node indices which are merged in each step of the clustering procedure.
  Negative entries -1...-N point to singleton nodes, while positive entries 1...(N-1) point to nodes which are themselves parents of other nodes.
  node_size : array of node sizes - makes it easier
  order     : output array of size N
  Runtime: Θ(N)
  */
    
    IntegerVector node_size = path["sizes"];
  int N = node_size.size() + 1 ;
  std::vector<pos_node> queue(N/2);
  
  int_fast32_t parent;
  int_fast32_t child;
  int_fast32_t pos = 0;
  order = IntegerVector (N);
  queue[0].pos = 0;
  queue[0].node = N-2;
  int_fast32_t idx = 1;
  
  do {
    --idx;
    pos = queue[idx].pos;
    parent = queue[idx].node;
    
    // First child
    child = merge(parent, 0);
    if (child<0) { // singleton node, write this into the 'order' array.
      order[pos] = -child;
      ++pos;
    }
    else { /* compound node: put it on top of the queue and decompose it in a later iteration. */
        queue[idx].pos = pos;
        queue[idx].node = child-1; // convert index-1 based to index-0 based
        ++idx;
        pos += node_size[child-1];
    }
    // Second child
    child = merge(parent,1);
    if (child<0) {
      order[pos] = -child;
    }
    else {
      queue[idx].pos = pos;
      queue[idx].node = child-1;
      ++idx;
    }
  } while (idx>0);
  
};
