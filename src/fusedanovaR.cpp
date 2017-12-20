#include <bits/stdc++.h>
#include <stdio.h>
#include "hclust_format.h"

using namespace Rcpp ;
using namespace std ;

# define ZERO 1e-16

struct node {
  double lambda         ; 
  double beta           ;
  double slope          ;
  int_fast32_t  label    ;
  int_fast32_t i_low    ;
  int_fast32_t i_split  ;
  int_fast32_t i_high   ;
  int_fast32_t grp_low  ;
  int_fast32_t grp_high ;
  int_fast32_t grp_size ;
  bool active           ;
  bool has_grp_low      ;
  bool has_grp_high     ;
};

class Fusion {
  
private:
  node *node1 ;
  node *node2 ;
  double lambda ;
  void set_lambda() {
    lambda = (node1->beta - node2->beta - node1->slope * node1->lambda + node2->slope * node2->lambda) / (node2->slope - node1->slope) ;
  } ;
  
public:
  Fusion(node *node1_, node *node2_) : node1(node1_),  node2(node2_) {
    set_lambda();
  } ;
  int_fast32_t label1() {return node1->label ;}  
  int_fast32_t label2() {return node2->label ;}  
  double get_lambda() const {return lambda ;}
  bool is_active() const {return(node1->active & node2->active);}
};

class UpcomingFusions {
public:
  double operator() (const Fusion& r1, const Fusion& r2) {
    return r1.get_lambda() > r2.get_lambda();
  }
};

// node merge_nodes(int_fast32_t k, Rule rule) {
//   
// };

//' @export
// [[Rcpp::export]]
List fuse(NumericVector beta0, NumericVector slope0, IntegerVector grp_size0) {

  // VARIABLES DECLARATION
  int_fast32_t n = grp_size0.size() ;
  vector<node> nodes  (2*n - 1)  ; // the vector of the successive nodes in the fusion tree
  priority_queue <Fusion, vector<Fusion>, UpcomingFusions> CandidateFusions ; // a heap to handle to fusion events
  IntegerMatrix merge (n - 1, 2) ; // a matrix to encode in the hclust format the final fusion tree

  NumericVector beta(n-1), lambda(n-1) ;
  IntegerVector down(n-1), high(n-1), split(n-1), sizes(n-1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  for (int_fast32_t k=0; k < n; k++) {
    nodes[k].lambda   = 0.0;
    nodes[k].beta     = beta0[k];
    nodes[k].slope    = slope0[k];
    nodes[k].label    = k;
    nodes[k].i_low    = k;
    nodes[k].i_split  = k;
    nodes[k].i_high   = k;
    nodes[k].grp_size = grp_size0[k] ;
    nodes[k].active   = true;
    if (k == 0) {
      nodes[k].has_grp_low = false ;
      nodes[k].grp_low = -1;
    } else {
      nodes[k].has_grp_low = true ;
      nodes[k].grp_low = k - 1;
    }
    if (k == (n-1)) {
      nodes[k].has_grp_high = false ;
      nodes[k].grp_high = -1 ;
    } else {
      nodes[k].has_grp_high = true ;
      nodes[k].grp_high = k + 1;
    }
  }
  
  // FIRST SET OF RULE BETWEEN THE SUCCSSIVE N-2 PAIRS OF NODES AT THE BOTTOM OF THE TREE
  for (int_fast32_t k=0; k < (n-1); k++) {
    Fusion candidate= Fusion(&nodes[k], &nodes[k+1]);
    if (candidate.get_lambda() > 0)
      CandidateFusions.push(candidate);
  }

  // n-1 FUSIONS WILL OCCUR  
  for (int_fast32_t k = n; k < (2*n-1);  k++) {
    
    // FIND THE NEXT ACTIVE RULE
    while (!CandidateFusions.top().is_active()) {
      CandidateFusions.pop() ;
    }
    Fusion candidate = CandidateFusions.top();
    CandidateFusions.pop() ;
    
    // MERGE THE TWO FUSING NODES
    // Updating the fields of the DataFrame of the fusing group
    int_fast32_t group1 = candidate.label1(), group2 = candidate.label2();
    double lambda_k = candidate.get_lambda();
    
    // merge(k, rule_.get_lambda(), rule_.get_group1(), rule_.get_group2(), table) ;
    
    nodes[k].label  = k ;
    nodes[k].lambda = lambda_k ;
    nodes[k].grp_size = nodes[group1].grp_size + nodes[group2].grp_size ;
    nodes[k].slope    = (nodes[group1].grp_size * nodes[group1].slope + nodes[group2].grp_size * nodes[group2].slope) / nodes[k].grp_size ;
    nodes[k].beta     = nodes[group1].beta + (lambda_k - nodes[group1].lambda) * nodes[group1].slope ;
    if(nodes[group1].i_low < nodes[group2].i_low) {
      nodes[k].i_split      = nodes[group1].i_high      ;
      nodes[k].i_low        = nodes[group1].i_low       ;
      nodes[k].i_high       = nodes[group2].i_high      ;
      nodes[k].grp_low      = nodes[group1].grp_low     ;
      nodes[k].grp_high     = nodes[group2].grp_high    ;
      nodes[k].has_grp_low  = nodes[group1].has_grp_low ;
      nodes[k].has_grp_high = nodes[group2].has_grp_high;
    } else {
      nodes[k].i_split      = nodes[group2].i_high      ; 
      nodes[k].i_low        = nodes[group2].i_low       ;
      nodes[k].i_high       = nodes[group1].i_high      ;
      nodes[k].grp_low      = nodes[group2].grp_low     ;
      nodes[k].grp_high     = nodes[group1].grp_high    ;
      nodes[k].has_grp_low  = nodes[group2].has_grp_low ;
      nodes[k].has_grp_high = nodes[group1].has_grp_high;
    }
    
    // Updating neighbors
    if(nodes[k].has_grp_low )  nodes[nodes[k].grp_low ].grp_high = k;
    if(nodes[k].has_grp_high)  nodes[nodes[k].grp_high].grp_low  = k;
    
    nodes[group1].active = false ;
    nodes[group2].active = false ;
    nodes[k].active = true ;
    
    // get new rules and add them to the queue
    if (nodes[k].has_grp_low & nodes[nodes[k].grp_low].active) {
      Fusion candidate= Fusion(&nodes[nodes[k].grp_low], &nodes[k]);
      if (candidate.get_lambda() > 0)
        CandidateFusions.push(candidate);
    }
    
    if (nodes[k].has_grp_high & nodes[nodes[k].grp_high].active) {
      Fusion candidate= Fusion(&nodes[k], &nodes[nodes[k].grp_high]);
      if (candidate.get_lambda() > 0)
        CandidateFusions.push(candidate);
    }

    // outputing in hclust merge format
    merge(k-n,_) = hc_merge(n, group1, group2) ;
      
    // outputing the path 
    beta  (k-n) = nodes[k].beta ;
    lambda(k-n) = nodes[k].lambda ;
    down  (k-n) = nodes[k].i_low + 1;
    high  (k-n) = nodes[k].i_high + 1;
    split (k-n) = nodes[k].i_split + 1;
    sizes (k-n) = nodes[k].grp_size;
  }
  
  DataFrame path = DataFrame::create(
        Named("beta"  ) = beta,
        Named("lambda") = lambda,
        Named("down"  ) = down,
        Named("high"  ) = high,
        Named("split" ) = split
  );
  
  IntegerVector order = hc_order(n, merge, sizes);
  
  return List::create( Named("path") = path, Named("merge") = merge, Named("order") = order);
  
}
