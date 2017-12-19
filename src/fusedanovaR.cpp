#include <bits/stdc++.h>
#include <stdio.h>
#include "order.h"

using namespace Rcpp ;
using namespace std ;

# define ZERO 1e-16

struct node {
  double lambda         ; 
  double beta           ;
  double slope          ;
  int label    ;
  int i_low    ;
  int i_split  ;
  int i_high   ;
  int grp_low  ;
  int grp_high ;
  int grp_size ;
  bool active           ;
  bool has_grp_low      ;
  bool has_grp_high     ;
};

double compute_lambda(const node& node1, const node& node2) {
  return((node1.beta - node2.beta - node1.slope * node1.lambda + node2.slope * node2.lambda) / (node2.slope - node1.slope));
}

class Rule {
  
private:
  node *node1 ;
  node *node2 ;
  double lambda ;
  
public:
  Rule(node *node1_, node *node2_) : node1(node1_),  node2(node2_) {
    lambda = compute_lambda(*node1_, *node2_);
  } ;
  int label1() {return node1->label ;}  
  int label2() {return node2->label ;}  
  double get_lambda() const {return lambda ;}
  bool is_active() const {return(node1->active & node2->active);}
};

class RuleComparator {
public:
  double operator() (const Rule& r1, const Rule& r2) {
    return r1.get_lambda() > r2.get_lambda();
  }
};

// node merge_nodes(int k, Rule rule) {
//   
// };

//' @export
// [[Rcpp::export]]
List fuse(NumericVector beta0, NumericVector slope0, IntegerVector grp_size0) {

  // VARIABLES DECLARATION
  int n = grp_size0.size() ;
  vector<node> nodes  (2*n - 1)  ; // the vector of the successive nodes in the fusion tree
  priority_queue <Rule, vector<Rule>, RuleComparator> myMinHeap ; // a heap to handle to fusion events
  IntegerMatrix merge (n - 1, 2) ; // a matrix to encode in the hclust format the final fusion tree

  NumericVector beta(n-1), lambda(n-1) ;
  IntegerVector down(n-1), high(n-1), split(n-1), sizes(n-1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  for (int k=0; k < n; k++) {
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
  for (int k=0; k < (n-1); k++) {
    if (compute_lambda(nodes[k], nodes[k+1]) > 0)
      myMinHeap.push(Rule(&nodes[k], &nodes[k+1]));
  }

  // n-1 kS WILL OCCUR  
  for (int k = n; k < (2*n-1);  k++) {
    
    // FIND THE NEXT ACTIVE RULE
    while (!myMinHeap.top().is_active()) {
      myMinHeap.pop() ;
    }
    Rule rule_ = myMinHeap.top();
    myMinHeap.pop() ;
    
    // MERGE THE TWO FUSING NODES
    // Updating the fields of the DataFrame of the fusing group
    int group1 = rule_.label1(), group2 = rule_.label2();
    double lambda_k = rule_.get_lambda();
    
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
      if (compute_lambda(nodes[nodes[k].grp_low], nodes[k]) > 0)
        myMinHeap.push(Rule(&nodes[nodes[k].grp_low], &nodes[k]));
    }
    
    if (nodes[k].has_grp_high & nodes[nodes[k].grp_high].active) {
      if (compute_lambda(nodes[k], nodes[nodes[k].grp_high]) > 0)
        myMinHeap.push(Rule(&nodes[k], &nodes[nodes[k].grp_high]));
    }
    
    // outputing in hclust format
    signed int merge1, merge2;
    if (group1 >= n) {
      merge1 = (group1+1)-n;
    } else {
      merge1 = -(group1+1);
    }
    if (group2 >= n) {
      merge2 = (group2+1)-n;
    } else {
      merge2 = -(group2+1);
    }
    if (merge1 < merge2) {
      merge(k-n,0) = merge1;
      merge(k-n,1) = merge2;
    }
    else {
      merge(k-n,0) = merge2;
      merge(k-n,1) = merge1;
    }

    // output the path 
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
  
  IntegerVector order = get_order(n, merge, sizes);
  
  return List::create( Named("path") = path, Named("merge") = merge, Named("order") = order);
  
}
