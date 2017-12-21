// [[Rcpp::plugins(cpp11)]]

#include <bits/stdc++.h>
#include <stdio.h>
#include "hclust_format.h"

using namespace Rcpp ;
using namespace std ;

# define ZERO 1e-16

class node {
public:
  double lambda ; 
  double beta   ;
  double slope  ;
  int size      ;
  int label     ;
  int idown     ;
  int isplit    ;
  int iup       ;
  int down      ;
  int up        ;
  bool active   ;
  
  node() {};
  
  node(int n, int label_, double beta_, double slope_, int size_) 
    : lambda(0.0), beta(beta_), slope(slope_), size(size_),
      label(label_), idown(label_), isplit(label_), iup(label_), active(true) {
      if (label == 0  ) down = -1; else  down = label - 1;
      if (label == n-1) up   = -1; else  up   = label + 1;
    } ;
  
  bool has_down () const {return (down != -1);} ;
  bool has_up   () const {return (up   != -1);} ;

};

class Fusion {
  
private:
  node *node1 ;
  node *node2 ;
  double lambda ;

public:
  Fusion(node *node1_, node *node2_) : node1(node1_),  node2(node2_) {
    lambda = (node1->beta - node2->beta - node1->slope * node1->lambda + node2->slope * node2->lambda) / (node2->slope - node1->slope) ;
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

class FusionTree {
public:
  vector<node> nodes ;

  // constructor
  FusionTree(const NumericVector beta0, const NumericVector slope0, const IntegerVector size0) {
    int_fast32_t n = size0.size() ;
    nodes = vector<node> (2 * n - 1) ;
    for (int_fast32_t k=0; k < n; k++) {
      node node_(n, k, beta0[k], slope0[k], size0[k]) ;
      nodes.push_back(node_) ;
    }
  } ;
  
};

// node merge_nodes(int_fast32_t k, Rule rule) {
//   
// };

//' @export
// [[Rcpp::export]]
List fusedanova_cpp(NumericVector beta0, NumericVector slope0, IntegerVector size0) {

  // VARIABLES DECLARATION
  int_fast32_t n = size0.size() ;
  vector<node> nodes ; // the vector of the successive nodes in the fusion tree
  nodes.reserve(2 * n - 1) ;
  priority_queue <Fusion, vector<Fusion>, UpcomingFusions> CandidateFusions ; // a heap to handle to fusion events
  IntegerMatrix merge (n - 1, 2) ; // a matrix to encode in the hclust format the final fusion tree

  // FusionTree myTree (beta0, slope0, size0) ;

  NumericVector beta(n-1), lambda(n-1) ;
  IntegerVector idown(n-1), iup(n-1), isplit(n-1), sizes(n-1) ;
  
  // INITIALIZATION OF THE FIRST N-1 NODES (INITIAL GROUPS)
  // std::cout << "starting initialization" << std::endl;
  for (int k=0; k < n; k++) {
    node node_(n, k, beta0[k], slope0[k], size0[k]) ;
    nodes.push_back(node_) ;
  }
  // std::cout << "initialization done" << std::endl;
    
  // FIRST SET OF RULE BETWEEN THE SUCCSSIVE N-2 PAIRS OF NODES AT THE BOTTOM OF THE TREE
  for (int_fast32_t k=0; k < (n-1); k++) {
    Fusion candidate = Fusion(&nodes[k], &nodes[k+1]);
    if (candidate.get_lambda() > 0)
      CandidateFusions.push(candidate);
  }

  // n-1 FUSIONS WILL OCCUR  
  for (int_fast32_t k = n; k < (2*n-1);  k++) {
    // std::cout << "Fusion #" << k-n+1 << std::endl;
    
    // FIND THE NEXT ACTIVE RULE
    while (!CandidateFusions.top().is_active()) {
      CandidateFusions.pop() ;
    }
    Fusion candidate = CandidateFusions.top();
    CandidateFusions.pop() ;
    
    // MERGE THE TWO FUSING NODES
    // Updating the fields of the DataFrame of the fusing group
    int_fast32_t group1 = candidate.label1(), group2 = candidate.label2();
    int size_ = nodes[group1].size + nodes[group2].size ;
    double lambda_ = candidate.get_lambda();
    double slope_    = (nodes[group1].size * nodes[group1].slope + nodes[group2].size * nodes[group2].slope) / size_ ;
    double beta_     = nodes[group1].beta + (lambda_ - nodes[group1].lambda) * nodes[group1].slope ;
    node node_(n, k, beta_, slope_, size_) ;
    nodes.push_back(node_) ;
    
    nodes[k].lambda = lambda_ ;
    if(nodes[group1].idown < nodes[group2].idown) {
      nodes[k].isplit      = nodes[group1].iup      ;
      nodes[k].idown        = nodes[group1].idown       ;
      nodes[k].iup       = nodes[group2].iup      ;
      nodes[k].down      = nodes[group1].down     ;
      nodes[k].up     = nodes[group2].up    ;
    } else {
      nodes[k].isplit      = nodes[group2].iup      ; 
      nodes[k].idown        = nodes[group2].idown       ;
      nodes[k].iup       = nodes[group1].iup      ;
      nodes[k].down      = nodes[group2].down     ;
      nodes[k].up     = nodes[group1].up    ;
    }
    
    // Updating neighbors
    if(nodes[k].has_down() )  nodes[nodes[k].down ].up   = k;
    if(nodes[k].has_up()   )  nodes[nodes[k].up   ].down = k;
    
    nodes[group1].active = false ;
    nodes[group2].active = false ;
    nodes[k].active = true ;
    
    // get new rules and add them to the queue
    if (nodes[k].has_down() & nodes[nodes[k].down].active) {
      Fusion candidate= Fusion(&nodes[nodes[k].down], &nodes[k]);
      if (candidate.get_lambda() > 0)
        CandidateFusions.push(candidate);
    }
    
    if (nodes[k].has_up() & nodes[nodes[k].up].active) {
      Fusion candidate= Fusion(&nodes[k], &nodes[nodes[k].up]);
      if (candidate.get_lambda() > 0)
        CandidateFusions.push(candidate);
    }

    // outputing in hclust merge format
    merge(k-n,_) = hc_merge(n, group1, group2) ;
      
    // outputing the path 
    beta   (k-n) = nodes[k].beta ;
    lambda (k-n) = nodes[k].lambda ;
    idown  (k-n) = nodes[k].idown + 1;
    iup    (k-n) = nodes[k].iup + 1;
    isplit (k-n) = nodes[k].isplit + 1;
    sizes (k-n) = nodes[k].size;
  }
  
  DataFrame path = DataFrame::create(
        Named("beta"  ) = beta  ,
        Named("lambda") = lambda,
        Named("down"  ) = idown ,
        Named("high"  ) = iup   ,
        Named("split" ) = isplit
  );
  
  IntegerVector order = hc_order(n, merge, sizes);
  
  return List::create( Named("path") = path, Named("merge") = merge, Named("order") = order);
  
}
