// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <bits/stdc++.h>
#include <stdio.h>

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
  std::vector<node> nodes ;

  // constructor
  FusionTree(const Rcpp::NumericVector beta0, const Rcpp::NumericVector slope0, const Rcpp::IntegerVector size0) {
    int_fast32_t n = size0.size() ;
    nodes = std::vector<node> (2 * n - 1) ;
    for (int_fast32_t k=0; k < n; k++) {
      node node_(n, k, beta0[k], slope0[k], size0[k]) ;
      nodes.push_back(node_) ;
    }
  } ;
  
};

// node merge_nodes(int_fast32_t k, Rule rule) {
//   
// };
