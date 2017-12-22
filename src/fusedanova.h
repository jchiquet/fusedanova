#ifndef _fusedanova_H
#define _fusedanova_H

#include <Rcpp.h>
#include <bits/stdc++.h>
#include <stdio.h>

# define ZERO 1e-16

class node {
public:
  double lambda ; 
  double beta   ;
  double slope  ;
  int range     ;
  int size      ;
  int nsize     ;
  int label     ;
  int idown     ;
  int isplit    ;
  int iup       ;
  int down      ;
  int up        ;
  bool active   ;

  // Constructors/Destructor
   node() ;
  ~node() ;
  node(int n, int label_, double beta_, double slope_, int size_) ;

  // Basic methods for acces
  bool has_down () const {return (down != -1);} ;
  bool has_up   () const {return (up   != -1);} ;

  // Define overloaded + for fusing two nodes
  node operator+ (const node& node_) ;
};

class Fusion {

public:
  node *node1 ;
  node *node2 ;
  double lambda ;

    // Constructor
  Fusion(node *node1_, node *node2_) ; 
  // Getter
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

// class FusionTree {
// public:
//   
//   std::vector<node> nodes ;
// 
//   // constructor
//   FusionTree(const Rcpp::NumericVector beta0, const Rcpp::NumericVector slope0, const Rcpp::IntegerVector size0) {} ;
// 
// };

#endif        
