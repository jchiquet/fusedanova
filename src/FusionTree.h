#ifndef _fusiontree_H
#define _fusiontree_H

#include "node.h"
#include <Rcpp.h>
#include <vector> 
#include <queue>  

class FusionTree {
  public:
    
  // the number of group to fuse
  int K ;
  
  // vector of the successive nodes in the fusion tree
  std::vector<node> nodes ;
  
  // heap to handle the fusion events
  std::priority_queue <Fusion, std::vector<Fusion>, UpcomingFusions> CandidateFusions ; 
  
  // data frame for outputing the path of fusions
  Rcpp::DataFrame tree ;
  
  // constructor
  FusionTree(const Rcpp::NumericVector beta0, const Rcpp::NumericVector slope0, const Rcpp::IntegerVector size0) ;
  
  // perform fusion between two nodes by creating a new node with label 'label'
  void fuse(int label, Fusion fusion) ;
  
  // update the list of upcoming fusion events
  void update() ;
  
  // methods to export the path of solution
  void export_tree () ;
  
};

#endif        
