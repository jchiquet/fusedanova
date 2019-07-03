#ifndef _fusiontree_ward1d_H
#define _fusiontree_ward1d_H

#include "node_ward1d.h"
#include <Rcpp.h>
#include <vector> 
#include <queue>  

class FusionTree_ward1d {
  public:
    
  // the number of group to fuse
  int K ;
  
  // vector of the successive nodes in the fusion tree
  std::vector<node_ward1d> nodes ;
  
  // heap to handle the fusion events
  std::priority_queue <Fusion_ward1d, std::vector<Fusion_ward1d>, UpcomingFusions_ward1d> CandidateFusions ; 
  
  // data frame for outputing the path of fusions
  Rcpp::DataFrame tree ;
  
  // constructor
  FusionTree_ward1d(const Rcpp::NumericVector sum_0, const Rcpp::NumericVector sum2_0, const Rcpp::IntegerVector size0) ;
  
  // perform fusion between two nodes by creating a new node with label 'label'
  void fuse(int label, Fusion_ward1d fusion) ;
  
  // update the list of upcoming fusion events
  void update() ;
  
  // methods to export the path of solution
  void export_tree () ;
  
};

#endif        
