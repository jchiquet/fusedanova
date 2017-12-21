#include "hclust_format.h"

using namespace Rcpp ;

struct pos_node {
  int_fast32_t pos;
  int_fast32_t node;
};

IntegerVector hc_order(const int_fast32_t N, const IntegerMatrix& merge, const IntegerVector& node_size) {
  /* Parameters:
  N         : number of data points
  merge     : (N-1)×2 array which specifies the node indices which are merged in each step of the clustering procedure.
  Negative entries -1...-N point to singleton nodes, while positive entries 1...(N-1) point to nodes which are themselves parents of other nodes.
  node_size : array of node sizes - makes it easier
  order     : output array of size N
  Runtime: Θ(N)
  */
  std::vector<pos_node> queue(N/2);
  
  int_fast32_t parent;
  int_fast32_t child;
  int_fast32_t pos = 0;
  IntegerVector order (N);
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
  
  return(order) ;
}

IntegerVector hc_merge(const int_fast32_t n, const int_fast32_t group1, const int_fast32_t group2) {
  /* 
   * encoding the output to the merge format found in hclust
   */
  int_fast32_t merge1, merge2;
  IntegerVector merge(2);
  
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
    merge(0) = merge1;
    merge(1) = merge2;
  }
  else {
    merge(0) = merge2;
    merge(1) = merge1;
  }
  return(merge);  
}