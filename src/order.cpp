#include "order.h"

using namespace Rcpp ;

struct pos_node {
  int pos;
  int node;
};

IntegerVector get_order(const int N, const IntegerMatrix& merge, const IntegerVector& node_size) {
  /* Parameters:
  N         : number of data points
  merge     : (N-1)×2 array which specifies the node indices which are merged in each step of the clustering procedure.
  Negative entries -1...-N point to singleton nodes, while positive entries 1...(N-1) point to nodes which are themselves parents of other nodes.
  node_size : array of node sizes - makes it easier
  order     : output array of size N
  Runtime: Θ(N)
  */
  std::vector<pos_node> queue(N/2);
  
  int parent;
  int child;
  int pos = 0;
  IntegerVector order (N);
  queue[0].pos = 0;
  queue[0].node = N-2;
  int idx = 1;
  
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
