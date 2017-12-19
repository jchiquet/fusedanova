#include <Rcpp.h>
using namespace Rcpp;

IntegerMatrix get_clustering (
    const NumericVector &set_lambdas,
    const NumericVector &all_lambdas,
    const IntegerVector &idown      ,
    const IntegerVector &iup        ,
    const int K) {

  // get the children from idown  
  IntegerVector lastAp(K, -1) ;
  IntegerVector child = clone(idown) ;
  for(int j = (idown.size()-1); j >= 0; j--) {
    int temp = child[j] ;
    child[j] = lastAp[child[j]-1] ;
    lastAp[temp-1] = j + 1 ;
  }

  // now build the clustering matrix along the tree
  IntegerMatrix clustering(K, set_lambdas.size()) ;
  IntegerVector class0(K, 1);

  int i = 0 ; // current class label
  
  for(int k = 0; k < set_lambdas.size(); k++) {
    
    while (set_lambdas[k] < all_lambdas[i] & i < all_lambdas.size()) {
      if (child[i] == -1) {
        class0[idown[i]-1] = i + 2 ;
      } else {
        for (int j2 = idown[i]; j2 <= iup[child[i]-1]; j2++) {
          class0[j2-1] = i + 2 ;
        }
      }
      i++;
    } 
    clustering(_, k) = class0 ;
  }
  
  return(clustering);
}

