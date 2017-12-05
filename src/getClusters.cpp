#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector vCalculeClassif (
    const NumericVector &mesLambdas,
    const NumericVector &lambdas   ,
    const IntegerVector &idown     ,
    const IntegerVector &iup       ,
    const IntegerVector &fils      ,
    const IntegerVector &class0,
    const int &longueurClassif) {

  int longueurMesLambdas = mesLambdas.size() ;
  int longueurLambdas = lambdas.size()       ;
  int longueurClassifFinale = longueurMesLambdas*longueurClassif;
  std::cout << "\n" << longueurClassif << std::endl;
  IntegerVector clusters(longueurMesLambdas*longueurClassif) ;
  IntegerVector class0_(class0);

  int i = 0, n = 0, o = 0 ;
  
  //On initialise le vecteur final qui contiendra toutes les classifs
  for(int l=0;l<longueurClassifFinale;l++) 
    clusters[l]=0 ;
  
  for(int k=0; k < longueurMesLambdas; k++) {
    
    while (mesLambdas[k] < lambdas[i] & i<longueurLambdas) {
      
      if(fils[i] != -1) {
        for(int j2=idown[i];j2<=iup[fils[i]-1];j2++) {
          class0_[j2-1] = i+2 ;
        }
      } else
        class0_[idown[i]-1] = i+2 ;
      
      i++;
    } 
    o = 0 ;
    
    for (int m=n;m<n+longueurClassif;m++){
      clusters[m]=class0_[o];
      o++;
    }
    n=n+longueurClassif;
  }
  return(clusters);
}

// [[Rcpp::export]]
IntegerVector CalculeFils (const IntegerVector &idown, int longueurClassif) {
  
  int longueurFils = idown.size();
  IntegerVector lastAp(longueurClassif), res(longueurFils);
  IntegerVector idown_out = idown;
  
  int temp ;

  for(int i=0;i<longueurClassif;i++)
    lastAp[i]=-1 ;	
  
  for(int j=longueurFils-1;j>=0;j--) {
    temp=idown_out[j] ;
    idown_out[j] = lastAp[idown_out[j]-1] ;
    lastAp[temp-1]=j+1 ;
    
  }
  return(idown_out);  
}
