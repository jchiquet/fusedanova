#include <Rcpp.h>
#include <stdlib.h>
#include <vector>
using namespace Rcpp;
using namespace std;

# define square(x) ((x)*(x))
# define NEW_EVENT_THRESH 1e-16

//' @export
// [[Rcpp::export]]
NumericVector get_slopes(NumericVector &xm    ,
                         IntegerVector &ngroup,
                         NumericVector &xv    , 
                         std::string weights       ,
                         double gamma         ,
                         NumericMatrix &W     ) {
  
  int n = xm.length()  ; 
  int N = sum(ngroup)  ; // number of individuals in all groups
  int sign             ; // some auxiliary variables
  double sum1=0,sum2=0 ;
  std::vector<double> slopes; // the vector of slopes (output of the function)
  
  if (weights=="default"){ // easiest
    
    if (N==n){ //slopes are known (clusterpath case)
      for (int i =0; i<n; i++){
        slopes.push_back(n-1-2*i);
      }
    } else { // nk.nl weights (O(n) calculation)
      for (int i=1;i<n;i++){
        sum1 += ngroup[i];
      }
      slopes.push_back(sum1);
      for (int i =0; i<(n-1);i++){
        slopes.push_back(slopes[i]-ngroup[i]-ngroup[i+1]);
      }
    } 
  } else if (weights == "laplace"){ // O(n) calculation
    sum1 = 0;
    sum2 = 0;
    vector<double> slopes1(n),slopes2(n);
    
    for (int i=n-1;i>0;i--){
      sum1 += ngroup[i]*exp(-gamma * xm[i]);
      slopes1[i-1] = sum1;
    }
    slopes1[n-1] =0;
    
    for (int i=1;i<n;i++){
      sum2 += ngroup[i-1]*exp(gamma * xm[i-1]);
      slopes2[i] = sum2;
    }
    slopes2[0] =0;	
    
    for (int i=0;i<n;i++){
      slopes.push_back(slopes1[i]*exp(gamma * xm[i])- slopes2[i]*exp(-gamma * xm[i]));
    }		
    
    // O(n2) calculation for all the rest
    // gaussian and adaptive trick : summing starting by the smallest element for each line
    // the smallest are the one away from the "diagonal"
    // create two vectors of sum and then add them
  } else if (weights == "gaussian"){
    // calculate the positive part
    for(int i=0;i<(n-1);i++){
      sum1=0;
      for(int j =n-1;j>i;j--){
        sum1 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
      }
      slopes.push_back(sum1);
    }
    slopes.push_back(0);
    // calculate the negative part
    
    for(int i=1;i<n;i++){
      sum2=0;
      for (int j = 0;j<i ;j++){
        sum2 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
      }
      slopes[i]-=sum2;
    }
    
  } else if (weights == "adaptive"){
    // calculate the positive part
    for(int i=0;i<(n-1);i++){
      sum1=0;
      for(int j =n-1;j>i;j--){
        sum1 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
      }
      slopes.push_back(sum1);
    }
    slopes.push_back(0);
    // calculate the negative part
    
    for(int i=1;i<n;i++){
      sum2=0;
      for (int j = 0;j<i ;j++){
        sum2 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
      }
      slopes[i]-=sum2;
    }
    
  } else if (weights == "naivettest"){ 
    for (int i=0; i<n; i++){
      sum1 = 0 ;
      sign = -1 ;
      for (int j=0; j<n; j++){
        if (i==j){
          sign = 1 ;
        }else{
          sum1+= sign / sqrt(1/ngroup[i]+1/ngroup[j]) ;
        }
      } 
      slopes.push_back(sum1/ ngroup[i]);
    } 	
  } else if (weights == "ttest"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1;
      for (int j=0; j<n; j++){
        if (i==j){
          sign= 1;
        }else{
          sum1+= sign * sqrt(((ngroup[i]-1)*xv[i] + (ngroup[j]-1)*xv[j])/(ngroup[i]+ngroup[j]-2)) /sqrt(1/ngroup[i]+1/ngroup[j]) ;
        }
      } 
      slopes.push_back(sum1/ ngroup[i]);
    } 	
  } else if (weights == "welch"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1 ;
      for (int j=0; j<n; j++){
        if (i==j){
          sign= 1 ;
        }else{
          sum1+= sign * sqrt(xv[i]/ngroup[i]+xv[j]/ngroup[j])/(1/ngroup[i]+1/ngroup[j]) ;
        }
      } 
      slopes.push_back(sum1/ngroup[i]);
    } 	
  } else if (weights == "personal"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1 ;
      for (int j=0; j<n; j++){
        if (i==j){
          sign= 1 ;
        }else{
          sum1 += sign * W(i,j);
        }
      } 
      slopes.push_back(sum1/ngroup[i]);
    } 	
  }
  
  return(wrap(slopes));
}