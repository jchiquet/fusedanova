#include <Rcpp.h>
#include <bits/stdc++.h>
#include "nosplit.h"
#include <stdio.h>

using namespace Rcpp ;

class Rule {
  
private:
  int group1 ;
  int group2 ;
  double lambda ;
  
public:
  Rule(int group1_, int group2_, double lambda_) : group1(group1_),  group2(group2_), lambda(lambda_) {} ;

  int get_group1() const {return group1 ;}  
  int get_group2() const {return group2 ;}  
  double get_lambda() const {return lambda ;}
  bool is_active(const LogicalVector& active) const {
    return(active[group1] & active[group2]);
  }
};

class RuleComparator {
public:
  int operator() (const Rule& r1, const Rule& r2) {
    return r1.get_lambda() > r2.get_lambda();
  }
};

double get_lambda(int group1, int group2, const DataFrame& table) {

  NumericVector slope  = table["slope" ] ;
  NumericVector beta   = table["beta"  ] ;
  NumericVector lambda = table["lambda"] ;

  return((beta[group1] - beta[group2] - slope[group1] * lambda[group1] + slope[group2] * lambda[group2]) / (slope[group2] - slope[group1]));
}

void merge(int fusion, double lambda_fusion, int group1, int group2, DataFrame& table) {

  // pointers to columns of the DataFrame
  NumericVector slope        = table["slope"   ] ;
  NumericVector beta         = table["beta"    ] ;
  NumericVector lambda       = table["lambda"  ] ;
  IntegerVector i_low        = table["i_low"   ] ;
  IntegerVector i_split      = table["i_split" ] ;
  IntegerVector i_high       = table["i_high"  ] ;
  IntegerVector grp_low      = table["grp_low" ] ;
  IntegerVector grp_high     = table["grp_high"] ;
  IntegerVector grp_size     = table["grp_size"] ;
  LogicalVector active       = table["active"  ] ;
  LogicalVector has_grp_low  = table["has_grp_low"] ;
  LogicalVector has_grp_high = table["has_grp_high"] ;

  // Updating the fields of the DataFrame of the fusing group
  lambda  [fusion] = lambda_fusion ;
  grp_size[fusion] = grp_size[group1] + grp_size[group2] ;
  slope   [fusion] = (grp_size[group1] * slope[group1] + grp_size[group2] * slope[group2]) / grp_size[fusion] ;
  beta    [fusion] = beta[group1] + (lambda_fusion - lambda[group1]) * slope[group1] ;
  if(i_low[group1] < i_low[group2]) {
    i_split     [fusion] = i_high      [group1];
    i_low       [fusion] = i_low       [group1];
    i_high      [fusion] = i_high      [group2];
    grp_low     [fusion] = grp_low     [group1];
    grp_high    [fusion] = grp_high    [group2];
    has_grp_low [fusion] = has_grp_low [group1];
    has_grp_high[fusion] = has_grp_high[group2];
  } else {
    i_split     [fusion] = i_high      [group2]; 
    i_low       [fusion] = i_low       [group2];
    i_high      [fusion] = i_high      [group1];
    grp_low     [fusion] = grp_low     [group2];
    grp_high    [fusion] = grp_high    [group1];
    has_grp_low [fusion] = has_grp_low [group2];
    has_grp_high[fusion] = has_grp_high[group1];
  }
  
  // Updating neighbors
  if(has_grp_low [fusion])  grp_high[  grp_low [fusion] ] = fusion;
  if(has_grp_high [fusion])  grp_low[ grp_high [fusion] ] = fusion;

  active[group1] = false ;
  active[group2] = false ;
  active[fusion] = true ;

  // std::cout << " beta     " << beta[fusion]     << "\t" 
  //           << " lambda   " << lambda[fusion]   << "\t" 
  //           << " i_low    " << i_low[fusion]    << "\t" 
  //           << " i_high   " << i_high[fusion]   << "\t" 
  //           << " i_split  " << i_split[fusion]  << "\t" 
  //           << " grp_low  " << grp_low[fusion]  << "\t" 
  //           << " grp_high " << grp_high[fusion] << std::endl ;
}

//' @export
// [[Rcpp::export]]
DataFrame fuse(NumericVector beta0, NumericVector slope0, IntegerVector grp_size0) {

  // initialize  
  int n = grp_size0.size() ;
  NumericVector lambda       (2 * n - 1) ; 
  NumericVector beta         (2 * n - 1) ;
  NumericVector slope        (2 * n - 1) ;
  IntegerVector i_low        (2 * n - 1) ;
  IntegerVector i_split      (2 * n - 1) ;
  IntegerVector i_high       (2 * n - 1) ;
  IntegerVector grp_low      (2 * n - 1) ;
  IntegerVector grp_high     (2 * n - 1) ;
  IntegerVector grp_size     (2 * n - 1) ;
  LogicalVector active       (2 * n - 1) ;
  LogicalVector has_grp_low  (2 * n - 1) ;
  LogicalVector has_grp_high (2 * n - 1) ;
  
  DataFrame table = DataFrame::create(Named("lambda")       = lambda      ,
                                      Named("beta")         = beta        ,
                                      Named("slope")        = slope       ,
                                      Named("i_low")        = i_low       ,
                                      Named("i_high")       = i_high      ,
                                      Named("i_split")      = i_split     ,
                                      Named("grp_low")      = grp_low     ,
                                      Named("grp_high")     = grp_high    ,
                                      Named("has_grp_low")  = has_grp_low ,
                                      Named("has_grp_high") = has_grp_high,
                                      Named("grp_size")     = grp_size    ,
                                      Named("active")       = active      ) ;
    
  for (int k=0; k < n; k++) {
    lambda[k]   = 0.0;
    beta[k]     = beta0[k];
    slope[k]    = slope0[k];
    i_low[k]    = k;
    i_split[k]  = k;
    i_high[k]   = k;
    grp_size[k] = grp_size0[k] ;
    active[k]   = true;
    if (k == 0) {
      has_grp_low[k] = false ;
      grp_low[k] = -1;
    } else {
      has_grp_low[k] = true ;
      grp_low[k] = k - 1;
    }
    if (k == (n-1)) {
      has_grp_high[k] = false ;
      grp_high[k] = -1 ;
    } else {
      has_grp_high[k] = true ;
      grp_high[k] = k + 1;
    }
  }
  
  priority_queue <Rule, vector<Rule>, RuleComparator> myMinHeap ;
  
  // first round of fusion between the n-1 pairs of successive groups
  for (int k=0; k < (n-1); k++) {
    myMinHeap.push(Rule(k, k + 1, get_lambda(k, k + 1, table)));
  }
  
  for (int k = n-1; k < (2*n-1);  k++) {
    
    // find the next active rule
    while (!myMinHeap.top().is_active(active) & !myMinHeap.empty()) {
      myMinHeap.pop() ;
      R_CheckUserInterrupt();
    }
    std::cout << "fusion = " << k - (n-1) + 1 
              << ", heapSize = " << myMinHeap.size() 
              << ", active groups = " << sum(active) << std::endl;
    if (myMinHeap.empty()) {
      std::cout << "DIANTRE: no more active rule (with both group actives...)" << std::endl;
      break;
    }
    Rule rule_ = myMinHeap.top();
    myMinHeap.pop() ;

    // merge the two groups and unactivate the fused groups
    merge(k, rule_.get_lambda(), rule_.get_group1(), rule_.get_group2(), table) ;

    // get new rules an add them to the queue
    if (has_grp_low[k]) { //& active[grp_low[k]]) {
      myMinHeap.push(Rule(grp_low[k], k, get_lambda(grp_low[k], k, table)));
    }
    
    if (has_grp_high[k]){ // & active[grp_high[k]]) {
      myMinHeap.push(Rule(k, grp_high[k], get_lambda(k, grp_high[k], table)));
    }
  }

  return(DataFrame::create(Named("beta"  , beta),
                           Named("lambda", lambda),
                           Named("slope" , slope),
                           Named("down"  , i_low + 1),
                           Named("high"  , i_high + 1),
                           Named("grp_down"  , grp_low + 1),
                           Named("grp_high"  , grp_high + 1),
                           Named("split" , i_split + 1)));
  
}

//' @export
// [[Rcpp::export]]
DataFrame fuse_old(NumericVector x, NumericVector slopes, NumericVector ngroup){

	Group *G = maketree(&x[0], x.length(), &slopes[0], &ngroup[0]);

  int nrow = 2*G->len -1; // nb infos = K init gpes + K-1 fusions 
  
  NumericVector beta(nrow), lambda(nrow), slope(nrow);
  IntegerVector idown(nrow), iup(nrow);
  
  int row=0;
  
  add_results(G, &beta[0], &lambda[0], &slope[0], &idown[0], &iup[0], &row);
  
  DataFrame L = DataFrame::create(Named("beta"  , beta),
                                  Named("lambda", lambda),
                                  Named("slope" , slope),
                                  Named("idown" , idown),
                                  Named("iup"   , iup));
	
	delete_tree(G); // delete the old tree for memory

	return L;
}

//             beta       lambda slope idown iup
// 1  -4.857226e-17 2.705679e-02     0     1  40
// 2  -1.463452e-02 1.973953e-02     2     1  39
// 3  -3.274782e-02 1.521120e-02     4     1  38
// 4  -7.122460e-02 9.714520e-03     7     1  37
// 5  -9.194006e-02 7.642974e-03    10     1  36
// 6  -1.143391e-01 7.322988e-03    70     1  30
// 65 -5.381602e-02 7.276397e-03  -104    31  36
// 7  -1.388752e-01 6.982209e-03    72     1  28
// 8  -2.276016e-01 5.985283e-03    89     1  25
// 9  -3.157023e-01 5.027667e-03    92     1  24
// 69  2.557720e-01 4.674817e-03  -119    33  36
// 57 -7.606477e-02 4.251323e-03   -23    26  28
// 58 -7.219435e-02 4.036300e-03   -18    26  27
// 10 -5.262880e-01 2.810975e-03    95     1  23
// 11 -6.421649e-01 1.717796e-03   106     1  22
// 62  1.232100e-01 1.667057e-03   -42    29  30
// 12 -6.561043e-01 1.587522e-03   107     1  21
// 13 -6.755379e-01 1.407581e-03   108     1  20
// 14 -7.077331e-01 1.112213e-03   109     1  19
// 15 -7.325156e-01 8.869165e-04   110     1  18
// 66  3.236904e-01 8.779836e-04   -59    31  32
// 16 -7.412300e-01 8.084087e-04   111     1  17
// 18 -7.442028e-01 7.811355e-04   109     2  17
// 20 -7.467669e-01 7.573936e-04   108     3  17
// 70  7.234843e-01 6.428142e-04  -116    33  35
// 21 -7.669999e-01 5.783407e-04   113     3  16
// 22 -7.769732e-01 4.908558e-04   114     3  15
// 23 -7.911740e-01 3.673701e-04   115     3  14
// 25 -7.940296e-01 3.418741e-04   112     4  14
// 26 -7.965488e-01 3.207041e-04   119     4  12
// 72  7.634002e-01 3.016528e-04  -117    34  35
// 27 -8.028539e-01 2.690232e-04   122     4  11
// 39 -8.057345e-01 2.366568e-04    89    10  11
// 28 -8.079476e-01 2.323777e-04   139     4   9
// 30 -8.091022e-01 2.240108e-04   138     5   9
// 32 -8.241863e-01 1.079799e-04   130     6   9
// 33 -8.269726e-01 8.687101e-05   132     6   8
// 34 -8.334388e-01 4.068414e-05   140     6   7
// 43 -8.123267e-01 3.175347e-05    59    13  14
// 17 -8.883604e-01 0.000000e+00   182     1   1
// 19 -8.840260e-01 0.000000e+00   179     2   2
// 24 -8.554638e-01 0.000000e+00   175     3   3
// 29 -8.476842e-01 0.000000e+00   171     4   4
// 31 -8.453920e-01 0.000000e+00   162     5   5
// 35 -8.394600e-01 0.000000e+00   148     6   6
// 36 -8.388905e-01 0.000000e+00   134     7   7
// 37 -8.372234e-01 0.000000e+00   118     8   8
// 38 -8.358481e-01 0.000000e+00   108     9   9
// 40 -8.303468e-01 0.000000e+00   104    10  10
// 41 -8.263236e-01 0.000000e+00    87    11  11
// 42 -8.186774e-01 0.000000e+00    69    12  12
// 44 -8.142319e-01 0.000000e+00    60    13  13
// 45 -8.140096e-01 0.000000e+00    53    14  14
// 46 -8.020068e-01 0.000000e+00    51    15  15
// 47 -7.953386e-01 0.000000e+00    49    16  16
// 48 -7.793348e-01 0.000000e+00    43    17  17
// 49 -7.653316e-01 0.000000e+00    37    18  18
// 50 -7.466605e-01 0.000000e+00    35    19  19
// 51 -7.219881e-01 0.000000e+00    33    20  20
// 52 -7.053175e-01 0.000000e+00    31    21  21
// 53 -6.919810e-01 0.000000e+00    29    22  22
// 54 -5.740746e-01 0.000000e+00    17    23  23
// 55 -3.307853e-01 0.000000e+00     3    24  24
// 56 -2.096458e-01 0.000000e+00    -3    25  25
// 59 -3.577255e-03 0.000000e+00   -17    26  26
// 60  4.485834e-02 0.000000e+00   -29    27  27
// 61  7.273154e-02 0.000000e+00   -35    28  28
// 63  1.915594e-01 0.000000e+00   -41    29  29
// 64  1.948935e-01 0.000000e+00   -43    30  30
// 67  3.711015e-01 0.000000e+00   -54    31  31
// 68  3.842712e-01 0.000000e+00   -69    32  32
// 71  7.716954e-01 0.000000e+00   -75    33  33
// 73  7.875324e-01 0.000000e+00   -80    34  34
// 74  7.999002e-01 0.000000e+00  -121    35  35
// 75  1.008418e+00 0.000000e+00  -161    36  36
// 76  1.551100e+00 0.000000e+00  -167    37  37
// 77  2.598790e+00 0.000000e+00  -173    38  38
// 78  3.499001e+00 0.000000e+00  -178    39  39
// 79  4.924335e+00 0.000000e+00  -182    40  40

// we process each dimension individually using this function
// RcppExport SEXP noSplitcv(SEXP R_x,SEXP R_xv,SEXP R_ngroup, SEXP R_xtest,SEXP R_ngrouptest ,SEXP R_args) {
// 
// 	NumericVector x(R_x);
// 	NumericVector xv(R_xv);
// 	NumericVector xtest(R_xtest);
// 	NumericVector ngroup(R_ngroup);
// 	NumericVector ngrouptest(R_ngrouptest);
// 	List args(R_args);
//   
// 	std::string weights = Rcpp::as<std::string>(args["weights"]); 
// 	double gamma = Rcpp::as<double>(args["gamma"]); 
// 	NumericMatrix W = args["W"];
// 	NumericVector lambdalist = args["lambdalist"];
// 
// 	NumericVector error(lambdalist.length());
// 
// 	vector<double> sl = calculateSlope(x,ngroup,xv,weights,gamma,W,x.length());
//  
// 	Group *G = maketree(&x[0], x.length(), &sl[0],&ngroup[0]);
// 
// 	error_cv(G,&lambdalist[0],lambdalist.length(),&xtest[0], &ngrouptest[0],&error[0]);
// 
// 	delete_tree(G); 
// 
// 	return(error);
// }
// 
