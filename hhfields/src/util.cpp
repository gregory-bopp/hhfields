#include <Rcpp.h>
#include <numeric>
#include <math.h>

using namespace Rcpp;


// Calculate Euclidean distance between the rows of X
// [[Rcpp::export]]
NumericMatrix calc_distC(NumericMatrix X){
  int n = X.nrow();
  NumericMatrix D(n, n);
  for(int i=0; i<n; i++){
    for(int j=i; j<n; j++){
      D(i,j) = std::sqrt(sum(pow(X(i,_) - X(j,_), 2.0)));
      D(j,i) = D(i,j);
    }
  }
  return(D);
}

// Calculate Euclidean distance between rows of X and rows of Y
// D(i,j) contains distance between row i of X and row j of Y
// [[Rcpp::export]]
NumericMatrix calc_distXYC(NumericMatrix X, NumericMatrix Y){
  int n_xrows = X.nrow();
  int n_yrows = Y.nrow();
  NumericMatrix D(n_xrows, n_yrows);

  for(int i=0; i<n_xrows; i++){
    for(int j=0; j<n_yrows; j++){
      D(i,j) = std::sqrt(sum(pow(X(i,_) - Y(j,_), 2.0)));
    }
  }
  return(D);
}
