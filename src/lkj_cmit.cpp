#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>


inline double cpc(NumericVector numbers) {


  std::transform(numbers.begin(), numbers.end(),
                 numbers.begin(), [](double x) {
                   return std::sqrt(1 - x*x); // Square each element
                 });

  double product = std::accumulate(numbers.begin(),
                                   numbers.end(),
                                   1.0, std::multiplies<double>());
  return product;
}


// [[Rcpp::export(.lkj_cmit)]]
NumericMatrix lkj_cmit(int k, NumericVector y) {

  NumericVector y_map(k*(k-1)/2);

  std::transform(y.begin(), y.end(),
                 y_map.begin(), [](double x) {
                   return std::tanh(x); // tanh y
                 });

  // allocate the matrix we will return
  NumericMatrix z(k, k);
  NumericMatrix w(k, k);

  int counter = 0;

  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      if(i > j){
        z(i, j) = 0.0;
      } else if(i == j){
        z(i, j) = 1.0;
      }
      else if(i < j){
        z(i, j) = y_map[counter++];}
    }
  }

  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      if(i > j){
        w(i, j) = 0;
      } else if(i == 0 && i == j){
        w(i, j) = 1;
      } else if(0 < i && i == j){
        NumericMatrix subM = z(Range(0, i-1), Range(j, j));
        NumericVector vec_col = subM(_, 0);
        double res = cpc(vec_col);
        w(i, j) = res;
      } else if(0 == i && i < j ){
        w(i, j) = z(i, j);
      } else if(0 < i && i < j){
        NumericMatrix subM = z(Range(0, i-1), Range(j, j));
        NumericVector vec_col = subM(_, 0);
        double res = cpc(vec_col);
        w(i, j) = z(i, j) * res;
      }
    }
  }

  return(w);
}
