#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

inline double cpc_inv(NumericVector numbers) {


  std::transform(numbers.begin(), numbers.end(),
                 numbers.begin(), [](double x) {
                   return 1/std::sqrt(1 - x*x); //
                 });

  double product = std::accumulate(numbers.begin(),
                                   numbers.end(),
                                   1.0, std::multiplies<double>());
  return product;
}

// [[Rcpp::export(.lkj_par)]]
NumericVector lkj_par(NumericMatrix w) {

  int n_rows = w.nrow();
  int n_cols = w.ncol();

  // allocate the matrix we will return
  NumericMatrix z(w.nrow(), w.ncol());


  int k = 0;
  NumericVector output_vect(n_rows*(n_rows-1)/2);

  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {

      if(i >= j){
        z(i, j) = 0L;
      } else if(i == 0 && i < j){
        z(i, j) = w(i, j);
        output_vect[k++] = z(i, j);
      }
      else if(0 < i && i < j){
        NumericMatrix subM = z(Range(0, i-1), Range(j, j));
        NumericVector vec_col = subM(_, 0);
        double res = cpc_inv(vec_col);
        z(i, j) = w(i, j) * res;
        output_vect[k++] = z(i, j);
      }



    }
  }

  std::transform(output_vect.begin(), output_vect.end(),
                 output_vect.begin(), [](double x) {
                   return std::atanh(x); // Square each element
                 });
  return output_vect;
}
