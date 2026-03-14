#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// cvine_inverse
//
// Inverse of the C-vine LKJ Cholesky mapping: given a correlation matrix R,
// recover the unconstrained real vector v such that
//   cvine_cholesky(v, d, eta) %*% t(cvine_cholesky(v, d, eta)) == R.
//
// Algorithm
// ---------
// 1. Recover the table of partial correlations p[k, ell] (0-indexed k < ell)
//    from the off-diagonal entries of R using the inverse vine recursion:
//
//    For each column j = 1 .. d-1 (0-indexed):
//      p[0, j] = R[0, j]
//      For i = 1 .. j-1:
//        start with rij = R[i, j]
//        for m = 0 .. i-1 (undo the forward peeling from m=i-1 downto 0):
//          rij = (rij - p[m,i]*p[m,j]) / sqrt((1-p[m,i]^2)*(1-p[m,j]^2))
//        p[i, j] = rij
//
// 2. Map each partial correlation p[k, ell] back to an unconstrained real via:
//      u     = pbeta((p[k,ell] + 1) / 2, phi_k, phi_k)
//      v[idx] = logit(u)
//    where phi_k = eta + (d - k - 2) / 2  (matching the forward direction).
//
// Parameters
// ----------
// R   : d x d correlation matrix (symmetric, positive-definite, unit diagonal).
// d   : dimension (>= 1).
// eta : LKJ shape parameter (eta > 0).
//
// Returns
// -------
// A numeric vector of length d*(d-1)/2 of unconstrained reals v such that
// cvine_cholesky(v, d, eta) produces the Cholesky factor of R.
//
// Reference: Lewandowski D, Kurowicka D, Joe H (2009). "Generating random
// correlation matrices based on vines and extended onion method." Journal of
// Multivariate Analysis, 100(9):1989-2001.

// [[Rcpp::export(cvine_inverse_cpp)]]
NumericVector cvine_inverse(NumericMatrix R, int d, double eta = 1.0) {

  // --- Input validation (basic; full validation done in the R wrapper) --------
  if (d < 1)    Rcpp::stop("d must be >= 1");
  if (eta <= 0) Rcpp::stop("eta must be > 0 (got %g)", eta);

  // Trivial 1 x 1 case
  if (d == 1) return NumericVector(0);

  const double eps   = 1e-12;  // clamping for the unit interval
  const double eps15 = 1e-15;  // clamping for the open interval (-1, 1)

  // --- Step 1: recover partial correlation table p from R --------------------
  // p(k, ell) = partial correlation of variable k+1 and variable ell+1
  //             given {variable 1, ..., variable k}  (0-indexed k < ell).
  NumericMatrix p(d, d); // initialised to 0

  for (int j = 1; j < d; j++) {
    // Level-0 partial equals the unconditional (marginal) correlation.
    p(0, j) = std::min(std::max(R(0, j), -1.0 + eps15), 1.0 - eps15);

    for (int i = 1; i < j; i++) {
      // Unconditional correlation between variable i and variable j.
      double rij = R(i, j);

      // Peel back one conditioning variable at a time, undoing the forward
      // vine recursion that was applied in order m = i-1, i-2, ..., 0.
      // The inverse applies m in the reverse order: m = 0, 1, ..., i-1.
      for (int m = 0; m < i; m++) {
        double rim = p(m, i);
        double rjm = p(m, j);
        double denom = std::sqrt(
          std::max(0.0, (1.0 - rim * rim) * (1.0 - rjm * rjm)));
        if (denom < eps15) {
          // Degenerate case: one of the marginals is exactly ±1, so the
          // conditional correlation is not well-defined; clamp to 0.
          rij = 0.0;
        } else {
          rij = (rij - rim * rjm) / denom;
        }
        rij = std::min(std::max(rij, -1.0 + eps15), 1.0 - eps15);
      }
      p(i, j) = rij;
    }
  }

  // --- Step 2: map partial correlations to unconstrained reals ---------------
  // Mirrors the forward direction: level-major order (k then ell within k).
  int m = d * (d - 1) / 2;
  NumericVector v(m);
  int idx = 0;

  for (int k = 0; k < d - 1; k++) {
    // Shape parameter for the symmetric Beta distribution at level k.
    // Matches the forward formula: phi_k = eta + (d - k - 2) / 2.
    double phi_k = eta + (d - k - 2) / 2.0;

    for (int ell = k + 1; ell < d; ell++) {
      double pval = p(k, ell);
      // Map (-1, 1) -> (0, 1) then apply inverse Beta CDF.
      double u = R::pbeta((pval + 1.0) / 2.0, phi_k, phi_k, 1, 0);
      u = std::min(std::max(u, eps), 1.0 - eps);
      // Logit: maps (0, 1) -> R.
      v[idx++] = std::log(u / (1.0 - u));
    }
  }

  return v;
}
