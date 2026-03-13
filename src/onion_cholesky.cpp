#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// onion_cholesky
//
// Build the Cholesky factor L of a correlation matrix using the extended
// onion method (Ghosh & Henderson 2003; Lewandowski, Kurowicka & Joe 2009).
//
// The construction proceeds sequentially: starting from a 1x1 trivial case,
// each step adds one row/column by drawing a radial component y ~ Beta(k/2, beta)
// and a uniform direction on the k-sphere (normalised standard normals).
//
// Parameters
// ----------
// v   : numeric vector of length d*(d+1)/2 - 2 (one unconstrained real per
//       parameter in level-major order). For d=1 use numeric(0); for d=2 use
//       a single real.
// d   : dimension of the target correlation matrix (>= 1).
// eta : LKJ shape parameter (eta > 0). eta = 1 gives the uniform distribution
//       over correlation matrices; larger eta concentrates mass near identity.
//
// Returns
// -------
// A d x d lower-triangular matrix L such that R = L %*% t(L) is a valid
// correlation matrix with det(R) proportional to det(R)^(eta-1).
//
// Reference: Lewandowski D, Kurowicka D, Joe H (2009). "Generating random
// correlation matrices based on vines and extended onion method." Journal of
// Multivariate Analysis, 100(9):1989-2001.

// [[Rcpp::export(onion_cholesky_cpp)]]
NumericMatrix onion_cholesky_impl(NumericVector v, int d, double eta = 1.0) {

  // --- Input validation -----------------------------------------------------
  if (d < 1)    Rcpp::stop("d must be >= 1");
  if (eta <= 0) Rcpp::stop("eta must be > 0 (got %g)", eta);

  // Trivial 1 x 1 case
  if (d == 1) {
    NumericMatrix L(1, 1);
    L(0, 0) = 1.0;
    return L;
  }

  // Expected length of v
  int m_needed = (d == 2) ? 1 : d * (d + 1) / 2 - 2;
  if ((int)v.size() != m_needed) {
    Rcpp::stop("v must have length %d for d=%d (got %d).",
               m_needed, d, (int)v.size());
  }

  const double eps   = 1e-12;  // clamping for the unit interval
  const double eps15 = 1e-15;  // clamping for the open interval

  // Inline logistic (sigmoid) function: maps R -> (0,1)
  auto sigmoid = [&](double x) -> double {
    double u = 1.0 / (1.0 + std::exp(-x));
    return std::min(std::max(u, eps), 1.0 - eps);
  };

  NumericMatrix L(d, d); // initialised to 0

  // L[0,0] = 1 (trivial 1x1 block)
  L(0, 0) = 1.0;

  int idx = 0; // linear index into v

  // Initial beta = eta + (d-2)/2
  double beta = eta + (d - 2) / 2.0;

  // --- Step 1: build the 2x2 block ------------------------------------------
  // r12 from a symmetric Beta on (-1,1) with shape beta
  double u1  = sigmoid(v[idx++]);
  double r12 = 2.0 * R::qbeta(u1, beta, beta, 1, 0) - 1.0;
  r12 = std::min(std::max(r12, -1.0 + eps15), 1.0 - eps15);

  L(1, 0) = r12;
  L(1, 1) = std::sqrt(std::max(0.0, 1.0 - r12 * r12));

  if (d == 2) return L;

  // --- Sequential steps: add row n+1 for n = 1, ..., d-2 -------------------
  for (int step = 1; step < d - 1; step++) {
    beta -= 0.5;                  // decrement beta for each new dimension
    int k    = step + 1;          // dimension of the existing submatrix
    int row_j = step + 1;         // 0-indexed row being filled (= k in 0-indexed)

    // Radial component: y ~ Beta(k/2, beta)
    double u_y = sigmoid(v[idx++]);
    double y   = R::qbeta(u_y, k / 2.0, beta, 1, 0);
    y = std::min(std::max(y, eps15), 1.0 - eps15);

    // Direction on the k-sphere: normalise k standard-normal draws
    // The normals are obtained via the inverse Gaussian CDF applied to
    // logistic-transformed reals, i.e. qnorm(sigmoid(v[...])).
    NumericVector z_raw(k);
    for (int i = 0; i < k; i++) {
      double u_z = sigmoid(v[idx++]);
      z_raw[i]   = R::qnorm(u_z, 0.0, 1.0, 1, 0);
    }

    // Compute norm of z_raw
    double nrm = 0.0;
    for (int i = 0; i < k; i++) nrm += z_raw[i] * z_raw[i];
    nrm = std::sqrt(nrm);

    NumericVector u_vec(k);
    if (nrm < 1e-12) {
      // Degenerate: pick first unit vector
      u_vec[0] = 1.0;
    } else {
      for (int i = 0; i < k; i++) u_vec[i] = z_raw[i] / nrm;
    }

    // New row of L: first k entries are sqrt(y) * u_vec,
    // diagonal entry is sqrt(1 - y).
    double sqrt_y = std::sqrt(y);
    for (int col = 0; col < k; col++) {
      L(row_j, col) = sqrt_y * u_vec[col];
    }
    L(row_j, row_j) = std::sqrt(std::max(0.0, 1.0 - y));
  }

  return L;
}
