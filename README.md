# lkj

<!-- badges: start -->
<!-- badges: end -->

An R package for generating and working with correlation matrices from the
**Lewandowski–Kurowicka–Joe (LKJ) distribution**, using both the **C-vine** and
the **extended onion** construction methods.

The LKJ distribution has density proportional to det(**R**)^(η − 1), so η = 1
gives the uniform distribution over all correlation matrices and larger η
concentrates mass near the identity.

---

## Goals

- Provide pure-R and (later) C++ implementations of the LKJ distribution with
  both the **C-vine** and **onion** methods.
- Support **deterministic mappings** (unconstrained reals → Cholesky factor)
  and **random generators** for each method.
- Design a clean, unified API: `lkj()` and `rlkj()`.
- Document the mathematics clearly alongside the code (see `METHODS.md` and
  `roadmap.md`).

---

## Installation

```r
# Install from GitHub (development version)
# install.packages("remotes")
remotes::install_github("alrobles/lkj")
```

---

## Functions

| Function | Method | Type | Description |
|---|---|---|---|
| `lkj(v, d, eta, method, output)` | cvine / onion | Deterministic | Unified dispatcher: maps reals to a Cholesky factor or correlation matrix |
| `rlkj(n, d, eta, method, output)` | cvine / onion | Random | Unified dispatcher: generates n random Cholesky factors or correlation matrices |
| `cvine_cholesky(v, d, eta)` | C-vine | Deterministic | Maps d(d−1)/2 reals to a Cholesky factor L (C++ backend) |
| `rlkj_cvine(n, d, eta)` | C-vine | Random | Generates n random correlation matrices |
| `onion_cholesky(v, d, eta)` | Onion | Deterministic | Maps d(d+1)/2−2 reals to a Cholesky factor L (C++ backend) |
| `onion_corr(d, eta)` | Onion | Random | Generates a single random correlation matrix |
| `ronion(n, d, eta)` | Onion | Random | Generates n random correlation matrices |
| `L_matrix(y)` | C-vine (Rcpp) | Deterministic | Cholesky factor from parameter vector |
| `R_matrix(k, y)` | C-vine (Rcpp) | Deterministic | Correlation matrix from parameter vector |
| `L_par(L)` | C-vine (Rcpp) | Inverse | Parameter vector from Cholesky factor |

---

## Quick Examples

### Unified API (recommended)

```r
library(lkj)

# Deterministic mapping: C-vine, d=4 (need 6 parameters)
v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
R <- lkj(v, d = 4, eta = 1, method = "cvine")
print(round(R, 4))

# Deterministic mapping: Onion, d=4 (need 8 parameters)
v2 <- rep(0, 8)
R2 <- lkj(v2, d = 4, eta = 1, method = "onion")

# Return Cholesky factor directly
L <- lkj(v, d = 4, eta = 1, method = "cvine", output = "L")

# Random generation: single 3x3 C-vine matrix
R <- rlkj(1, d = 3, eta = 2)

# Random generation: ten 4x4 Onion matrices
arr <- rlkj(10, d = 4, eta = 1, method = "onion")
```

### Lower-level functions

### Deterministic mapping (C-vine)

```r
library(lkj)

# d = 4: need d*(d-1)/2 = 6 parameters
v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
L <- cvine_cholesky(v, d = 4, eta = 1)
R <- L %*% t(L)
print(round(R, 4))
```

### Deterministic mapping (Onion)

```r
# d = 4: need d*(d+1)/2 - 2 = 8 parameters
v <- rep(0, 8)
L <- onion_cholesky(v, d = 4, eta = 1)
R <- L %*% t(L)
print(round(R, 4))
```

### Random generation (C-vine)

```r
# Single 3×3 matrix
R <- rlkj_cvine(1, d = 3, eta = 2)

# Ten 4×4 matrices returned as a 4×4×10 array
arr <- rlkj_cvine(10, d = 4, eta = 1)
```

### Random generation (Onion)

```r
# Single 4×4 matrix
R <- onion_corr(d = 4, eta = 2)

# Five 3×3 matrices
arr <- ronion(5, d = 3, eta = 1)
```

---

## Mathematical Background

See [`METHODS.md`](METHODS.md) for the full derivation of both methods,
including the sequential onion construction, the Beta-distribution argument,
the Cholesky factorisation identities, and the deterministic mapping via
logistic/Beta quantile transforms.

---

## Development Roadmap

See [`roadmap.md`](roadmap.md) for the three-stage plan:

| Stage | Goal | Status |
|---|---|---|
| 1 | Pure-R implementations + tests | ✅ Complete |
| 2 | C++ / Rcpp migration | ✅ Complete |
| 3 | Unified `lkj()` / `rlkj()` API | ✅ Complete |

---

## References

- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
  correlation matrices based on vines and extended onion method.
  *Journal of Multivariate Analysis*, 100(9), 1989–2001.
- Joe, H. (2006). Generating random correlation matrices based on partial
  correlations. *Journal of Multivariate Analysis*, 97, 2177–2189.
- Ghosh, S., & Henderson, S. G. (2003). Behavior of the NORTA method for
  correlated random vector generation as the dimension increases.
  *ACM Transactions on Modeling and Computer Simulation*, 13(3), 276–294.

---

## License

GPL (≥ 3) — see [`LICENSE.md`](LICENSE.md).
