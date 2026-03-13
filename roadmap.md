# lkj — Development Roadmap

This document records the three-stage implementation plan for the `lkj`
R package and tracks the status of each milestone.

---

## Stage 1 — Pure-R Implementation ✅

Goal: implement and test all four core functions in base R so that the
algorithms are correct, readable, and fully covered by unit tests before any
performance work begins.

### Functions

| Function | Method | Type | Status |
|---|---|---|---|
| `cvine_cholesky(v, d, eta)` | C-vine | Deterministic | ✅ |
| `rlkj_cvine(n, d, eta)` | C-vine | Random | ✅ |
| `onion_cholesky(v, d, eta)` | Onion | Deterministic | ✅ |
| `onion_corr(d, eta)` | Onion | Random | ✅ |
| `ronion(n, d, eta)` | Onion | Random (batch) | ✅ |

### Tasks

- [x] `cvine_cholesky` — map d(d−1)/2 reals to a Cholesky factor via partial
  correlations and the Yule–Kendall recursion (LKJ 2009, Section 2.4).
- [x] `rlkj_cvine` — draw random correlation matrices using `cvine_cholesky`
  with i.i.d. standard-logistic inputs.
- [x] `onion_cholesky` — deterministic onion mapping: r₁₂ via symmetric Beta
  quantile; each subsequent step k uses one Beta-quantile real (radial) and k
  standard-normal-projection reals (directional). Total input length:
  d(d+1)/2 − 2 for d ≥ 2.
- [x] `onion_corr` — single random correlation matrix via the extended onion
  method (Ghosh & Henderson 2003; LKJ 2009).
- [x] `ronion` — batch wrapper around `onion_corr`; returns a d×d×n array.
- [x] Tests in `tests/testthat/`:
  - `test-cvine-cholesky.R` — structure, validity, edge cases, input
    validation, and `rlkj_cvine` for the C-vine functions.
  - `test-onion.R` — structure, validity, edge cases, input validation, and
    determinism for all onion functions.
- [x] `NAMESPACE` updated to export all Stage 1 functions.
- [x] `README.md` with goals, function table, and quick examples.
- [x] `METHODS.md` with full mathematical derivation of both methods.

### Mathematics summary

**C-vine (deterministic)**

For a d×d matrix, d(d-1)/2 real parameters v₁,…,vₘ are converted to
partial correlations ρ_{k,l|1:(k-1)} via

```
u = logistic(vᵢ),   ρ = 2 · qBeta(u, φₖ, φₖ) − 1,
```

where φₖ = η + (d − k − 1)/2.  The Cholesky factor L is then built
row by row via the Yule–Kendall recursion and a forward solve.

**Onion (deterministic)**

For a d×d matrix, d(d+1)/2 − 2 real parameters are consumed sequentially:

1. r₁₂ = 2 · qBeta(logistic(v₁), β, β) − 1,  β = η + (d−2)/2.
2. For k = 2, …, d−1 (adding row k+1):
   - y = qBeta(logistic(vₚ), k/2, β),  β ← β − 1/2.
   - u = normalise(qNorm(logistic(vₚ₊₁)), …, qNorm(logistic(vₚ₊ₖ))).
   - New row of L: [√y · u,  √(1−y)].

The key identity is that the new Cholesky row equals w = √y · u directly
(the A⁻¹z cancellation), ensuring diag(L L ᵀ) = 1 exactly.

---

## Stage 2 — C++ / Rcpp Migration ✅

Goal: port the pure-R algorithms to efficient C++ via Rcpp, keeping the R
wrappers as thin pass-throughs.

### Tasks

- [x] Port `cvine_cholesky` logic to `src/cvine_cholesky.cpp`.
- [x] Port `onion_cholesky` logic to `src/onion_cholesky.cpp`.
- [x] Regenerate `R/RcppExports.R` with new C++ exports.
- [x] Refactor `R/onion_cholesky.R` to delegate to the C++ backend.
- [x] Ensure all Stage 1 tests still pass against the C++ back-end.

---

## Stage 3 — Unified API ✅

Goal: provide a single `lkj()` / `rlkj()` entry point that dispatches to the
chosen method and returns either L or R.

### Proposed API

```r
# Deterministic mapping
lkj(v, d, eta = 1, method = c("cvine", "onion"), output = c("R", "L"))

# Random generation
rlkj(n, d, eta = 1, method = c("cvine", "onion"), output = c("R", "L"))
```

### Tasks

- [x] Implement `lkj()` dispatcher (`R/lkj.R`).
- [x] Implement `rlkj()` dispatcher (`R/rlkj.R`).
- [x] Add `output = "L"` option to return the Cholesky factor directly.
- [x] Write integration tests for the unified API (`tests/testthat/test-lkj-rlkj.R`).
- [x] Write R-vs-C++ consistency tests (`tests/testthat/test-rcpp-consistency.R`).
- [x] Update `README.md` and `man/` documentation.
- [ ] Tag v0.1.0 release.
