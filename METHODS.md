# The extended onion method for generating random correlation matrices
## Lewandowski, Kurowicka & Joe (2009)

This document provides a step‑by‑step mathematical explanation of the **extended onion method** introduced by Ghosh & Henderson (2003) and extended by Lewandowski, Kurowicka & Joe (2009).  
Like the C‑vine method, it produces correlation matrices whose density is proportional to  

\[
\det(\mathbf{R})^{\eta-1}, \qquad \eta>0,
\]

with the uniform distribution (\(\eta = 1\)) as a special case.  
The method builds the matrix sequentially, adding one variable at a time, and can be used both for **random generation** and as a **deterministic mapping** from a vector of real numbers to a correlation matrix.

---

## 1. Sequential construction of a correlation matrix

Let \(\mathbf{R}_k\) denote a \(k\times k\) correlation matrix (positive definite, unit diagonal).  
We want to extend it to a \((k+1)\times(k+1)\) correlation matrix

\[
\mathbf{R}_{k+1} = \begin{pmatrix}
\mathbf{R}_k & \mathbf{z} \\
\mathbf{z}^\top & 1
\end{pmatrix},
\]

where \(\mathbf{z} = (z_1,\dots,z_k)^\top\) is a vector of correlations between the new variable and the existing \(k\) variables.  
For \(\mathbf{R}_{k+1}\) to be a valid correlation matrix, we need \(\mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z} < 1\).

A key identity links the determinants:

\[
\det(\mathbf{R}_{k+1}) = \det(\mathbf{R}_k)\;\bigl(1 - \mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z}\bigr).
\tag{1}
\]

Thus the “new information” brought by the \((k+1)\)-st variable is the scalar  

\[
y = \mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z} \in (0,1).
\]

---

## 2. Conditional distribution of \(\mathbf{z}\) given \(\mathbf{R}_k\)

The onion method is based on an **elliptically contoured** construction.  
Suppose \(\mathbf{R}_k\) has density proportional to \(\det(\mathbf{R}_k)^{\beta_k-1}\) for some \(\beta_k>0\).  
If we then choose \(\mathbf{z}\) from the conditional density

\[
f(\mathbf{z} \mid \mathbf{R}_k) \;\propto\; \det(\mathbf{R}_k)^{-1/2}\,
      \bigl(1 - \mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z}\bigr)^{\beta-1},
\tag{2}
\]

where \(\beta>0\) is a parameter, then by (1) the joint density of \((\mathbf{R}_k,\mathbf{z})\) becomes

\[
f(\mathbf{R}_k,\mathbf{z}) \propto
      \det(\mathbf{R}_k)^{\beta_k-3/2}\,
      \bigl(1 - \mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z}\bigr)^{\beta-1}.
\]

Now define the new matrix \(\mathbf{R}_{k+1}\) via the partition above.  
Using (1) we can replace \(\bigl(1 - \mathbf{z}^\top\mathbf{R}_k^{-1}\mathbf{z}\bigr)\) by \(\det(\mathbf{R}_{k+1})/\det(\mathbf{R}_k)\).  
Hence

\[
f(\mathbf{R}_{k+1}) \propto
      \det(\mathbf{R}_k)^{\beta_k-3/2}\,
      \bigl(\det(\mathbf{R}_{k+1})/\det(\mathbf{R}_k)\bigr)^{\beta-1}
      = \det(\mathbf{R}_k)^{\beta_k-\beta-1/2}\,
        \det(\mathbf{R}_{k+1})^{\beta-1}.
\]

To make the final density depend **only** on \(\det(\mathbf{R}_{k+1})\), we set

\[
\beta_k = \beta + \tfrac12.
\]

Then

\[
f(\mathbf{R}_{k+1}) \propto \det(\mathbf{R}_{k+1})^{\beta-1}.
\]

Thus by choosing the parameters appropriately at each step, we maintain the same power‑law form for the determinant as we increase the dimension.

---

## 3. Generating \(\mathbf{z}\) via a spherical representation

Condition (2) is exactly the density of an elliptically contoured vector.  
A standard stochastic representation is

\[
\mathbf{z} = \mathbf{A}\,\mathbf{w},
\]

where \(\mathbf{A}\) is the lower‑triangular Cholesky factor of \(\mathbf{R}_k\) (so \(\mathbf{A}\mathbf{A}^\top = \mathbf{R}_k\)), and \(\mathbf{w}\) is a vector drawn from a **spherically symmetric** distribution with radial part \(r = \|\mathbf{w}\|\) having density

\[
f_r(r) \propto r^{k-1} (1-r^2)^{\beta-1}, \qquad 0<r<1.
\]

Equivalently, let \(y = r^2 = \|\mathbf{w}\|^2\). Then \(y\) follows a **Beta** distribution:

\[
y \sim \mathrm{Beta}\!\left(\frac{k}{2},\,\beta\right).
\tag{3}
\]

The direction \(\mathbf{u} = \mathbf{w} / \|\mathbf{w}\|\) is independent of \(y\) and uniformly distributed on the surface of the \(k\)-dimensional sphere.

Thus we can generate \(\mathbf{z}\) in three simple steps:

1. Generate \(y\) from \(\mathrm{Beta}(k/2,\beta)\).
2. Generate a uniform direction \(\mathbf{u}\) on the sphere (e.g., by normalising a vector of i.i.d. standard normals).
3. Set \(\mathbf{w} = \sqrt{y}\,\mathbf{u}\) and then \(\mathbf{z} = \mathbf{A}\mathbf{w}\).

---

## 4. The extended onion algorithm

We now assemble the complete algorithm for generating a \(d\times d\) correlation matrix with density proportional to \(\det(\mathbf{R})^{\eta-1}\).

**Initialisation** (for \(d=2\)):  
Set \(\beta = \eta + \frac{d-2}{2}\).  
Generate the single correlation \(r_{12}\) from a symmetric Beta distribution on \((-1,1)\) with shape \(\beta\):

\[
r_{12} = 2\,\mathrm{Beta}(\beta,\beta) - 1.
\]

Form \(\mathbf{R}_2 = \begin{pmatrix}1 & r_{12}\\ r_{12} & 1\end{pmatrix}\).

**Sequential steps** for \(k = 2,3,\dots,d-1\) (current size \(k\)):

- Update \(\beta \leftarrow \beta - \frac12\).
- Generate \(y \sim \mathrm{Beta}\!\left(\frac{k}{2},\,\beta\right)\).
- Generate a uniform direction \(\mathbf{u}\) on the \(k\)-sphere.
- Compute \(\mathbf{w} = \sqrt{y}\,\mathbf{u}\).
- Compute the lower Cholesky factor \(\mathbf{A}\) of the current \(\mathbf{R}_k\) (so that \(\mathbf{A}\mathbf{A}^\top = \mathbf{R}_k\)).
- Set \(\mathbf{z} = \mathbf{A}\mathbf{w}\) – this is the vector of correlations with the new variable.
- Form the augmented matrix
  \[
  \mathbf{R}_{k+1} = \begin{pmatrix}
  \mathbf{R}_k & \mathbf{z} \\
  \mathbf{z}^\top & 1
  \end{pmatrix}.
  \]

After \(d-1\) steps we obtain the desired \(d\times d\) correlation matrix \(\mathbf{R}\).

---

## 5. Deterministic mapping from reals to a correlation matrix

The onion method can also be used as a **deterministic transformation** from a vector of real numbers \(\mathbf{v}\) of length \(d(d-1)/2\) to a correlation matrix.  
This is useful when one needs an invertible mapping (e.g., for optimisation or MCMC on an unconstrained space).

The mapping follows the same sequential construction, but replaces random draws by quantile functions applied to transformed reals.  
Let \(\Phi\) be the logistic function \(\Phi(x)=1/(1+e^{-x})\) mapping \(\mathbb{R}\) to \((0,1)\).

- **For the first correlation** \(r_{12}\) we use the first real \(v_1\):
  \[
  u_1 = \Phi(v_1), \qquad
  r_{12} = 2\,q_{\mathrm{Beta}(\beta_1,\beta_1)}(u_1) - 1,
  \]
  where \(\beta_1 = \eta + (d-2)/2\).

- **For each subsequent step** \(k = 2,\dots,d-1\) (current size \(k\)), we need:
  * one real for the Beta‑distributed \(y\),
  * \(k\) reals for the direction \(\mathbf{u}\).

  The reals are consumed in order: first the one for \(y\), then the \(k\) for the direction.

  Let the current index into \(\mathbf{v}\) be \(p\).  
  - Set \(u_y = \Phi(v_{p+1})\) and compute
    \[
    y = q_{\mathrm{Beta}(k/2,\,\beta)}(u_y),
    \]
    where \(\beta\) is the current value (\(\beta = \eta + (d-2)/2 - (k-1)/2\)).
  - For the direction, take the next \(k\) reals \(v_{p+2},\dots,v_{p+k+1}\), transform them to independent standard normals via the inverse Gaussian CDF:
    \[
    z_i = \Phi^{-1}_{\mathcal{N}(0,1)}\bigl(\Phi(v_{p+1+i})\bigr),\quad i=1,\dots,k,
    \]
    then normalise:
    \[
    \mathbf{u} = \frac{(z_1,\dots,z_k)^\top}{\sqrt{z_1^2+\cdots+z_k^2}}.
    \]
  - Proceed as before: \(\mathbf{w} = \sqrt{y}\,\mathbf{u}\), \(\mathbf{z} = \mathbf{A}\mathbf{w}\), augment \(\mathbf{R}\).

Because each step uses disjoint sets of reals, the overall mapping is injective and differentiable, with a Jacobian that can be derived from the densities involved.  
In fact, if the input reals are i.i.d. standard logistic, then the output correlation matrix follows the desired LKJ(\(\eta\)) distribution.

---

## 6. Parameter \(\eta\) and the uniform case

From the construction, the final density of \(\mathbf{R}\) is proportional to \(\det(\mathbf{R})^{\eta-1}\).  
The parameter \(\eta\) enters only through the initial value of \(\beta\) and its decrements:

\[
\beta_{\text{start}} = \eta + \frac{d-2}{2}, \qquad
\beta \leftarrow \beta - \tfrac12 \text{ at each step}.
\]

For \(\eta = 1\) (uniform distribution over correlation matrices) we have:

- Initial \(\beta = 1 + \frac{d-2}{2} = \frac{d}{2}\).
- At step \(k\) (current size \(k\)), \(\beta = \frac{d-k}{2}\).

Thus the Beta shapes for \(y\) are \(\mathrm{Beta}(k/2,\,(d-k)/2)\), and the marginal distribution of any single off‑diagonal correlation \(r_{ij}\) is \(\mathrm{Beta}(d/2,d/2)\) on \((-1,1)\) – a well‑known property of the uniform LKJ distribution.

---

## 7. Comparison with the vine method

Both the C‑vine and the onion method generate the same family of distributions.  
The C‑vine method works with **partial correlations** assigned to a fixed vine structure and then recursively converts them to ordinary correlations.  
The onion method builds the matrix **sequentially** by adding variables, using a spherical construction for the new correlation vector.

**Computational aspects:**  
- The onion method requires a Cholesky decomposition at each step (or an incremental update).  
- The C‑vine method avoids matrix factorisations and can be faster, especially in compiled code.  
- Both methods can be used deterministically by feeding transformed real numbers.

---

## 8. References

- Ghosh, S., & Henderson, S.G. (2003). Behavior of the NORTA method for correlated random vector generation as the dimension increases. *ACM Transactions on Modeling and Computer Simulation*, 13(3), 276–294.
- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. *Journal of Multivariate Analysis*, 100(9), 1989–2001.
- Joe, H. (2006). Generating random correlation matrices based on partial correlations. *Journal of Multivariate Analysis*, 97, 2177–2189.
