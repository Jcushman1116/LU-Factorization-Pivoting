# LU Factorization with Pivoting

## Overview

Implements and validates LU factorization via Gaussian elimination using three
pivoting strategies: no pivoting $A = LU$, partial pivoting $PA = LU$, and complete
pivoting $PAQ = LU$. The algorithm is tested against eight structured matrix types
of increasing dimension, evaluating factorization error, residual error, growth factor,
and condition number to assess numerical stability across methods.

## Methodology

All three factorization methods are implemented within a single function `LUdense`,
controlled by a `routinenum` flag. L and U are stored in-place within the same matrix
for memory efficiency — possible because L is unit diagonal, so the ones and zeros
need not be stored explicitly. Row and column permutations are tracked as index
vectors rather than full permutation matrices.

For each elimination step k, a pivot element is selected and the trailing submatrix
is updated via a rank-one operation. An error guard exits early if the pivot magnitude
falls below $10^{-12}$, protecting against ill-conditioned or singular matrices.

**No pivoting (flag = 1):** Pivot is always $a_{kk}$. Fastest but numerically
fragile. Complexity: $O(n^3)$

**Partial pivoting (flag = 2):** Searches column k for the largest magnitude
element and swaps rows. Reduces risk of small pivot division.
Complexity: $O(n^3) + O(n^2)$ search cost.

**Complete pivoting (flag = 3):** Searches the entire active submatrix $A(k{:}n,\, k{:}n)$
for the largest element and swaps both rows and columns. Most numerically stable,
highest computational cost. Complexity: $O(n^3) + O(n^3)$ search cost.

Accuracy is measured using the 2-norm across three metrics:

$$
E_{\text{fac}} = \frac{\|P_r A P_c - LU\|_2}{\|A\|_2}, \quad
\gamma = \frac{\||L||U|\|_2}{\|A\|_2}, \quad
E_{\text{res}} = \frac{\|b - Ax\|_2}{\|b\|_2}
$$

The condition number $K_2(A) = \|A\|_2 \|A^{-1}\|_2$ is computed for each matrix
type to relate conditioning to factorization stability.

## Matrix Types Tested

- **Diagonal** — near-zero error, growth factor constant at 1 across all methods
- **Anti-diagonal** — no pivoting fails (zero on first diagonal); partial and
  complete pivoting succeed
- **Diagonal + anti-diagonal** — singular structure; error guard triggers for all
  methods, no factorization produced
- **Unit lower triangular** — low error across all methods; growth factor grows
  modestly due to floating point limits
- **Lower triangular** — low error; no pivoting holds growth factor at 1, other
  methods introduce minor growth from unnecessary comparisons
- **Tridiagonal, diagonally dominant** — well-conditioned; low error and growth
  factor of 1 across all methods
- **Growth factor matrix** — growth factor blows up for no pivoting and partial
  pivoting; complete pivoting remains stable, demonstrating the limitation of
  row-only permutations
- **Symmetric positive definite** — stable across all methods; U = DL^T where D
  is the diagonal of U, with $L_{\text{tilde}} = LD^{1/2}$ recovering the
  Cholesky-like structure

## Language

MATLAB

## How to Run

1. Ensure all four `.m` files are in the same directory: `LUdense`,
   `TestDriverProgram3`, `Lvsolve`, `Uvsolve`
2. Run `TestDriverProgram3` — executes all matrix structure tests across all
   three pivoting methods
3. Results are printed as tables showing condition number, mean/max factorization
   error, mean/max growth factor, and mean/max residual error per dimension
4. Full test suite completes in under 30 seconds
