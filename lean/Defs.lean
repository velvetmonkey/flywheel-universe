/-
# Budgeted Hebbian Kuramoto Model — Definitions

This file contains the core definitions for the Budgeted Hebbian Kuramoto model:
- Constraint set C for coupling weights
- Cosine difference matrix
- Lyapunov function L
- Phase and weight dynamics
-/
import Mathlib

noncomputable section

open Matrix Finset Real

namespace KuramotoHebbian

variable {n : ℕ}

/-! ## Index pairs for upper triangular sums -/

/-- The set of strictly upper triangular index pairs {(i,j) : i < j}. -/
def upperTriPairs (n : ℕ) : Finset (Fin n × Fin n) :=
  Finset.univ.filter fun p => p.1 < p.2

/-! ## Cosine and sine difference matrices -/

/-- The cosine difference matrix: `cosDiff θ i j = cos(θ i - θ j)`. -/
def cosDiff (θ : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => Real.cos (θ i - θ j)

/-- The sine difference matrix: `sinDiff θ i j = sin(θ i - θ j)`. -/
def sinDiff (θ : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => Real.sin (θ i - θ j)

/-! ## Constraint set -/

/-- The constraint set C for the coupling weight matrix W:
  - W is symmetric
  - W has nonnegative entries
  - W is supported on A (W_ij = 0 when A_ij = 0)
  - W has zero diagonal
  - Each row sum is bounded by the budget B -/
def constraintSet (n : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) (B : ℝ) :
    Set (Matrix (Fin n) (Fin n) ℝ) :=
  {W | W.IsSymm ∧
       (∀ i j, 0 ≤ W i j) ∧
       (∀ i j, A i j = 0 → W i j = 0) ∧
       (∀ i, W i i = 0) ∧
       (∀ i, ∑ j, W i j ≤ B)}

/-! ## Lyapunov function -/

/-- The Lyapunov function:
  `L(θ, W) = -s · Σ_{i<j} W_ij cos(θ_i - θ_j) + (λ/2) · Σ_{i<j} W_ij²`
  where s = sign(K). -/
def lyapunovFn (n : ℕ) (s lam : ℝ) (θ : Fin n → ℝ) (W : Matrix (Fin n) (Fin n) ℝ) : ℝ :=
  -s * ∑ p ∈ upperTriPairs n, W p.1 p.2 * Real.cos (θ p.1 - θ p.2) +
  lam / 2 * ∑ p ∈ upperTriPairs n, W p.1 p.2 ^ 2

/-! ## Dynamics -/

/-- Kuramoto phase dynamics (zero detuning):
  `dθ_i/dt = K · Σ_j W_ij sin(θ_j - θ_i)`. -/
def phaseDot (K : ℝ) (W : Matrix (Fin n) (Fin n) ℝ) (θ : Fin n → ℝ) (i : Fin n) : ℝ :=
  K * ∑ j, W i j * Real.sin (θ j - θ i)

/-- The Hebbian weight gradient (before projection):
  `G(θ, W) = η · (s · cosDiff(θ) - λ · W)`. -/
def weightGrad (s η lam : ℝ) (θ : Fin n → ℝ) (W : Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  η • (s • cosDiff θ - lam • W)

/-- The gradient of L with respect to W (for the upper triangle entries):
  `∂L/∂W_ij = -s · cos(θ_i - θ_j) + λ · W_ij`. -/
def lyapunovGradW (s lam : ℝ) (θ : Fin n → ℝ) (W : Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => -s * Real.cos (θ i - θ j) + lam * W i j

/-! ## Key algebraic identities -/

/-- The per-node coupling force: `f_i(θ, W) = Σ_j W_ij sin(θ_j - θ_i)`. -/
def couplingForce (W : Matrix (Fin n) (Fin n) ℝ) (θ : Fin n → ℝ) (i : Fin n) : ℝ :=
  ∑ j, W i j * Real.sin (θ j - θ i)

/-! ## Properties of cosDiff and sinDiff -/

theorem cosDiff_symm (θ : Fin n → ℝ) (i j : Fin n) :
    cosDiff θ i j = cosDiff θ j i := by
  simp only [cosDiff, Matrix.of_apply]
  rw [show θ j - θ i = -(θ i - θ j) from by ring, Real.cos_neg]

theorem sinDiff_antisymm (θ : Fin n → ℝ) (i j : Fin n) :
    sinDiff θ i j = -sinDiff θ j i := by
  simp only [sinDiff, Matrix.of_apply]
  rw [show θ i - θ j = -(θ j - θ i) from by ring, Real.sin_neg]

theorem sinDiff_diag (θ : Fin n → ℝ) (i : Fin n) :
    sinDiff θ i i = 0 := by
  simp [sinDiff, Matrix.of_apply, sub_self, Real.sin_zero]

theorem cosDiff_diag (θ : Fin n → ℝ) (i : Fin n) :
    cosDiff θ i i = 1 := by
  simp [cosDiff, Matrix.of_apply, sub_self, Real.cos_zero]

/-! ## The weight gradient is the negative gradient of L -/

/-- The weight gradient direction `s · cosDiff - λW` equals `-∇_W L`. -/
theorem weightGrad_eq_neg_lyapunovGradW (s lam : ℝ) (θ : Fin n → ℝ)
    (W : Matrix (Fin n) (Fin n) ℝ) :
    s • cosDiff θ - lam • W = -lyapunovGradW s lam θ W := by
  ext i j
  simp [cosDiff, lyapunovGradW, Matrix.of_apply, Matrix.sub_apply, Matrix.smul_apply,
        Pi.smul_apply, smul_eq_mul, Matrix.neg_apply]
  ring

end KuramotoHebbian
