/-!
# Future Work: Tangent-Cone Formulation of Weight Dynamics

This file sketches the stronger `Trajectory'` structure that would replace
`hWeight_dyn` (variational inequality on `C`) with the canonical projected-
gradient-flow characterisation via tangent/normal cones.

The proof of `hW_descent_derived` currently uses a Fermat argument that
requires everywhere two-sided differentiability (`hWeight_diff`). The cleaner
path replaces this with a one-line Moreau decomposition, but requires
tangent-cone calculus for `constraintSet` — infrastructure not yet in
Mathlib (v4.28.0, the pinned version) for this constraint shape.

## What is needed in Mathlib

1. `tangentCone_constraintSet`: `T_C(W)` for the budget polytope
   (nonneg + symmetric + zero-diagonal + row-sum ≤ B, masked support).
2. `normalCone_constraintSet`: `N_C(W) = (T_C(W))^⊥`.
3. Moreau decomposition: `∀ v, v = proj_{T_C(W)} v + proj_{N_C(W)} v`,
   with the two components orthogonal.

## Aspirational structure (not a proof — do not uncomment without (1)–(3))

```
structure Trajectory' (n : ℕ) (cfg : SystemConfig n) where
  phases  : ℝ → Fin n → ℝ
  weights : ℝ → Matrix (Fin n) (Fin n) ℝ
  Wdot    : ℝ → Matrix (Fin n) (Fin n) ℝ
  hPhase_dyn   : ∀ t i, HasDerivAt (fun s => phases s i)
                   (phaseDot cfg.coupling (weights t) (phases t) i) t
  hW_in_C      : ∀ t, weights t ∈ constraintSet n cfg.support cfg.budget
  hWeight_ac   : -- weights absolutely continuous; Wdot is a.e. velocity
  hWeight_dyn' : ∀ t, Wdot t ∈ tangentCone_constraintSet (weights t) ∧
                   ∀ V ∈ constraintSet n cfg.support cfg.budget,
                     ⟪Wdot t + cfg.learnRate • lyapunovGradW (-1) cfg.decayRate
                       (phases t) (weights t), V - weights t⟫_F ≥ 0
```

With `hWeight_dyn'`, `hW_descent_derived` follows from Moreau (one line):

    ⟨Ẇ, ∇L⟩ = ⟨Ẇ, Ẇ + η∇L⟩ − ‖Ẇ‖² ≤ 0

(since `⟨Ẇ, Ẇ + η∇L⟩ ≤ 0` by tangent-cone VI with `V = W − Ẇ`).

This file is documentation only; it deliberately defines nothing and
imports nothing. Uncommenting the structure above without first developing
items (1)–(3) will not type-check.
-/
