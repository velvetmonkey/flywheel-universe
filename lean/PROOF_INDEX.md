# Proof Index — Budgeted Hebbian Kuramoto Lean Formalisation

Status legend:
- `[OK]` proved, sorry-free
- `[TODO]` open, sorry-bearing stub
- `[HYP]` taken as a Trajectory hypothesis (a *field*, not a *theorem*)

## Algebraic core

| Result | Status | Location |
|---|---|---|
| `cosDiff_symm`, `sinDiff_antisymm`, `cosDiff_diag`, `sinDiff_diag` | `[OK]` | `Defs.lean` |
| `weightGrad_eq_neg_lyapunovGradW` | `[OK]` | `Defs.lean` |
| `phase_contribution_identity` | `[OK]` | `LyapunovDescent.lean` |
| `phase_descent` (`K < 0 ⟹ K · Σ cf² ≤ 0`) | `[OK]` | `LyapunovDescent.lean` |
| `lyapunovGradW_isSymm` (private helper) | `[OK]` | `LyapunovDescent.lean` |
| `sum_full_eq_twice_upperTri` (public, reusable) | `[OK]` | `LyapunovDescent.lean` |

## Constraint set

| Result | Status | Location |
|---|---|---|
| `constraintSet_convex` | `[OK]` | `LyapunovDescent.lean` |
| `constraintSet_isClosed` | `[OK]` | `LyapunovDescent.lean` |
| `zero_mem_constraintSet` | `[OK]` | `LyapunovDescent.lean` |

## Trajectory and weight dynamics

| Result | Status | Notes |
|---|---|---|
| `Trajectory.Wdot`, `Trajectory.hWeight_diff` | `[OK]` (Issue #2) | Pointwise weight derivative + componentwise `HasDerivAt`. |
| `Trajectory.hWeight_dyn` | `[HYP]` (Issue #1) | Variational-inequality hypothesis intended to support projected gradient flow `Ẇ = P_{T_C(W)}(-η∇_W L)`. Does **not** characterise PGF on its own — `Wdot ∈ T_C(W)` is not encoded as a separate usable hypothesis. Combined with `hWeight_diff` + `hW_in_C`, suffices for Issue #3 discharge via Fermat. |
| `Trajectory.Wdot_isSymm` (private helper) | `[OK]` | Derived from `hW_in_C` symmetry + `hWeight_diff` uniqueness. |
| `Trajectory.Wdot_diag_zero` (private helper) | `[OK]` | Derived from `hW_in_C` diagonal-zero + `hWeight_diff` uniqueness. |
| `Trajectory.lyapunovAlong_hasDerivAt` | `[OK]` (Issue #2) | Chain rule on the explicit sum form of `lyapunovFn`; uses `phase_contribution_identity` to fold `−K · Σ W·sin·(cf_i−cf_j)` into `K · Σ cf²`. Sorry-free. |

## Descent and stationarity

| Result | Status | Notes |
|---|---|---|
| `Trajectory.hW_descent_derived` | `[OK]` (Issue #3 closed) | Weight contribution non-positivity proved via Fermat's interior extremum applied to `f(s) := ⟨W(s) − W(t), Ẇ(t) + η∇L(t)⟩`. From `hWeight_dyn` we get `f ≥ 0 ∀s` with `f(t) = 0`, so Fermat gives `f'(t) = 0`, yielding `⟨Ẇ, ∇L⟩_F = −‖Ẇ‖²/η ≤ 0`. Upper-tri restriction via `sum_full_eq_twice_upperTri` + symmetry / zero-diag of `Ẇ` and `∇L`. |
| `lyapunov_descent` | `[OK]` | Now takes no extra hypotheses beyond `Trajectory`. `dL/dt ≤ 0` from `lyapunovAlong_hasDerivAt` + `phase_descent` + `hW_descent_derived`; monotonicity via `antitone_of_deriv_nonpos`. |
| `limit_point_mem_constraintSet` | `[OK]` | Set-membership at a weight limit; companion to the full KKT theorem below. |
| `lyapunovFn_convex_in_W` (private helper) | `[OK]` | Convexity of `L(θ, ·)` in `W`: convexity residual `= (λ/2) Σ (V−W)² ≥ 0`. Used by `limit_point_isKKTStationary`. |
| `limit_point_isKKTStationary` (full KKT) | `[OK]` (Issue #4 closed) | At a weight fixed-point `Wdot t_star = 0`, the state `(phases t_star, weights t_star)` is KKT-stationary: `W_star` globally minimises `L(θ_star, ·)` over `C`. Argument: variational inequality from `hWeight_dyn` at the fixed point (Wdot = 0 collapses to `⟨V − W, η · ∇_W L⟩_F ≥ 0`) → upper-triangle restriction via `sum_full_eq_twice_upperTri` → convexity of `L` in `W` (`lyapunovFn_convex_in_W`) → global optimality. Sorry-free. |

## GitHub issue mapping

| Issue | Title | Status |
|---|---|---|
| #1 | Weight dynamics field missing from `Trajectory` | **Closed** — `hWeight_dyn` + `Wdot` + `hWeight_diff` added. Doc explicitly states this is a variational-inequality *hypothesis*, not a full PGF characterisation. |
| #2 | dL/dt decomposition taken as hypothesis, not derived | **Closed** — `lyapunovAlong_hasDerivAt` proves the chain rule sorry-free. |
| #3 | Weight contribution non-positivity assumed via variational inequality | **Closed** — `hW_descent_derived` proves it sorry-free via Fermat's interior extremum applied to `hWeight_dyn` + `hWeight_diff` + `hW_in_C`. (Surfaced by Gemini during roundtable review; original plan was to leave it open.) |
| #4 | KKT corollary proves set membership only, not full KKT stationarity | **Closed sorry-free** — `limit_point_isKKTStationary` proves full KKT at a weight fixed-point `Wdot t_star = 0` via variational inequality from `hWeight_dyn` + upper-triangle restriction + convexity of `L` in `W`. (Originally planned as a tracked stub; converted to a full proof following the same Fermat-stationary template used for Issue #3.) |

## Caveats and out-of-scope (for this Lean snapshot)

- **Caveat — `hWeight_diff` (two-sided differentiability):**
  The current proof of `hW_descent_derived` (`LyapunovDescent.lean` lines 533–537) uses a Fermat-at-interior-extremum argument that requires `hWeight_diff` to assert everywhere two-sided (`HasDerivAt`) differentiability of the weight trajectory. This is stronger than canonical projected-gradient flow, where velocity can jump at boundary contact with `∂C`.

  The clean resolution requires tangent/normal cone calculus for the specific polyhedral `constraintSet` (row-sum budget + nonneg + support + symmetry + zero diagonal) plus a Moreau decomposition lemma. Neither exists in Mathlib v4.28.0 (the version pinned in `lakefile.toml`) for this constraint shape. Estimated new infrastructure: ~500–1500 lines. Tracked for future work — see `RequestProject/FutureWork.lean`.
- Existence and uniqueness of the projected-gradient ODE solution.
- The Moreau decomposition for tangent / normal cones of `constraintSet`.
- Continuous dependence of `lyapunovGradW` on `(θ, W)` propagated through `Filter.Tendsto` (needed for Issue #4 full-KKT).
- Numerical solver fidelity (handled in `phase2_benchmark.py` validation).
