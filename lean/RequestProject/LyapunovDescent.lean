/-
# Lyapunov Descent for Budgeted Hebbian Kuramoto

This file formalises the main results:

1. The constraint set C is convex and closed.
2. The key algebraic identity: the phase contribution to dL/dt equals
   `K · Σᵢ fᵢ²` where `fᵢ = Σⱼ Wᵢⱼ sin(θⱼ - θᵢ)`, which is ≤ 0 when K < 0.
3. The weight update direction is the negative gradient of L projected onto C.
4. The Lyapunov descent theorem: dL/dt ≤ 0.
5. The KKT stationarity corollary for limit points.
-/

import Mathlib
import RequestProject.Defs

noncomputable section

open Matrix Finset Real KuramotoHebbian

namespace KuramotoHebbian

variable {n : ℕ}

/-! ## Convexity and closedness of the constraint set -/

/-- The constraint set C is convex. -/
theorem constraintSet_convex (A : Matrix (Fin n) (Fin n) ℝ) {B : ℝ} (hB : 0 < B) :
    Convex ℝ (constraintSet n A B) := by
  intro W hW W' hW' a b ha hb hab
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · simp_all +decide [Matrix.IsSymm]
    exact congr_arg₂ _ (congr_arg _ hW.1) (congr_arg _ hW'.1)
  · exact fun i j =>
      add_nonneg (mul_nonneg ha (hW.2.1 i j)) (mul_nonneg hb (hW'.2.1 i j))
  · intro i j hij
    have := hW.2.2.1 i j hij; have := hW'.2.2.1 i j hij
    simp_all +decide [Matrix.add_apply, Matrix.smul_apply]
  · exact fun i => by simp +decide [hW.2.2.2.1, hW'.2.2.2.1]
  · intro i
    have := hW.2.2.2.2 i; have := hW'.2.2.2.2 i
    simp_all +decide [Finset.sum_add_distrib]
    simpa only [← Finset.mul_sum _ _ _, ← Finset.sum_add_distrib] using by nlinarith

/-- The constraint set C is closed. -/
theorem constraintSet_isClosed (A : Matrix (Fin n) (Fin n) ℝ) (B : ℝ) :
    IsClosed (constraintSet n A B) := by
  unfold constraintSet
  simp +decide only [Set.setOf_and, Set.setOf_forall]
  refine IsClosed.inter ?_ ?_
  · exact isClosed_eq continuous_id.matrix_transpose continuous_id
  · refine IsClosed.inter ?_ ?_
    · exact isClosed_iInter fun i => isClosed_iInter fun j =>
        isClosed_le continuous_const <|
          continuous_apply _ |> Continuous.comp <| continuous_apply _
    · refine IsClosed.inter ?_ ?_
      · exact isClosed_iInter fun i => isClosed_iInter fun j =>
          isClosed_iInter fun _ =>
            isClosed_eq (continuous_apply _ |> Continuous.comp <| continuous_apply _)
              continuous_const
      · refine IsClosed.inter ?_ ?_
        · exact isClosed_iInter fun i =>
            isClosed_eq (continuous_apply _ |> Continuous.comp <| continuous_apply _)
              continuous_const
        · exact isClosed_iInter fun i =>
            isClosed_le (continuous_finset_sum _ fun j _ =>
              continuous_apply _ |> Continuous.comp <| continuous_apply _) continuous_const

/-- The zero matrix is in C (so C is nonempty when B ≥ 0). -/
theorem zero_mem_constraintSet (A : Matrix (Fin n) (Fin n) ℝ) {B : ℝ} (hB : 0 ≤ B) :
    (0 : Matrix (Fin n) (Fin n) ℝ) ∈ constraintSet n A B := by
  constructor <;> aesop

/-! ## Key algebraic identity for phase contribution -/

/-- When W is symmetric, the sum `Σ_{i<j} Wᵢⱼ sin(θᵢ-θⱼ)(fᵢ-fⱼ)` where
  `fᵢ = Σⱼ Wᵢⱼ sin(θⱼ-θᵢ)` equals `-Σᵢ fᵢ²`.

  The proof uses symmetry of W and antisymmetry of sin to convert the upper-triangle
  sum to half the full sum, then expands using `Σⱼ Wᵢⱼ sin(θᵢ-θⱼ) = -fᵢ`. -/
theorem phase_contribution_identity
    (W : Matrix (Fin n) (Fin n) ℝ) (theta : Fin n → ℝ)
    (hW_symm : W.IsSymm) :
    ∑ p ∈ upperTriPairs n,
      W p.1 p.2 * Real.sin (theta p.1 - theta p.2) *
        (couplingForce W theta p.1 - couplingForce W theta p.2) =
    -∑ i : Fin n, couplingForce W theta i ^ 2 := by
  have h_full : ∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
      (couplingForce W theta i - couplingForce W theta j) =
      -2 * ∑ i, (couplingForce W theta i) ^ 2 := by
    have hA : ∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
        couplingForce W theta i = -∑ i, couplingForce W theta i ^ 2 := by
      have : ∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
          couplingForce W theta i =
          ∑ i, (∑ j, W i j * Real.sin (theta i - theta j)) *
            couplingForce W theta i := by
        simp +decide only [Finset.sum_mul _ _ _]
      rw [this, ← Finset.sum_neg_distrib]
      exact Finset.sum_congr rfl fun i _ => by
        rw [show ∑ j, W i j * Real.sin (theta i - theta j) =
            -couplingForce W theta i from by
          rw [show couplingForce W theta i =
              ∑ j, W i j * Real.sin (theta j - theta i) from rfl]
          rw [← Finset.sum_neg_distrib]
          exact Finset.sum_congr rfl fun j _ => by
            rw [← mul_neg, ← Real.sin_neg]; ring]
        ring
    have hB : ∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
        couplingForce W theta j =
        -∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
          couplingForce W theta i := by
      simp +decide only [← sum_neg_distrib]
      rw [Finset.sum_comm]; congr; ext i; congr; ext j; ring
      rw [← hW_symm.apply]
      rw [show -theta j + theta i = -(theta j - theta i) by ring]
      rw [Real.sin_neg]; ring
    simp_all +decide [mul_sub]; linarith
  have h_split : ∑ i, ∑ j, W i j * Real.sin (theta i - theta j) *
      (couplingForce W theta i - couplingForce W theta j) =
      ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 < p.2),
        W p.1 p.2 * Real.sin (theta p.1 - theta p.2) *
          (couplingForce W theta p.1 - couplingForce W theta p.2) +
      ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 > p.2),
        W p.1 p.2 * Real.sin (theta p.1 - theta p.2) *
          (couplingForce W theta p.1 - couplingForce W theta p.2) := by
    rw [← Finset.sum_union]
    · rw [← Finset.sum_product', ← Finset.sum_subset]
      · grind
      · grind
    · exact Finset.disjoint_filter.mpr fun _ _ _ _ => lt_asymm ‹_› ‹_›
  have h_swap : ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 > p.2),
      W p.1 p.2 * Real.sin (theta p.1 - theta p.2) *
        (couplingForce W theta p.1 - couplingForce W theta p.2) =
      ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 < p.2),
        W p.2 p.1 * Real.sin (theta p.2 - theta p.1) *
          (couplingForce W theta p.2 - couplingForce W theta p.1) := by
    apply Finset.sum_bij (fun p _ => (p.2, p.1))
    · grind +extAll
    · grind
    · exact fun p hp => ⟨(p.2, p.1), by simpa using hp, rfl⟩
    · grind
  have h_swap_sign : ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 < p.2),
      W p.2 p.1 * Real.sin (theta p.2 - theta p.1) *
        (couplingForce W theta p.2 - couplingForce W theta p.1) =
      ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.1 < p.2),
        W p.1 p.2 * Real.sin (theta p.1 - theta p.2) *
          (couplingForce W theta p.1 - couplingForce W theta p.2) := by
    exact Finset.sum_congr rfl fun p _ => by
      rw [← hW_symm.apply, ← neg_sub, Real.sin_neg]; ring
  linarith!

/-- The phase contribution to dL/dt is non-positive when K < 0. -/
theorem phase_descent {K : ℝ} (hK : K < 0)
    (W : Matrix (Fin n) (Fin n) ℝ) (theta : Fin n → ℝ) :
    K * ∑ i : Fin n, couplingForce W theta i ^ 2 ≤ 0 := by
  exact mul_nonpos_of_nonpos_of_nonneg hK.le <|
    Finset.sum_nonneg fun _ _ => sq_nonneg _

/-! ## Weight contribution: projected gradient descent -/

/-- The Frobenius projection onto a closed convex set. We take this as an abstract
  projection operator satisfying the variational inequality. -/
structure FrobeniusProjection (n : ℕ) (C : Set (Matrix (Fin n) (Fin n) ℝ)) where
  proj : Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ
  proj_mem : ∀ x, proj x ∈ C
  variational : ∀ x y, y ∈ C →
    ∑ i : Fin n, ∑ j : Fin n, (proj x i j - x i j) * (y i j - proj x i j) ≥ 0

/-! ## System configuration and trajectories -/

/-- The full system configuration. -/
structure SystemConfig (n : ℕ) where
  /-- Support mask A, symmetric with zero diagonal -/
  support : Matrix (Fin n) (Fin n) ℝ
  /-- Budget B > 0 -/
  budget : ℝ
  /-- Coupling constant K < 0 -/
  coupling : ℝ
  /-- Learning rate, positive -/
  learnRate : ℝ
  /-- Decay parameter, positive -/
  decayRate : ℝ
  hK_neg : coupling < 0
  heta_pos : 0 < learnRate
  hlam_pos : 0 < decayRate
  hB_pos : 0 < budget
  hA_symm : support.IsSymm
  hA_binary : ∀ i j, support i j = 0 ∨ support i j = 1
  hA_diag : ∀ i, support i i = 0

/-- A trajectory of the Kuramoto-Hebbian system.

  Encodes phase and weight evolution together with explicit differentiability
  data. The weight dynamics is specified by a variational-inequality
  hypothesis (`hWeight_dyn`) intended to model projected gradient flow on `L`
  over the constraint set `C`; see that field's docstring for the precise
  semantics and the load-bearing caveats. -/
structure Trajectory (n : ℕ) (cfg : SystemConfig n) where
  /-- Phase trajectory θ(t) -/
  phases : ℝ → Fin n → ℝ
  /-- Weight trajectory W(t) -/
  weights : ℝ → Matrix (Fin n) (Fin n) ℝ
  /-- Pointwise time derivative of the weight matrix (Issue #1, #2). -/
  Wdot : ℝ → Matrix (Fin n) (Fin n) ℝ
  /-- Weights stay in C -/
  hW_in_C : ∀ t, weights t ∈ constraintSet n cfg.support cfg.budget
  /-- Phase dynamics: dθᵢ/dt = K · Σⱼ Wᵢⱼ sin(θⱼ - θᵢ) -/
  hPhase_dyn : ∀ t i,
    HasDerivAt (fun s => phases s i)
      (phaseDot cfg.coupling (weights t) (phases t) i) t
  /-- (Issue #2) Each weight entry is differentiable in time, with derivative
    `Wdot t i j`. Required for the chain-rule derivation of `dL/dt`. -/
  hWeight_diff : ∀ t i j,
    HasDerivAt (fun s => weights s i j) (Wdot t i j) t
  /-- (Issue #1) Variational-inequality hypothesis intended to support projected
    gradient flow `Ẇ = P_{T_C(W)}(-η ∇_W L)` (Option B in the Lean-issues design
    note).

    For every `V ∈ C` and every time `t`,
      `0 ≤ ∑_{i,j} (V_ij - W(t)_ij) · (Ẇ_ij(t) + η · ∂L/∂W_ij(θ(t),W(t)))`.

    Semantics. This is the normal-cone-membership condition
      `Ẇ + η ∇_W L ∈ N_C(W(t))`
    of projected-gradient flow. It does **not** fully characterise the
    projected-gradient ODE on its own: tangent feasibility `Ẇ(t) ∈ T_C(W(t))`
    is not encoded here as a separate, usable hypothesis. Combined with
    `hWeight_diff` (everywhere two-sided differentiability) and `hW_in_C`,
    this field is nevertheless strong enough to discharge `hW_descent` via
    Fermat's interior extremum (see `Trajectory.hW_descent_derived`).

    Caveat: requiring everywhere two-sided `hWeight_diff` is *stronger* than
    the canonical projected-gradient ODE, whose velocity can jump on contact
    with `∂C`. This restricts `Trajectory` to curves with no such transversal
    boundary contact — adequate for the present descent argument but not a
    full formalisation of measure-theoretic projected flow. -/
  hWeight_dyn : ∀ t, ∀ V ∈ constraintSet n cfg.support cfg.budget,
    0 ≤ ∑ i : Fin n, ∑ j : Fin n,
      (V i j - weights t i j) *
      (Wdot t i j +
       cfg.learnRate * lyapunovGradW (-1) cfg.decayRate (phases t) (weights t) i j)

/-- The Lyapunov function evaluated along a trajectory (using s = -1). -/
def Trajectory.lyapunovAlong {cfg : SystemConfig n}
    (traj : Trajectory n cfg) (t : ℝ) : ℝ :=
  lyapunovFn n (-1) cfg.decayRate (traj.phases t) (traj.weights t)

/-! ## The derivative of L along the trajectory

The time derivative of L splits into two terms:
  dL/dt = (phase contribution) + (weight contribution)

Phase contribution: From the algebraic identity `phase_contribution_identity`,
  this equals `K · Σᵢ fᵢ² ≤ 0` (since K < 0).

Weight contribution: From projected gradient descent, this is ≤ 0.

We express this decomposition as follows.
-/

/-- The phase contribution to dL/dt. With s = -1 and the Kuramoto dynamics,
  this equals `K · Σᵢ (Σⱼ Wᵢⱼ sin(θⱼ-θᵢ))²`, which is ≤ 0 when K < 0. -/
def phaseContribution (K : ℝ) (W : Matrix (Fin n) (Fin n) ℝ) (theta : Fin n → ℝ) : ℝ :=
  K * ∑ i : Fin n, couplingForce W theta i ^ 2

/-- The weight contribution to dL/dt. When the weight dynamics is projected
  gradient descent on L(θ,·), this is non-positive. -/
def weightContribution (lam : ℝ) (theta : Fin n → ℝ)
    (W : Matrix (Fin n) (Fin n) ℝ) (dW : Matrix (Fin n) (Fin n) ℝ) : ℝ :=
  ∑ p ∈ upperTriPairs n,
    dW p.1 p.2 * (-(-1) * Real.cos (theta p.1 - theta p.2) + lam * W p.1 p.2)

/-! ### Chain-rule decomposition of `dL/dt` (Issue #2)

  The Lyapunov function evaluated along a trajectory is differentiable in time,
  and its derivative decomposes into the phase and weight contributions defined
  above. This is the chain-rule application underlying the `hL_deriv` hypothesis
  that was previously taken as an axiom in `lyapunov_descent`. -/

/-- (Issue #2) The Lyapunov function evaluated along a trajectory has the
  expected time derivative, given by the sum of the phase and weight
  contributions.

  Proven by the chain rule on the explicit sum form of `lyapunovFn` (with
  `s = -1`), combining `HasDerivAt` for each weight entry (from
  `hWeight_diff`), for each phase (from `hPhase_dyn`), and `Real.cos`,
  together with `phase_contribution_identity` (already proved) to identify
  the phase part with `K · Σᵢ fᵢ²`. -/
theorem Trajectory.lyapunovAlong_hasDerivAt
    {cfg : SystemConfig n} (traj : Trajectory n cfg) (t : ℝ) :
    HasDerivAt traj.lyapunovAlong
      (phaseContribution cfg.coupling (traj.weights t) (traj.phases t) +
       weightContribution cfg.decayRate (traj.phases t) (traj.weights t)
         (traj.Wdot t)) t := by
  -- Per-pair HasDerivAt for `W_ij(s) * cos(θ_i(s) - θ_j(s))`
  have h_cos_pair : ∀ p ∈ upperTriPairs n,
      HasDerivAt (fun s => traj.weights s p.1 p.2 *
                            Real.cos (traj.phases s p.1 - traj.phases s p.2))
        (traj.Wdot t p.1 p.2 *
            Real.cos (traj.phases t p.1 - traj.phases t p.2) +
         traj.weights t p.1 p.2 *
           (-Real.sin (traj.phases t p.1 - traj.phases t p.2) *
            (phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.1 -
             phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.2))) t := by
    intro p _
    have hW := traj.hWeight_diff t p.1 p.2
    have hθ : HasDerivAt (fun s => traj.phases s p.1 - traj.phases s p.2)
        (phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.1 -
         phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.2) t :=
      (traj.hPhase_dyn t p.1).sub (traj.hPhase_dyn t p.2)
    exact hW.mul hθ.cos
  -- Per-pair HasDerivAt for `(W_ij(s))^2`
  have h_sq_pair : ∀ p ∈ upperTriPairs n,
      HasDerivAt (fun s => (traj.weights s p.1 p.2) ^ 2)
        (2 * traj.weights t p.1 p.2 * traj.Wdot t p.1 p.2) t := by
    intro p _
    have hW := traj.hWeight_diff t p.1 p.2
    have hPow := hW.pow 2
    convert hPow using 1
    push_cast
    ring
  -- Sum the per-pair derivatives.
  have h_cos_sum := HasDerivAt.fun_sum h_cos_pair
  have h_sq_sum := HasDerivAt.fun_sum h_sq_pair
  -- Assemble the HasDerivAt for the unfolded Lyapunov function:
  --   `(-(-1) * cos_sum) + (lam/2) * sq_sum`.
  have h_full : HasDerivAt
      (fun s => -(-1 : ℝ) *
          (∑ p ∈ upperTriPairs n, traj.weights s p.1 p.2 *
            Real.cos (traj.phases s p.1 - traj.phases s p.2))
        + cfg.decayRate / 2 *
          (∑ p ∈ upperTriPairs n, (traj.weights s p.1 p.2) ^ 2))
      (-(-1 : ℝ) *
          (∑ p ∈ upperTriPairs n,
            (traj.Wdot t p.1 p.2 *
              Real.cos (traj.phases t p.1 - traj.phases t p.2) +
             traj.weights t p.1 p.2 *
               (-Real.sin (traj.phases t p.1 - traj.phases t p.2) *
                (phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.1 -
                 phaseDot cfg.coupling (traj.weights t) (traj.phases t) p.2))))
        + cfg.decayRate / 2 *
          (∑ p ∈ upperTriPairs n, 2 * traj.weights t p.1 p.2 * traj.Wdot t p.1 p.2)) t :=
    (h_cos_sum.const_mul _).add (h_sq_sum.const_mul _)
  -- `traj.lyapunovAlong` definitionally equals the unfolded sum form.
  have hfn_eq : traj.lyapunovAlong = fun s =>
      -(-1 : ℝ) *
        (∑ p ∈ upperTriPairs n, traj.weights s p.1 p.2 *
          Real.cos (traj.phases s p.1 - traj.phases s p.2))
      + cfg.decayRate / 2 *
        (∑ p ∈ upperTriPairs n, (traj.weights s p.1 p.2) ^ 2) := rfl
  rw [hfn_eq]
  -- Identify the chain-rule output with `phaseContribution + weightContribution`.
  have hW_symm : (traj.weights t).IsSymm := (traj.hW_in_C t).1
  have hPCI := phase_contribution_identity (traj.weights t) (traj.phases t) hW_symm
  -- `phaseDot K W θ i = K * couplingForce W θ i` by definition.
  have hpd : ∀ i, phaseDot cfg.coupling (traj.weights t) (traj.phases t) i
                = cfg.coupling * couplingForce (traj.weights t) (traj.phases t) i :=
    fun _ => rfl
  convert h_full using 1
  -- Goal: `phaseContribution + weightContribution = chain_rule_output`.
  simp only [phaseContribution, weightContribution, hpd]
  -- Step A: split the chain-rule sum (RHS) per-pair into a "weight part" and
  -- a "phase part" using `ring` at the per-term level, then expose both as
  -- separate sums.
  have h_split :
      -(-1 : ℝ) *
          (∑ p ∈ upperTriPairs n,
            (traj.Wdot t p.1 p.2 *
                Real.cos (traj.phases t p.1 - traj.phases t p.2) +
              traj.weights t p.1 p.2 *
                (-Real.sin (traj.phases t p.1 - traj.phases t p.2) *
                  (cfg.coupling * couplingForce (traj.weights t) (traj.phases t) p.1 -
                   cfg.coupling * couplingForce (traj.weights t) (traj.phases t) p.2))))
        + cfg.decayRate / 2 *
          (∑ p ∈ upperTriPairs n, 2 * traj.weights t p.1 p.2 * traj.Wdot t p.1 p.2) =
      (∑ p ∈ upperTriPairs n,
          traj.Wdot t p.1 p.2 *
            (-(-1 : ℝ) * Real.cos (traj.phases t p.1 - traj.phases t p.2) +
             cfg.decayRate * traj.weights t p.1 p.2))
      + (-cfg.coupling) *
          (∑ p ∈ upperTriPairs n,
            traj.weights t p.1 p.2 *
              Real.sin (traj.phases t p.1 - traj.phases t p.2) *
              (couplingForce (traj.weights t) (traj.phases t) p.1 -
               couplingForce (traj.weights t) (traj.phases t) p.2)) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib,
        Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun p _ => ?_)
    ring
  rw [h_split, hPCI]
  ring

/-! ### Issue #3 discharge: weight contribution non-positivity via Fermat -/

/-- The time derivative of the weight trajectory is symmetric.
  Follows from pointwise symmetry of `weights s` (via `hW_in_C s`) by
  uniqueness of derivatives. -/
private lemma Trajectory.Wdot_isSymm {cfg : SystemConfig n}
    (traj : Trajectory n cfg) (t : ℝ) : (traj.Wdot t).IsSymm := by
  ext i j
  have hij := traj.hWeight_diff t i j
  have hji := traj.hWeight_diff t j i
  -- pointwise: weights s j i = weights s i j by symmetry
  have h_eq : (fun s => traj.weights s j i) = (fun s => traj.weights s i j) := by
    funext s
    exact (traj.hW_in_C s).1.apply i j
  rw [h_eq] at hji
  -- both hij and hji are HasDerivAt of the same function fun s => weights s i j;
  -- uniqueness of the derivative gives Wd t i j = Wd t j i
  exact (hij.unique hji).symm

/-- The time derivative of the weight trajectory has zero diagonal,
  inherited from `weights s i i = 0` for all `s`. -/
private lemma Trajectory.Wdot_diag_zero {cfg : SystemConfig n}
    (traj : Trajectory n cfg) (t : ℝ) (i : Fin n) : traj.Wdot t i i = 0 := by
  have hii := traj.hWeight_diff t i i
  have h_const : (fun s => traj.weights s i i) = (fun _ : ℝ => (0 : ℝ)) := by
    funext s
    exact (traj.hW_in_C s).2.2.2.1 i
  rw [h_const] at hii
  exact hii.unique (hasDerivAt_const t 0)

/-- `lyapunovGradW (-1) lam θ W` is symmetric in its index arguments when `W`
  is symmetric (since `cos(θ_i - θ_j) = cos(θ_j - θ_i)`). -/
private lemma lyapunovGradW_isSymm
    (lam : ℝ) (theta : Fin n → ℝ) {W : Matrix (Fin n) (Fin n) ℝ}
    (hW : W.IsSymm) :
    (lyapunovGradW (-1 : ℝ) lam theta W).IsSymm := by
  ext i j
  simp only [Matrix.transpose_apply, lyapunovGradW, Matrix.of_apply]
  rw [show theta j - theta i = -(theta i - theta j) from by ring, Real.cos_neg,
      hW.apply i j]

/-- For symmetric matrices `A` and `B` over `Fin n × Fin n` with `A` zero-
  diagonal, the full Frobenius inner product equals twice the upper-triangle
  pairing. Generic helper, made `public` because it is reusable beyond this
  file. -/
lemma sum_full_eq_twice_upperTri
    {A B : Matrix (Fin n) (Fin n) ℝ}
    (hAsymm : A.IsSymm) (hBsymm : B.IsSymm) (hAdiag : ∀ i, A i i = 0) :
    ∑ i, ∑ j, A i j * B i j =
    2 * ∑ p ∈ upperTriPairs n, A p.1 p.2 * B p.1 p.2 := by
  -- Step 1: split into upper + lower (diagonal contributes zero since A i i = 0).
  -- We mirror the pattern from `phase_contribution_identity`.
  have h_split : ∑ i, ∑ j, A i j * B i j =
      ∑ p ∈ upperTriPairs n, A p.1 p.2 * B p.1 p.2 +
      ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.2 < p.1),
        A p.1 p.2 * B p.1 p.2 := by
    rw [← Finset.sum_union]
    · rw [← Finset.sum_product', ← Finset.sum_subset]
      · -- inclusion: (upper ∪ lower) ⊆ univ
        intros p _
        exact Finset.mem_univ p
      · -- complement: the diagonal. A i i = 0 kills the term.
        intros p _hp hp'
        simp only [upperTriPairs, Finset.mem_union, Finset.mem_filter,
                   Finset.mem_univ, true_and, not_or, not_lt] at hp'
        have heq : p.1 = p.2 := le_antisymm hp'.2 hp'.1
        rw [heq, hAdiag, zero_mul]
    · -- disjointness of upper and lower
      rw [Finset.disjoint_left]
      intros p hp hp'
      simp only [upperTriPairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp hp'
      omega
  rw [h_split]
  -- Step 2: lower-triangle sum equals upper-triangle sum via swap + symmetry.
  have h_lower_eq : ∑ p ∈ Finset.univ.filter (fun p : Fin n × Fin n => p.2 < p.1),
      A p.1 p.2 * B p.1 p.2 =
      ∑ p ∈ upperTriPairs n, A p.1 p.2 * B p.1 p.2 := by
    refine Finset.sum_bij (fun p _ => (p.2, p.1)) ?_ ?_ ?_ ?_
    · -- maps into upperTriPairs
      intros p hp
      simp only [upperTriPairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp ⊢
      exact hp
    · -- injective on the filtered set
      intros p₁ _ p₂ _ h
      have hpair : p₁.2 = p₂.2 ∧ p₁.1 = p₂.1 := Prod.mk.inj h
      exact Prod.ext hpair.2 hpair.1
    · -- surjective onto upperTriPairs
      intros q hq
      simp only [upperTriPairs, Finset.mem_filter, Finset.mem_univ, true_and] at hq
      refine ⟨(q.2, q.1), ?_, rfl⟩
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      exact hq
    · -- preserves value via symmetry
      intros p _
      rw [hAsymm.apply p.1 p.2, hBsymm.apply p.1 p.2]
  rw [h_lower_eq]
  ring

/-- (Issue #3 — discharged by Fermat's interior extremum.)
  The weight contribution to `dL/dt` is non-positive along any trajectory.

  **Argument.** Define `f(s) := ∑_{i,j} (W(s)_ij − W(t)_ij) · (Ẇ(t)_ij + η · ∂L/∂W_ij(t))`.
  By `hWeight_dyn` applied at time `t` with `V := W(s)` (admissible because
  `hW_in_C s`), `f(s) ≥ 0` for all `s`. At `s = t`, `f(t) = 0`, so `t` is a
  global minimum of `f`. `f` is differentiable at `t` (via `hWeight_diff`)
  with derivative `∑_{i,j} Ẇ(t)_ij · (Ẇ(t)_ij + η · ∂L/∂W_ij(t))`. By Fermat's
  interior-extremum theorem this derivative vanishes:
    `‖Ẇ(t)‖²_F + η · ⟨Ẇ(t), ∇_W L(t)⟩_F = 0`.
  Since `η > 0` and `‖Ẇ‖²_F ≥ 0`, this forces `⟨Ẇ, ∇_W L⟩_F ≤ 0`. Symmetry
  and zero-diagonal of `Ẇ(t)`, together with symmetry of `∇_W L(t)`, give
  `⟨Ẇ, ∇_W L⟩_F = 2 · weightContribution`, so `weightContribution ≤ 0`. -/
theorem Trajectory.hW_descent_derived {cfg : SystemConfig n}
    (traj : Trajectory n cfg) (t : ℝ) :
    weightContribution cfg.decayRate (traj.phases t) (traj.weights t)
      (traj.Wdot t) ≤ 0 := by
  -- Abbreviations
  set Wd := traj.Wdot t with hWd_def
  set W := traj.weights t with hW_def
  set θ := traj.phases t with hθ_def
  set gL := lyapunovGradW (-1 : ℝ) cfg.decayRate θ W with hgL_def
  have heta : 0 < cfg.learnRate := cfg.heta_pos
  -- Test function f(s) := ⟨W(s) - W(t), Wd + eta * gL⟩_full
  let f : ℝ → ℝ := fun s =>
    ∑ i, ∑ j, (traj.weights s i j - W i j) * (Wd i j + cfg.learnRate * gL i j)
  -- f is differentiable at t with the expected derivative
  have hf_deriv : HasDerivAt f
      (∑ i, ∑ j, Wd i j * (Wd i j + cfg.learnRate * gL i j)) t := by
    refine HasDerivAt.fun_sum (fun i _ => ?_)
    refine HasDerivAt.fun_sum (fun j _ => ?_)
    exact ((traj.hWeight_diff t i j).sub_const _).mul_const _
  -- f(t) = 0
  have hf_t : f t = 0 := by
    show (∑ i, ∑ j, (traj.weights t i j - W i j) *
            (Wd i j + cfg.learnRate * gL i j)) = 0
    refine Finset.sum_eq_zero fun i _ => ?_
    refine Finset.sum_eq_zero fun j _ => ?_
    simp [hW_def]
  -- f(s) ≥ 0 for all s — direct from hWeight_dyn at time t with V := traj.weights s
  have hf_nonneg : ∀ s, 0 ≤ f s := fun s =>
    traj.hWeight_dyn t (traj.weights s) (traj.hW_in_C s)
  -- t is a local (in fact global) minimum of f
  have hf_min : IsLocalMin f t :=
    Filter.Eventually.of_forall (fun s => hf_t ▸ hf_nonneg s)
  -- Fermat: f'(t) = 0
  have hf'_zero : (∑ i, ∑ j, Wd i j * (Wd i j + cfg.learnRate * gL i j)) = 0 :=
    hf_min.hasDerivAt_eq_zero hf_deriv
  -- Algebra: ∑ Wd² + η · ⟨Wd, gL⟩_F = 0 ⟹ ⟨Wd, gL⟩_F ≤ 0
  have h_split :
      (∑ i, ∑ j, Wd i j * (Wd i j + cfg.learnRate * gL i j)) =
      (∑ i, ∑ j, Wd i j ^ 2) + cfg.learnRate * (∑ i, ∑ j, Wd i j * gL i j) := by
    rw [show (fun i : Fin n => ∑ j, Wd i j * (Wd i j + cfg.learnRate * gL i j)) =
        (fun i => (∑ j, Wd i j ^ 2) + cfg.learnRate * (∑ j, Wd i j * gL i j)) from ?_]
    · rw [Finset.sum_add_distrib, ← Finset.mul_sum]
    · funext i
      rw [show (fun j : Fin n => Wd i j * (Wd i j + cfg.learnRate * gL i j)) =
          (fun j => Wd i j ^ 2 + cfg.learnRate * (Wd i j * gL i j)) from
        funext fun j => by ring]
      rw [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [h_split] at hf'_zero
  have h_sq_nn : 0 ≤ ∑ i, ∑ j, Wd i j ^ 2 :=
    Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => sq_nonneg _
  have h_full_nonpos : (∑ i, ∑ j, Wd i j * gL i j) ≤ 0 := by
    nlinarith
  -- Convert full sum to upper-triangle: ⟨Wd, gL⟩_F = 2 · weightContribution
  have hWd_symm : Wd.IsSymm := traj.Wdot_isSymm t
  have hWd_diag : ∀ i, Wd i i = 0 := traj.Wdot_diag_zero t
  have hW_symm_t : W.IsSymm := (traj.hW_in_C t).1
  have hgL_symm : gL.IsSymm := lyapunovGradW_isSymm cfg.decayRate θ hW_symm_t
  have h_double : (∑ i, ∑ j, Wd i j * gL i j) =
      2 * ∑ p ∈ upperTriPairs n, Wd p.1 p.2 * gL p.1 p.2 :=
    sum_full_eq_twice_upperTri hWd_symm hgL_symm hWd_diag
  rw [h_double] at h_full_nonpos
  -- Identify ∑_{i<j} Wd · gL with weightContribution
  have h_id : (∑ p ∈ upperTriPairs n, Wd p.1 p.2 * gL p.1 p.2) =
      weightContribution cfg.decayRate θ W Wd := by
    show (∑ p ∈ upperTriPairs n, Wd p.1 p.2 * gL p.1 p.2) =
        ∑ p ∈ upperTriPairs n,
          Wd p.1 p.2 * (-(-1) * Real.cos (θ p.1 - θ p.2) + cfg.decayRate * W p.1 p.2)
    refine Finset.sum_congr rfl fun p _ => ?_
    simp only [hgL_def, lyapunovGradW, Matrix.of_apply]
  rw [h_id] at h_full_nonpos
  linarith

/-! ### Main descent theorem -/

/-
**Main Theorem (Lyapunov Descent)**: The Lyapunov function `L` is non-increasing
  along trajectories of the Budgeted Hebbian Kuramoto system.

  The chain-rule decomposition `lyapunovAlong_hasDerivAt` (Issue #2) gives
  `dL/dt = phaseContribution + weightContribution`. The phase term is `≤ 0`
  by `phase_descent`. The weight term is `≤ 0` by `hW_descent_derived`
  (Issue #3 discharged via Fermat's interior extremum applied to
  `hWeight_dyn`).
-/
theorem lyapunov_descent (cfg : SystemConfig n) (traj : Trajectory n cfg)
    (t₁ t₂ : ℝ) (h12 : t₁ ≤ t₂) :
    traj.lyapunovAlong t₂ ≤ traj.lyapunovAlong t₁ := by
  have h_diff : Differentiable ℝ traj.lyapunovAlong := fun s =>
    (traj.lyapunovAlong_hasDerivAt s).differentiableAt
  have h_deriv_nonpos : ∀ s, deriv traj.lyapunovAlong s ≤ 0 := fun s => by
    rw [(traj.lyapunovAlong_hasDerivAt s).deriv]
    exact add_nonpos (phase_descent cfg.hK_neg _ _) (traj.hW_descent_derived s)
  exact antitone_of_deriv_nonpos h_diff h_deriv_nonpos h12

/-! ## KKT Stationarity Corollary (Issue #4 — closed sorry-free)

  Both halves proved: `limit_point_mem_constraintSet` for primal feasibility,
  `limit_point_isKKTStationary` for full global optimality at a weight
  fixed-point (`Wdot t_star = 0`). The latter uses the variational inequality
  `hWeight_dyn`, the new `lyapunovFn_convex_in_W` convexity lemma, and
  `sum_full_eq_twice_upperTri` from the Issue #3 discharge. -/

/-- A state (θ, W) is KKT-stationary for the weight optimisation if W minimises
  L(θ, ·) over C. -/
def IsKKTStationary (n : ℕ) (cfg : SystemConfig n)
    (theta : Fin n → ℝ) (W : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  W ∈ constraintSet n cfg.support cfg.budget ∧
  ∀ W' ∈ constraintSet n cfg.support cfg.budget,
    lyapunovFn n (-1) cfg.decayRate theta W ≤ lyapunovFn n (-1) cfg.decayRate theta W'

/--
**Set-membership corollary.** Any limit point `W*` of the weight trajectory
belongs to the constraint set `C`. Companion to `limit_point_isKKTStationary`
below, which establishes the full optimality half (closed sorry-free, see
Issue #4).
-/
theorem limit_point_mem_constraintSet (cfg : SystemConfig n)
    (traj : Trajectory n cfg)
    (W_star : Matrix (Fin n) (Fin n) ℝ)
    (h_limit : Filter.Tendsto traj.weights Filter.atTop (nhds W_star)) :
    W_star ∈ constraintSet n cfg.support cfg.budget := by
  exact IsClosed.mem_of_tendsto ( constraintSet_isClosed _ _ ) h_limit ( Filter.Eventually.of_forall fun x => traj.hW_in_C x )

/-- The Lyapunov function `L(θ, ·)` is convex in the weight matrix `W`:
  for any `V`, the first-order Taylor expansion at `W` gives a lower bound.

  Equivalently, `L(θ, V) ≥ L(θ, W) + ⟨∇_W L(θ, W), V − W⟩_{upperTri}` for the
  upper-triangle inner product. The convexity residual equals
  `(λ/2) · Σ_{i<j} (V_ij − W_ij)² ≥ 0` since `λ > 0` (from `cfg.hlam_pos`).

  This is the standard convexity certificate used in the KKT-stationarity
  corollary: at a fixed point `Ẇ = 0`, the variational inequality from
  `hWeight_dyn` plus convexity yields global minimisation. -/
private lemma lyapunovFn_convex_in_W (cfg : SystemConfig n) (theta : Fin n → ℝ)
    (W V : Matrix (Fin n) (Fin n) ℝ) :
    lyapunovFn n (-1) cfg.decayRate theta W +
      (∑ p ∈ upperTriPairs n,
        lyapunovGradW (-1) cfg.decayRate theta W p.1 p.2 *
        (V p.1 p.2 - W p.1 p.2)) ≤
    lyapunovFn n (-1) cfg.decayRate theta V := by
  unfold lyapunovFn
  simp only [Finset.mul_sum]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_le_sum fun p _ => ?_
  simp only [lyapunovGradW, Matrix.of_apply]
  have hlam : 0 < cfg.decayRate := cfg.hlam_pos
  have h_sq : 0 ≤ (V p.1 p.2 - W p.1 p.2)^2 := sq_nonneg _
  nlinarith [h_sq, hlam]

/--
**(Issue #4 — closed sorry-free, via Fermat-stationary fixed-point argument.)**
Full KKT stationarity at a *weight fixed point* of the trajectory.

If at some time `t_star` the weight derivative vanishes (`Wdot t_star = 0`),
then the state `(phases t_star, weights t_star)` is KKT-stationary: the weight
matrix minimises `L(phases t_star, ·)` over the constraint set `C`.

This is the *limit-point* interpretation in the sense that a converged
Lyapunov-descending trajectory has `Ẇ → 0`; once `Ẇ = 0` is attained at any
specific time, the surrounding state is globally optimal for the weight
sub-problem.

**Argument.**
1. *Primal feasibility:* `weights t_star ∈ C` directly from `hW_in_C`.
2. *Variational inequality at the fixed point.* From `hWeight_dyn` at time
   `t_star` with `V` arbitrary in `C`, substituting `Wdot t_star = 0` gives
     `0 ≤ ∑_{i,j} (V_ij − W_ij) · (η · ∇_W L_ij)`,
   so `0 ≤ ⟨V − W, η · ∇_W L⟩_F`. Since `η > 0`, `⟨V − W, ∇_W L⟩_F ≥ 0`.
3. *Upper-triangle restriction.* `V − W` is symmetric and zero-diagonal
   (linear inheritance from `V, W ∈ C`); `∇_W L` is symmetric
   (`lyapunovGradW_isSymm`). By `sum_full_eq_twice_upperTri`,
   `⟨V − W, ∇_W L⟩_F = 2 · Σ_{i<j} (V_ij − W_ij) · ∇_W L_ij`, so
   `Σ_{i<j} (V_ij − W_ij) · ∇_W L_ij ≥ 0`.
4. *Convexity of L in W* (`lyapunovFn_convex_in_W`) gives
   `L(V) ≥ L(W) + ⟨∇_W L(W), V − W⟩_{upperTri} ≥ L(W)`. -/
theorem limit_point_isKKTStationary (cfg : SystemConfig n)
    (traj : Trajectory n cfg) (t_star : ℝ)
    (h_Wdot_zero : traj.Wdot t_star = 0) :
    IsKKTStationary n cfg (traj.phases t_star) (traj.weights t_star) := by
  refine ⟨traj.hW_in_C t_star, fun V hV => ?_⟩
  set W := traj.weights t_star with hW_def
  set θ := traj.phases t_star with hθ_def
  set gL := lyapunovGradW (-1 : ℝ) cfg.decayRate θ W with hgL_def
  -- Step 1: variational inequality at the fixed point.
  -- hWeight_dyn at t_star with V, simplified by Wdot t_star = 0.
  have h_zero : ∀ i j, traj.Wdot t_star i j = 0 := fun i j => by
    rw [h_Wdot_zero]; rfl
  have h_VI_raw : 0 ≤ ∑ i, ∑ j, (V i j - W i j) *
      (cfg.learnRate * gL i j) := by
    have h := traj.hWeight_dyn t_star V hV
    simp only [h_zero, zero_add] at h
    exact h
  -- Factor out η to get ⟨V − W, ∇L⟩_F ≥ 0.
  have h_eta : 0 < cfg.learnRate := cfg.heta_pos
  have h_VI_full : 0 ≤ ∑ i, ∑ j, (V i j - W i j) * gL i j := by
    have h_factor : (∑ i, ∑ j, (V i j - W i j) * (cfg.learnRate * gL i j)) =
                    cfg.learnRate * ∑ i, ∑ j, (V i j - W i j) * gL i j := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j _ => ?_
      ring
    rw [h_factor] at h_VI_raw
    nlinarith [h_VI_raw]
  -- Step 2: upper-triangle restriction.
  have hVsymm : V.IsSymm := hV.1
  have hWsymm : W.IsSymm := (traj.hW_in_C t_star).1
  have hVdiag : ∀ i, V i i = 0 := hV.2.2.2.1
  have hWdiag : ∀ i, W i i = 0 := (traj.hW_in_C t_star).2.2.2.1
  -- (V - W) is symmetric and zero-diagonal.
  have hVWsymm : (V - W).IsSymm := by
    ext i j
    simp only [Matrix.transpose_apply, Matrix.sub_apply]
    rw [hVsymm.apply i j, hWsymm.apply i j]
  have hVWdiag : ∀ i, (V - W) i i = 0 := fun i => by
    simp only [Matrix.sub_apply, hVdiag i, hWdiag i, sub_zero]
  have hgL_symm : gL.IsSymm := lyapunovGradW_isSymm cfg.decayRate θ hWsymm
  -- Bridge entry-wise to matrix-sub form.
  have h_full_eq_VW : (∑ i, ∑ j, (V i j - W i j) * gL i j) =
                       ∑ i, ∑ j, (V - W) i j * gL i j := by
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [Matrix.sub_apply]
  rw [h_full_eq_VW] at h_VI_full
  have h_full_eq_twice :
      ∑ i, ∑ j, (V - W) i j * gL i j =
      2 * ∑ p ∈ upperTriPairs n, (V - W) p.1 p.2 * gL p.1 p.2 :=
    sum_full_eq_twice_upperTri hVWsymm hgL_symm hVWdiag
  rw [h_full_eq_twice] at h_VI_full
  -- Upper-tri sum is non-negative (after halving).
  have h_uptri_VW_nn : 0 ≤ ∑ p ∈ upperTriPairs n, (V - W) p.1 p.2 * gL p.1 p.2 := by
    linarith
  -- Convert back from (V - W) p.1 p.2 to V p.1 p.2 - W p.1 p.2.
  have h_uptri_entry_eq :
      (∑ p ∈ upperTriPairs n, (V - W) p.1 p.2 * gL p.1 p.2) =
      ∑ p ∈ upperTriPairs n, (V p.1 p.2 - W p.1 p.2) * gL p.1 p.2 := by
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Matrix.sub_apply]
  rw [h_uptri_entry_eq] at h_uptri_VW_nn
  -- Step 3: convexity of L in W, plus the upper-triangle non-negativity.
  have h_convex := lyapunovFn_convex_in_W cfg θ W V
  have h_swap :
      (∑ p ∈ upperTriPairs n, gL p.1 p.2 * (V p.1 p.2 - W p.1 p.2)) =
      ∑ p ∈ upperTriPairs n, (V p.1 p.2 - W p.1 p.2) * gL p.1 p.2 :=
    Finset.sum_congr rfl fun _ _ => mul_comm _ _
  rw [h_swap] at h_convex
  linarith

end KuramotoHebbian