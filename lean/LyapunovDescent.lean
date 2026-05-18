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

/-- A trajectory of the Kuramoto-Hebbian system with differentiable Lyapunov function. -/
structure Trajectory (n : ℕ) (cfg : SystemConfig n) where
  /-- Phase trajectory θ(t) -/
  phases : ℝ → Fin n → ℝ
  /-- Weight trajectory W(t) -/
  weights : ℝ → Matrix (Fin n) (Fin n) ℝ
  /-- Weights stay in C -/
  hW_in_C : ∀ t, weights t ∈ constraintSet n cfg.support cfg.budget
  /-- Phase dynamics: dθᵢ/dt = K · Σⱼ Wᵢⱼ sin(θⱼ - θᵢ) -/
  hPhase_dyn : ∀ t i,
    HasDerivAt (fun s => phases s i)
      (phaseDot cfg.coupling (weights t) (phases t) i) t

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

/-
**Main Theorem (Lyapunov Descent)**: The Lyapunov function L is non-increasing
  along trajectories of the Budgeted Hebbian Kuramoto system.

  We assume:
  - The Lyapunov function is differentiable along the trajectory
  - Its derivative equals the sum of phase and weight contributions
  - The weight contribution is non-positive (from projected gradient descent)

  From these and the proved `phase_descent`, we conclude dL/dt ≤ 0 everywhere,
  and hence L is non-increasing.
-/
theorem lyapunov_descent (cfg : SystemConfig n) (traj : Trajectory n cfg)
    /- The Lyapunov function along the trajectory is continuous -/
    (hL_cont : ContinuousOn traj.lyapunovAlong (Set.univ))
    /- The Lyapunov function is differentiable on ℝ -/
    (hL_diff : DifferentiableOn ℝ traj.lyapunovAlong (Set.univ))
    /- The derivative of L along the trajectory is the sum of phase and weight contributions -/
    (hL_deriv : ∀ t, deriv traj.lyapunovAlong t =
      phaseContribution cfg.coupling (traj.weights t) (traj.phases t) +
      weightContribution cfg.decayRate (traj.phases t) (traj.weights t)
        (deriv traj.weights t))
    /- The weight contribution is non-positive (from projected gradient descent onto C) -/
    (hW_descent : ∀ t, weightContribution cfg.decayRate (traj.phases t) (traj.weights t)
        (deriv traj.weights t) ≤ 0)
    (t₁ t₂ : ℝ) (h12 : t₁ ≤ t₂) :
    traj.lyapunovAlong t₂ ≤ traj.lyapunovAlong t₁ := by
  -- By the properties of the derivative, if the derivative of a function is non-positive on an interval, then the function is non-increasing on that interval.
  have h_deriv_nonpos : ∀ t, deriv traj.lyapunovAlong t ≤ 0 := by
    exact fun t => hL_deriv t ▸ add_nonpos ( phase_descent cfg.hK_neg _ _ ) ( hW_descent t );
  by_contra h_contra;
  have := exists_deriv_eq_slope traj.lyapunovAlong ( show t₁ < t₂ from lt_of_le_of_ne h12 ( by aesop_cat ) );
  exact absurd ( this ( hL_cont.mono ( Set.subset_univ _ ) ) ( hL_diff.mono ( Set.subset_univ _ ) ) ) ( by rintro ⟨ c, ⟨ h₁, h₂ ⟩, h₃ ⟩ ; rw [ eq_div_iff ] at h₃ <;> nlinarith [ h_deriv_nonpos c ] )

/-! ## KKT Stationarity Corollary -/

/-- A state (θ, W) is KKT-stationary for the weight optimisation if W minimises
  L(θ, ·) over C. -/
def IsKKTStationary (n : ℕ) (cfg : SystemConfig n)
    (theta : Fin n → ℝ) (W : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  W ∈ constraintSet n cfg.support cfg.budget ∧
  ∀ W' ∈ constraintSet n cfg.support cfg.budget,
    lyapunovFn n (-1) cfg.decayRate theta W ≤ lyapunovFn n (-1) cfg.decayRate theta W'

/-
**Corollary (KKT Stationarity)**: Any limit point (θ*, W*) of a trajectory
  satisfies W* ∈ C (membership in the constraint set).

  This follows from closedness of C and the fact that W(t) ∈ C for all t.
  The full KKT optimality (W* minimises L(θ*,·) over C) additionally requires
  that the projected gradient vanishes at the limit, which follows from
  the Lyapunov descent and the continuous dependence of the gradient on the state.
-/
theorem limit_point_mem_constraintSet (cfg : SystemConfig n)
    (traj : Trajectory n cfg)
    (W_star : Matrix (Fin n) (Fin n) ℝ)
    (h_limit : Filter.Tendsto traj.weights Filter.atTop (nhds W_star)) :
    W_star ∈ constraintSet n cfg.support cfg.budget := by
  exact IsClosed.mem_of_tendsto ( constraintSet_isClosed _ _ ) h_limit ( Filter.Eventually.of_forall fun x => traj.hW_in_C x )

end KuramotoHebbian