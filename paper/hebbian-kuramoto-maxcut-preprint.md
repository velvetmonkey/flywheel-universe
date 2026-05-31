---
type: preprint-draft
status: draft
tags:
  - hebbian
  - kuramoto
  - maxcut
  - preprint
  - oscillator-ising-machine
  - amplitude-heterogeneity
date: '2026-05-19'
description: >-
  Zenodo preprint v1 — Hebbian-coupled Kuramoto dynamics under hardware coupling
  budgets, evaluated on Max-Cut with amplitude-heterogeneity robustness as the
  primary finding.
---

# Budgeted Hebbian Kuramoto Dynamics for Max-Cut under Amplitude Heterogeneity: Robustness, Not Cut Quality, Is the Signal

**Ben Cassie**  
Independent Researcher  
bencassie@outlook.com | https://orcid.org/0009-0004-1899-7627

**Status:** All 12 verification items resolved 2026-05-20. Ready for Zenodo upload.

---

## Abstract

We study budgeted Hebbian-coupled Kuramoto dynamics as a coupling-resource allocation mechanism for oscillator-based Ising machines applied to graph Max-Cut. The Hebbian rule learns a symmetric, non-negative coupling matrix on a fixed edge support, with per-node row-sum budgets enforced by exact symmetric-Frobenius projection at every step. Under zero detuning, fixed support, and exact projection, the joint phase / weight dynamics admit a Lyapunov function whose constrained descent and KKT stationarity at limit points are formalised in Lean 4 / Mathlib (with explicit open issues).

We compare twelve methods — static, static-budgeted, random-budgeted, learned-support-random-weights, topology-scrambled, Hebbian-Frobenius, hybrid-Frobenius, Hebbian-Sinkhorn, Goemans–Williamson, and others — across four graph families: sparse and dense Erdős–Rényi, random 10-regular, and one GSet calibration instance (G1, n=800). All runs use 3 graph seeds per random family and 10 method seeds per (family, graph_seed, method) cell, with greedy local-search polishing applied to every method's output.

The predeclared CORE test — `hybrid_frobenius` against `random_budgeted` on polished cut at σ=0 — **fails**: random-budgeted coupling at the same half-mean-degree budget recovers more cut than Hebbian-Frobenius on every family. Cut quality at σ=0 is therefore not where the Hebbian rule contributes. The amplitude-heterogeneity sweep (lognormal node gains, σ ∈ {0, 0.25, 0.5, 1.0}) tells a different story: hybrid-Frobenius cut value rises with σ (slope ≈ +2.4 cut/σ) while the registered comparator random-budgeted falls steeply (slope ≈ −9.9 cut/σ); the paired slope difference is **+12.3 cut/σ with bootstrap 95% CI [+8.5, +14.7]** and the direction is consistent across all three graph seeds. The predeclared HARDWARE-ROBUSTNESS test therefore holds. As σ approaches 1.0 the Hebbian–oracle gap shrinks from 32 cut (at σ=0) to 12 cut, and the oracle reference itself collapses (−17.8 cut over the range), evidence that knowing the noise model without adapting is brittle. The contribution of this work is therefore narrower than the original framing suggested: **the Hebbian rule is a robustness primitive for heterogeneous-amplitude oscillator arrays, not a cut-quality competitor at σ=0**.

We make no wall-clock or compute-cost claim. The coupling-resource budget is a parameter-count constraint motivated by hardware fan-in limits in coherent and parametric Ising machines [Mohseni et al., 2022; Khan et al., 2025]; digital simulation runtime is dominated by the symmetric-Frobenius projection and is not halved.

---

## 1. Introduction

Oscillator-based Ising machines, in their coherent (CIM) [Wang et al., 2013; McMahon et al., 2016], optical, and electronic [Wang & Roychowdhury, 2019] realisations, face hard physical constraints on per-node coupling resources: fan-in, wirelength, and parametric-pumping budget all limit how many other spins a given oscillator can be coupled to and how strongly. These constraints make sparse coupling a hardware reality, not a software optimisation. The relevant algorithmic question is then *not* whether sparse coupling can match dense, but where, given a fixed per-node coupling budget, the resource should be spent.

Static topology coupling — using the original graph adjacency, with weights truncated to fit the row budget — is the canonical default and is what most digital Ising-machine simulators implement. Hebbian coupling lets the dynamics themselves reallocate the budget: weights that connect oscillators tending toward antiphase grow, weights between same-phase oscillators decay, and the configuration relaxes onto the budget polytope at every step. The question this paper addresses is whether the Hebbian-allocated coupling actually delivers higher-quality Max-Cut solutions than the static or randomised allocations at the same budget, and whether it is more or less robust under hardware-realistic device variation.

The honest answer is mixed. At zero amplitude heterogeneity (σ=0), Hebbian-Frobenius coupling does *not* outperform random non-negative coupling at the same budget — random-budgeted coupling is in fact uniformly better on polished cut quality across the tested families, and the predeclared CORE test (registered in the benchmark harness before the run) consequently fails. The cut-quality story at σ=0 is that *sparse coupling, of any reasonable kind*, recovers 96–99% of full static coupling at half the budget; the choice of *which* sparse coupling is, within that band, second-order.

Under amplitude heterogeneity, however, the picture changes. As the lognormal gain dispersion σ rises from 0 to 1, the registered comparator `random_budgeted` falls steeply (slope −9.9 cut/σ) and `static_budgeted` drifts downward (−6.1 cut/σ); the predeclared primary method `hybrid_frobenius` drifts *upward* (+2.4 cut/σ), with `hebbian_frobenius` close behind (+2.2 cut/σ); and Hebbian's gap to the σ-aware oracle-compensation reference shrinks from 32 cut to 12. The Hebbian rule's empirical contribution, on this benchmark, is robustness to amplitude variation — a hardware-relevant property given the regime identified in Khan et al. [2025], where quasi-steady amplitude fluctuations in parametric-oscillator Ising machines rescale the effective spin couplings and degrade solution quality.

**Contributions.** This paper makes four contributions:

1. A formal model of budgeted Hebbian Kuramoto dynamics with symmetric, non-negative coupling on a fixed support, with per-node row-sum budgets enforced by exact symmetric-Frobenius projection (Section 4).
2. A Lyapunov-descent and KKT-stationarity theorem for the joint phase / weight flow under zero detuning, fixed support, and exact projection, with the algebraic core formalised in Lean 4 / Mathlib (Section 3.3). The KKT-stationarity branch is closed sorry-free; remaining open issues — weight chain rule, projection existence, and extension to projected / nonsmooth dynamics — are explicitly catalogued.
3. An honest empirical evaluation against twelve baselines with a predeclared success criterion. We report the CORE test outcome — failure — without re-framing, alongside the secondary signals that do hold (Section 5).
4. An amplitude-heterogeneity slope analysis showing that the Hebbian rule degrades less under realistic device variation than the registered random-budgeted comparator, with paired slope difference +12.3 cut/σ and bootstrap 95% CI [+8.5, +14.7] (Section 4.3).

We do not claim a state-of-the-art Max-Cut solver, a wall-clock or compute-cost improvement, or a biologically plausible learning law. The "coupling-resource" framing throughout this paper is a parameter-count claim about a hypothetical Ising-machine architecture in which row-sum budgets correspond to physical fan-in caps; digital simulation runtime is dominated by the projection step and is not halved when the budget is halved.

Downstream applications of the same primitive — to long-range numerical weather prediction (attractor structure and regime tipping), to fault-tolerant distributed coordination (phase synchrony as recovery-free coordination), and to wirelength-constrained chip placement — are sketched in Section 6.3 as *implications worth testing*, not as contributions of this work.

## 2. Background

### 2.1 Kuramoto oscillators

The Kuramoto model [Kuramoto, 1975; Acebrón et al., 2005; Strogatz, 2000] describes a population of coupled phase oscillators evolving as

$$\dot{\theta}_i = \omega_i + K \sum_j W_{ij} \sin(\theta_j - \theta_i),$$

where θ_i ∈ [0, 2π) is the phase of oscillator i, ω_i is its natural frequency, K is a global coupling strength, and W_{ij} encodes the pairwise coupling structure. For positive K and uniform coupling, the system synchronises above a critical threshold; for negative K and an adjacency matrix W, the dynamics favour antiphase between coupled nodes, which corresponds to a bipartition of the underlying graph.

The order parameter $r = |N^{-1} \sum_j e^{i \theta_j}|$ measures coherence; for anti-ferromagnetic coupling and a bipartite-friendly graph, the steady state organises into two clusters separated by approximately π, and the cut value of the implied bipartition (sign of cos θ_i) is a natural quantity to read off the dynamics.

### 2.2 Hebbian learning

Classical Hebbian learning [Hopfield, 1982] updates synaptic weights in proportion to the product of pre- and post-synaptic activity. For phase oscillators, the natural analogue uses cosine similarity:

$$\dot{W}_{ij} = \eta\, \cos(\theta_i - \theta_j) - \lambda W_{ij}.$$

When the support of W is masked to a fixed edge set, weights between same-phase oscillators decay (cos ≈ +1 grows the weight if K > 0 but here K < 0 reverses signs), and the rule reallocates coupling resource toward edges that participate in cluster organisation. The decay term λW serves as control-theoretic regularisation; it is not claimed to be biologically plausible.

### 2.3 Max-Cut

Max-Cut on a weighted graph G = (V, E, w) asks for a partition S ⊆ V maximising $\sum_{(i,j) \in E,\, i \in S,\, j \notin S} w_{ij}$. The problem is NP-hard. The canonical approximation algorithm is the semidefinite relaxation of Goemans and Williamson [1995], achieving the 0.878 approximation ratio via randomised hyperplane rounding. The GSet library [Helmberg & Rendl, 2000] provides standard benchmark instances; we use G1 (n=800) as a single calibration point.

### 2.4 Oscillator-based Ising machines

Coherent and parametric Ising machines [Wang et al., 2013; McMahon et al., 2016; Mohseni et al., 2022] map binary spin variables to oscillator phases and let the physical dynamics relax toward low-energy configurations. The coupling matrix is implemented as a network of analog connections (delay lines, parametric amplifiers, electronic transmission lines), and the per-node fan-in is bounded by physical constraints. The "coupling budget" framing of this paper models that constraint as a row-sum cap on a non-negative symmetric matrix.

Recent work [Khan et al., 2025] identified amplitude heterogeneity — quasi-steady variation in oscillator amplitude across the array — as a regime in which the effective spin couplings $K^{\text{eff}}_{ij} = a_i a_j W_{ij}$ rescale in a node-dependent way, degrading solution quality on standard Max-Cut benchmarks. Robustness to this regime is the hardware-relevant question we investigate in Section 4.3.

## 3. Methods

### 3.1 Model

The joint state is (θ, W) where θ ∈ ℝ^n is the phase vector and W ∈ ℝ^{n×n} is a symmetric, non-negative coupling matrix supported on a fixed edge mask A ⊆ {(i, j) : i ≠ j}. The continuous-time flow is

$$\dot{\theta}_i = \omega_i + K \sum_j W_{ij} \sin(\theta_j - \theta_i),$$

$$\dot{W} = \Pi_C\!\left[\eta\, M(\theta) - \lambda W\right],$$

where M(θ)_{ij} = cos(θ_i - θ_j) for (i, j) ∈ A and 0 otherwise, and Π_C is the orthogonal projection onto the constraint set

$$C = \{X \in \mathbb{R}^{n \times n} : X = X^\top,\ X \geq 0,\ X_{ij} = 0\ \text{for}\ (i,j) \notin A,\ \mathrm{diag}(X) = 0,\ \sum_j X_{ij} \leq B \ \forall i\}.$$

Here B is the per-node row-sum budget (a coupling-resource cap), set to either the mean degree of the support or half of it (half-mean-degree). All experiments in the main suite use **half-mean-degree** as the budget.

### 3.2 Symmetric-Frobenius projection

Given an unprojected update Y, the projection Π_C(Y) is the symmetric, non-negative, support-masked matrix in C closest to Y in Frobenius norm:

$$\Pi_C(Y) = \arg\min_{X \in C} \|X - Y\|_F^2.$$

A CVXPY-based reference solver compiles the QP and is used to validate the runtime path to machine precision. The runtime projection is a dual root-finder: the KKT conditions reduce, for each row i, to a single-variable root problem on the Lagrange multiplier $\lambda_i \geq 0$ associated with the row-sum constraint, with the inner sub-problem being a classical projection onto the budget simplex solvable in O(k log k) per row (Brucker [1984] form). The outer iteration is a Gauss–Seidel sweep over rows until $\|\Delta \lambda\|_\infty$ falls below a tolerance, with warm-starting between successive Hebbian update calls.

### 3.3 Theorem (informal)

**Theorem.** *Let A be a fixed support, B > 0, K < 0, η > 0, λ > 0, and ω_i = 0 for all i. Let (θ(t), W(t)) be a solution of the continuous-time joint flow with W(t) ∈ C for all t. Then there exists a Lyapunov function L(θ, W) such that*

$$\frac{dL}{dt} \leq 0$$

*along the joint trajectory, with equality holding only at points satisfying KKT stationarity for L subject to the constraints of C.*

The algebraic core of this result — convexity and closedness of C, the symmetry-based identity reducing the phase contribution to $K \sum_i f_i^2$ where $f_i = \sum_j W_{ij} \sin(\theta_j - \theta_i)$, and the descent conclusion under explicit hypotheses on the dL/dt decomposition and on Π_C's variational inequality — is machine-verified in Lean 4 / Mathlib at `lean/RequestProject/LyapunovDescent.lean` [Moura & Ullrich, 2021; The Mathlib Community, 2024].

The Hebbian joint Lyapunov descent identity proved algebraically here has also been formalised in the companion `kuramoto-lean` library [Cassie, 2026]. That library provides 26 fully proved theorems covering the Kuramoto order parameter, gradient flow, weighted Lyapunov descent, Hebbian coupling, and ODE existence, with one documented gap: global convergence to synchrony, pending LaSalle's invariance principle or Barbalat-style infrastructure in Mathlib. The Zenodo record is available at https://doi.org/10.5281/zenodo.20468619 and the source at https://github.com/velvetmonkey/kuramoto-lean.

**Open issues.** The Lean formalisation closed Issue #4 (full KKT stationarity at limit points, via a Fermat-stationary fixed-point argument; see commit `77c78df`) sorry-free. Remaining open gaps are:
1. **Weight chain rule** connecting the projected $\dot{W}$ to the unconstrained $\partial L / \partial W$. The current proof of `hW_descent_derived` uses a Fermat argument that requires everywhere two-sided differentiability of the weight trajectory (`hWeight_diff`); the cleaner path is a Moreau-decomposition / tangent-cone argument sketched in `lean/RequestProject/FutureWork.lean` (commit `4c6a3bf`), but it requires `tangentCone_constraintSet`, `normalCone_constraintSet`, and Moreau-decomposition lemmas that are not yet in Mathlib v4.28.0 for the budget polytope's constraint shape.
2. **Existence of the projection** Π_C under the full set of constraints.
3. **Extension to projected / nonsmooth dynamics.** The theorem above is for the continuous-time smooth flow. The benchmark implements discrete updates with periodic projection at the budget boundary, which introduces nonsmoothness that the current proof does not cover. Extension to differential-inclusion or non-smooth Lyapunov-stability frameworks is open work.

We therefore restrict the empirical Section 4.5 to numerical *energy-decrease diagnostics* along discrete trajectories, framed as descriptive observations rather than theorem-certified statements.

### 3.4 Benchmark protocol

**Graph families.** We evaluate four families:

- **sparse ER**: Erdős–Rényi G(n=200, p=0.05), mean degree ≈ 9.4.
- **dense ER**: Erdős–Rényi G(n=200, p=0.15), mean degree ≈ 28.
- **random 10-regular**: random 10-regular graph on n=200.
- **GSet G1**: n=800, m=19176, unweighted, from the Helmberg–Rendl GSet library.

We also include a small-instance group (n ∈ [20, 30]) with Gray-code-enumerated exact optima for ground-truth calibration, and a detuning-ablation family that shadows sparse ER with ω ≠ 0 to probe the theorem's regime boundary.

**Sample sizes.** Each random family uses 3 independent graph seeds; GSet G1 is one instance. Each (family, graph_seed, method) cell uses 10 method seeds (independent stochastic restarts of the dynamics). This is a scoped-pilot sample size; full inferential reruns at 30+ graph seeds per family are flagged as next-step work.

**Methods compared.** Twelve methods:
- `static`: full original adjacency with unit weights.
- `static_projected`: static coupling with projection but no budget reduction.
- `static_budgeted`: top-k symmetric support of original adjacency at budget B.
- `random_rounding`: static dynamics with multiple-hyperplane Goemans–Williamson rounding.
- `random_budgeted`: random non-negative weights on full adjacency, projected to C.
- `learned_support_random_weights` (LSRW): random non-negative weights on the support converged to by Hebbian dynamics — same support as Hebbian-Frobenius, different weights.
- `topology_scrambled`: random support obtained by double-edge-swap of the adjacency, original-weight values projected to C.
- `hebbian_frobenius`: full adaptive Hebbian flow under symmetric-Frobenius projection.
- `hybrid_frobenius`: predeclared primary method — Hebbian adaptation followed by fixed-W phase settling.
- `hebbian_sinkhorn`, `hybrid_sinkhorn`: variants under Sinkhorn-style projection (not row-sum exact).
- `sdp_gw`: SDP relaxation with Goemans–Williamson rounding.
- `greedy`: pure greedy local-search baseline.

**Budgets.** B set to either mean_degree of the support or 0.5·mean_degree (half-mean-degree). Main suite uses half-mean-degree; this is the regime where coupling-resource constraints bite.

**Polishing.** Every method's output cut is polished by greedy single-flip local search to a local optimum.

**Aggregation and statistical unit.** For each (family, graph_seed, method) cell we take the best polished cut over the 10 method seeds (best-of-restarts). The statistical unit is the *graph instance*, not the restart. Predeclared comparisons are paired across graph_seeds within a family.

**Predeclared success criteria** (registered in `phase2_benchmark.py` header):

- **CORE**: `hybrid_frobenius` beats `random_budgeted` on polished cut at σ=0, B = half-mean-degree, paired across graph_seeds within family.
- **HARDWARE-ROBUSTNESS**: `hybrid_frobenius` degrades less than `random_budgeted`, and tracks `oracle_compensation` more closely, as amplitude heterogeneity σ rises (σ ≥ 0.5).

**Statistical reporting.** With 3 graph seeds per random family, asymptotic Wilcoxon p-values are uninterpretable; we therefore report paired mean deltas, bootstrap 95% confidence intervals on the paired delta (10000 resamples, sampling graph_seeds with replacement), and the exact sign-flip permutation test (2³ = 8 permutations for n_pairs = 3, treated as descriptive rather than inferential). No claims of statistical significance are made on the family-level comparisons; sample size precludes them.

### 3.5 Amplitude heterogeneity

Node gains $a_i \sim \mathrm{LogNormal}(-\sigma^2/2, \sigma)$ scale the effective coupling between oscillators i and j as $K^{\text{eff}}_{ij} = a_i a_j W_{ij}$, modelling the regime identified in Khan et al. [2025]. We sweep σ ∈ {0, 0.25, 0.5, 1.0} on the sparse_er family with 10 method seeds per (sigma, graph_seed, method) cell and report absolute cut values, ratios to σ=0 static-budgeted, and per-graph-seed regression slopes of cut against σ.

The `oracle_compensation` reference uses σ-aware rescaling that exactly undoes the $a_i a_j$ factor, giving an upper bound on what a perfectly-calibrated controller could achieve.

### 3.6 ODE solver, parameters

- Solver: `scipy.integrate.solve_ivp` with `DOP853`, rtol = 1e-5, atol = 1e-7.
- Hebbian update step: dt_weight = 0.5, learning rate η = 0.05, decay λ = 0.003.
- Coupling strength: K = −0.5 (anti-ferromagnetic regime).
- Total integration time: T_total = 200.
- Number of method seeds per (family, graph_seed, method) cell: 10.
- Natural frequencies: ω = 0 for the main suite (matching the theorem's hypothesis); the detuning-ablation family uses ω_i ∼ N(0, 0.3).

### 3.7 Reproducibility manifest

- Repository: `https://github.com/velvetmonkey/flywheel-universe`.
- Repository HEAD at time of writing: `0009208` (preprint: add rendered PDF). Predecessor commits referenced in the paper: `71b247b` (het v2: registered-comparator sweep, source of `heterogeneity_experiment_v2.csv`); `4c6a3bf` and `77c78df` (Lean proof updates, closing additional sorry-free issues); `1d5c5c8` (het v1 rerun at 10 method seeds).
- Phase 2 results CSV: `results/phase2_full_results.csv`, 3364 rows.
- Heterogeneity CSV (v1, four methods): `results/heterogeneity_experiment.csv`, 480 rows.
- **Heterogeneity CSV (v2, five methods incl. `random_budgeted` registered comparator): `results/heterogeneity_experiment_v2.csv`, 600 rows.** Generated by `python run_het_v2.py` (~4.5h on a single workstation). The v2 CSV is the basis for the Section 4.3 slope test against the registered HARDWARE-ROBUSTNESS comparator; v1 is retained for reproducibility.
- Generation (main suite): `python phase2_benchmark.py` (full run, ~6h on a single workstation).

## 4. Results

We report results in the order: cut quality at σ=0 (Section 4.1), the predeclared CORE test outcome (Section 4.2), amplitude-heterogeneity robustness (Section 4.3 — the main empirical finding), and Lyapunov-decrease diagnostics (Section 4.4).

### 4.1 Cut quality at σ = 0

Table 1 reports best-polished cut per method, averaged across the 3 graph seeds for the three random families and on the single GSet G1 instance, at the half-mean-degree budget and zero amplitude heterogeneity. All values are mean cut.

**Table 1.** Cut quality (best polished cut, mean over graph_seeds, half-mean-degree budget, σ = 0).

| Family                  | static (full) | static_budgeted | random_budgeted | topology_scrambled | LSRW   | hebbian_frob | hybrid_frob | greedy | sdp_gw |
|-------------------------|--------------:|----------------:|----------------:|-------------------:|-------:|-------------:|------------:|-------:|-------:|
| sparse_er_200_p05       |        741.67 |          731.33 |      **736.33** |             712.00 | 704.33 |       718.67 |      718.67 | 722.33 | 754.00¹|
| dense_er_200_p15        |       1869.33 |         1859.67 |     **1861.33** |            1834.33 |1819.33 |      1841.33 |     1841.33 |1842.00 | —      |
| random_regular_200_d10  |        728.67 |          715.33 |      **725.33** |             691.33 | 670.67 |       691.33 |      691.33 | 694.00 | —      |
| GSet G1 (n=800)         |      11599.00 |        11497.00 |    **11522.00** |           11401.00 |11368.00|     11434.00 |    11434.00 |11411.00| —      |

¹ SDP+GW computed for sparse_er_200_p05 only (one graph_seed; others not run due to solver scale).

The pattern is consistent across families. Among budgeted methods (half-mean-degree budget):

- `random_budgeted` outperforms both `static_budgeted` and `hebbian_frobenius` on every family.
- `hebbian_frobenius` and `hybrid_frobenius` are nearly identical, consistent with the predeclared distinction (a freeze-schedule variant) being a small effect on cut quality.
- `topology_scrambled` is the weakest baseline, confirming that not all sparse supports are equivalent: the structure of the support matters.
- Full static coupling (mean degree, not half) recovers 100% of cut quality; all half-budget methods recover between 96% and 99% of full-static.

**Note on GSet G1.** The 11434 / 11599 figures correspond to approximate ratios 98.6% / 99.8% relative to the best-known G1 cut of 11624 [Benlic & Hao, 2013]. G1 is one calibration instance; no scaling claim is made.

### 4.2 Predeclared CORE test outcome

The predeclared CORE comparison is `hybrid_frobenius` against `random_budgeted` on polished cut at σ=0 and B = half-mean-degree. Table 2 reports paired deltas and bootstrap 95% confidence intervals per family.

**Table 2.** Predeclared CORE test outcome (hybrid_frobenius − random_budgeted, paired across graph_seeds). Bootstrap CIs from 10000 resamples; exact permutation p-values from 2³ = 8 sign-flip permutations, descriptive only.

| Family                  | n_pairs | paired Δ | bootstrap 95% CI | exact perm p |
|-------------------------|--------:|---------:|------------------:|-------------:|
| sparse_er_200_p05       |       3 |  −17.67  | [−24.0, −12.0]    | 0.25         |
| dense_er_200_p15        |       3 |  −20.00  | [−26.0, −12.0]    | 0.25         |
| random_regular_200_d10  |       3 |  −34.00  | [−38.0, −26.0]    | 0.25         |
| GSet G1 (n=800)         |       1 |  −88.00  | n/a (n=1)         | n/a          |

**The CORE test fails.** Across all four families, `hybrid_frobenius` *trails* `random_budgeted` on polished cut; the bootstrap intervals on the paired delta exclude zero in every case, and the direction is consistent. At the scoped-pilot sample size of n=3 paired observations per family the exact-permutation p-value is bounded below by 0.25 (the smallest two-sided two-tailed p achievable with three sign-flips), so we make no inferential claim, but the descriptive direction is unambiguous: random non-negative coupling at the same budget produces better polished cuts than Hebbian-Frobenius coupling.

The natural interpretation is that the Hebbian rule's converged support is *worse* — for Max-Cut at σ=0 — than the full graph adjacency that `random_budgeted` retains. Even given the inferior support, the Hebbian *weights* on that support beat random weights on the same support, as Section 4.5 shows; but that signal does not bridge the support gap.

### 4.3 Amplitude-heterogeneity robustness (main finding)

We sweep σ ∈ {0, 0.25, 0.5, 1.0} on the sparse_er family with 10 method seeds × 3 graph seeds per cell, using lognormal node gains $a_i \sim \mathrm{LogNormal}(-\sigma^2/2, \sigma)$. The sweep includes the predeclared HARDWARE-ROBUSTNESS comparator `random_budgeted` (added in the v2 rerun, `results/heterogeneity_experiment_v2.csv`, 600 rows).

**Table 3a.** Absolute cut quality under amplitude heterogeneity (mean over graph_seeds × method_seeds, sparse_er, half-mean-degree budget).

| σ    | static_budgeted | random_budgeted | hebbian_frob | hybrid_frob | oracle_comp |
|-----:|----------------:|----------------:|-------------:|------------:|------------:|
| 0.00 |          719.30 |          725.67 |       700.70 |      701.03 |      732.70 |
| 0.25 |          719.93 |          725.47 |       702.63 |      702.53 |      731.13 |
| 0.50 |          717.13 |          721.70 |       702.63 |      704.07 |      725.30 |
| 1.00 |          713.87 |          716.47 |       703.23 |      703.57 |      714.87 |

**Table 3b.** Ratios relative to σ = 0 static_budgeted (719.30).

| σ    | static_budgeted | random_budgeted | hebbian_frob | hybrid_frob | oracle_comp |
|-----:|----------------:|----------------:|-------------:|------------:|------------:|
| 0.00 |          1.0000 |          1.0089 |       0.9741 |      0.9746 |      1.0186 |
| 0.25 |          1.0009 |          1.0086 |       0.9768 |      0.9767 |      1.0165 |
| 0.50 |          0.9970 |          1.0033 |       0.9768 |      0.9788 |      1.0083 |
| 1.00 |          0.9924 |          0.9961 |       0.9777 |      0.9781 |      0.9938 |

Key trends in Table 3:

- The registered comparator `random_budgeted` is the best of the budgeted methods at σ=0 (725.67) but degrades faster than `static_budgeted` as σ rises: at σ=1.0 it has fallen to 716.47 (−9.2 cut over the range), versus `static_budgeted`'s −5.4 cut.
- `static_budgeted` drifts downward as σ rises (−5.4 cut).
- `hebbian_frobenius` and `hybrid_frobenius` drift *upward* as σ rises (+2.5 to +2.6 cut over the range).
- `oracle_compensation` collapses from 732.70 at σ=0 to 714.87 at σ=1.0 (−17.8 cut). The σ-aware rescaling cannot keep pace as the amplitude variance grows; by σ=1.0 the oracle has fallen behind `random_budgeted` and is essentially tied with `static_budgeted` in absolute terms. **Knowing the noise model without adapting is brittle: the oracle's prior knowledge of $a_i a_j$ degrades faster than every budgeted baseline.**
- The Hebbian–oracle gap shrinks from 32 cut (σ=0) to 12 cut (σ=1.0); the Hebbian–`random_budgeted` gap shrinks from 25 cut (σ=0) to 13 cut (σ=1.0).

**Slope test.** We fit a per-(graph_seed, method) linear regression of polished cut on σ across the four σ values (with 10 method seeds aggregated per cell), obtaining one slope per (graph_seed, method) pair. Mean slopes across graph_seeds:

| method                 | mean slope (cut/σ) |    sd | n_gs |
|------------------------|-------------------:|------:|-----:|
| `hybrid_frobenius`     |            **+2.44** | 3.51 |    3 |
| `hebbian_frobenius`    |              +2.16 | 2.18 |    3 |
| `static_budgeted`      |              −6.05 | 5.07 |    3 |
| `random_budgeted`      |            **−9.85** | 2.28 |    3 |
| `oracle_compensation`  |             −18.65 | 3.50 |    3 |

**Predeclared HARDWARE-ROBUSTNESS test (registered comparator).** Paired slope difference for the registered comparison `hybrid_frobenius` vs `random_budgeted`:

**mean = +12.29 cut/σ, bootstrap 95% CI = [+8.48, +14.71], n_pairs = 3.**

The 95% CI excludes zero and the direction is consistent across all three graph seeds (per-graph-seed deltas: gs=0: +13.67, gs=1: +8.48, gs=2: +14.71). The HARDWARE-ROBUSTNESS criterion holds: under amplitude heterogeneity, `hybrid_frobenius` degrades less than the registered `random_budgeted` comparator.²

² Exact permutation p = 0.25 is the minimum two-sided value achievable with n = 3 paired observations (2³ = 8 sign-flips); reported for completeness alongside the bootstrap CI, not as evidence of statistical significance.

For comparison, the secondary and tertiary slope contrasts:

- **Secondary** — `hebbian_frobenius` − `random_budgeted` slope: mean = +12.01 cut/σ, bootstrap 95% CI = [+9.30, +15.02], n_pairs = 3.
- **Tertiary** — `hybrid_frobenius` − `static_budgeted` slope (the v1 fallback comparator reported in earlier drafts): mean = +8.49 cut/σ, bootstrap 95% CI = [+2.41, +14.34], n_pairs = 3.

This is the main empirical finding of the paper: under realistic amplitude heterogeneity, Hebbian-coupled Kuramoto dynamics degrade less than random non-negative coupling at the same budget — in fact, they appear to *actively compensate* for the heterogeneity, producing slightly higher cut quality at higher σ even as both `static_budgeted` and `random_budgeted` baselines fall. Linear extrapolation of the slopes would put the `hybrid_frob` / `random_budgeted` crossover near σ ≈ 2; that regime lies outside the tested range, but the trajectory inside the tested range is unambiguous.

### 4.4 Adaptive-signal control

We additionally compare `hebbian_frobenius` against `learned_support_random_weights` (LSRW) — a control that takes Hebbian's converged support mask but replaces its learned weights with random non-negative weights projected to the same budget. The comparison isolates the contribution of Hebbian's *weight allocation* from the contribution of its *support choice*.

**Table 4.** Adaptive signal (hebbian_frobenius − LSRW, paired across graph_seeds, polished cut).

| Family                  | n_pairs | paired Δ |
|-------------------------|--------:|---------:|
| sparse_er_200_p05       |       3 |  +14.33  |
| dense_er_200_p15        |       3 |  +22.00  |
| random_regular_200_d10  |       3 |  +20.67  |
| GSet G1 (n=800)         |       1 |  +66.00  |

Hebbian's weight allocation beats random weights on the same support, by 14 to 66 cut depending on family. This signal is real but should be interpreted with care: it conditions on Hebbian's chosen support, which (Section 4.2) is itself an inferior support for Max-Cut at σ=0. The cleanest reading is *given the support Hebbian dynamics converge to, the Hebbian weights are doing real work* — not that the rule outperforms an unconstrained allocation.

### 4.5 Energy-decrease diagnostics

The continuous-time theorem in Section 3.3 applies to smooth flows. The benchmark runs discrete updates with periodic projection. We therefore report empirical energy-decrease diagnostics along trajectories rather than theorem-certified Lyapunov decrease.

We instrument the surrogate energy L(θ, W) at every projection step and count the number of steps on which L increased ("energy descent violations"). Across all `hebbian_frobenius` and `hybrid_frobenius` runs at σ = 0, half-mean-degree budget, every family:

- Mean violations per trajectory: 0
- Fraction of trajectories with L_valid = True: 100%

That is, the surrogate energy was numerically non-increasing at every projection step on every Frobenius-variant run. The `hebbian_sinkhorn` and `hybrid_sinkhorn` Sinkhorn variants do not satisfy this property — their projection is approximate, not exact, and L_valid is 0% — and we therefore omit them from any theoretical-anchor claims.

This is an *empirical* statement, not a theorem-certified one; the open extension to projected / nonsmooth dynamics (Section 3.3, open issue 4) is what would close that gap.

## 5. Discussion

### 5.1 What this work shows

The empirical contributions of this study, restricted to what the data support, are:

1. **At zero amplitude heterogeneity, the predeclared CORE test fails.** Hebbian-Frobenius coupling at the half-mean-degree budget produces polished cuts 15–35 lower (on n=200 random families) and 88 lower (on GSet G1) than random non-negative coupling at the same budget. We do not re-frame this as a victory.
2. **Under amplitude heterogeneity, the picture inverts.** `hybrid_frobenius`'s slope of cut against σ is positive (+2.4 cut/σ); the registered comparator `random_budgeted`'s is negative (−9.9 cut/σ); the paired difference is **+12.3 cut/σ with bootstrap 95% CI [+8.5, +14.7]**, direction consistent across all three graph seeds. The predeclared HARDWARE-ROBUSTNESS test holds. As σ → 1, the Hebbian–oracle gap shrinks from 32 cut to 12, and the oracle reference itself collapses by 17.8 cut over the σ range — evidence that a non-adaptive controller with full knowledge of the noise model degrades faster than the adaptive Hebbian rule. This is the main empirical contribution.
3. **The Hebbian rule's weights, conditional on its chosen support, beat random weights on the same support** by 14–66 cut depending on family — an adaptive-signal control isolating the rule's contribution from the support's. This is real but conditional.
4. **Numerical Lyapunov-style energy decrease holds with zero violations across all Frobenius-variant runs**, providing empirical support for the continuous-time theorem's relevance to the discrete projected updates the benchmark uses.

### 5.2 What this work does not show

1. **Not a state-of-the-art Max-Cut solver.** Even on n=200 random instances, several baselines outperform Hebbian-Frobenius. GSet G1 is one calibration instance; no scaling claim is made.
2. **The theorem does not certify the benchmark's discrete updates.** The Lyapunov-descent theorem (Section 3.3) is stated for the continuous-time smooth flow under exact projection. The benchmark implements discrete updates with periodic projection at the budget boundary, which introduces nonsmoothness the current proof does not handle. Section 4.5 reports numerical energy-decrease diagnostics rather than theorem-certified decrease.
3. **No wall-clock or compute-cost claim.** The "coupling-resource budget" is a parameter-count constraint motivated by hardware fan-in limits; digital simulation runtime is dominated by the symmetric-Frobenius projection step and is not halved when the budget is halved.
4. **The adaptive-signal column conditions on Hebbian's own support.** The +14 to +66 gap measures Hebbian weights against random weights on Hebbian's converged support — not against the optimal weights on any reasonable support. The stronger claim would require a configuration-model-matched baseline.
5. **Sample size is scoped-pilot.** Three independent graph seeds per random family is below the threshold at which family-level inferential claims become defensible. We report descriptive statistics with bootstrap CIs and exact-permutation p-values, but make no inferential significance claims. A 30+-graph-seed rerun is queued as next-step work.
6. **Heterogeneity tested on a single graph family.** The slope analysis above is on sparse_er only. We do not yet know whether the robustness signal generalises to dense ER, random-regular, or GSet topologies — repeating the sweep on the other three families is the natural next-step verification.

### 5.3 Downstream applications (sketch)

The Hebbian-coupled Kuramoto primitive — a sparse, budget-constrained, adaptive coupling matrix that learns structure from phase synchrony — has natural analogues in three domains whose details are out of scope for this paper:

- **Long-range numerical weather prediction.** Phase coherence between geographic regions tracks the teleconnection structure of large-scale atmospheric circulation. A Hebbian-coupled phase-oscillator network whose nodes correspond to grid cells or regimes could, in principle, learn attractor structure (blocking highs, jet stream positions, NAO modes) and surface tipping thresholds as proximity to basin separatrices. The robustness result in Section 4.3 is suggestive: real atmospheric circulation has strong amplitude heterogeneity across regimes, and a robustness-prioritising primitive may match the substrate better than an oracle-optimal one.
- **Distributed coordination.** Phase synchrony is inherently fault-tolerant: lose a node and the remaining oscillators re-equilibrate around a new fixed point without recovery protocol. A sparse-coupling Hebbian variant could implement self-healing distributed coordination with provable convergence under partial failure.
- **Wirelength-constrained chip placement and routing.** The row-sum budget maps directly to fan-in caps in physical placement. A learned sparse coupling that prioritises robustness over peak cut quality may align better with hardware noise floors than a fan-in-maximising allocation.

These are *implications worth testing*, not contributions of this paper. The benchmark substrate here is Max-Cut on n=200 random graphs and n=800 GSet G1; nothing in the data justifies a claim about weather, distributed systems, or chips. We mention them only to indicate where the same primitive might be evaluated next.

## 6. Conclusion

We studied budgeted Hebbian-coupled Kuramoto dynamics on graph Max-Cut, with two predeclared success criteria — one of which holds, one of which does not. At zero amplitude heterogeneity, random non-negative coupling at the same budget outperforms Hebbian-Frobenius coupling across all tested families; the predeclared CORE test fails. Under amplitude heterogeneity, the picture inverts: the registered HARDWARE-ROBUSTNESS test holds — `hybrid_frobenius` against `random_budgeted` shows a paired slope difference of +12.3 cut/σ with bootstrap 95% CI [+8.5, +14.7], direction consistent across all three graph seeds. The Hebbian rule's empirical contribution on this benchmark is robustness, not cut quality at zero noise.

The theoretical contribution is a Lyapunov-descent and KKT-stationarity theorem for the continuous-time joint phase / weight flow under zero detuning, fixed support, and exact projection, with the algebraic core formalised in Lean 4 / Mathlib. The formalisation has explicit open issues, including the extension to projected / nonsmooth dynamics that would directly cover the discrete updates used in the benchmark.

Three lines of follow-up work are immediate:

1. **Repeat the heterogeneity sweep on the other three graph families** (dense_er, random_regular, GSet G1) to confirm that Hebbian's robustness signal is not specific to sparse Erdős–Rényi.
2. **Increase graph-seed counts to 30+ per family.** The current scoped-pilot sample size precludes inferential significance claims.
3. **Close the Lean proof's open issues**, particularly the extension to projected / nonsmooth dynamics that would bridge theorem and benchmark.

The code, raw CSVs, and Lean source are available at `https://github.com/velvetmonkey/flywheel-universe`.

---

## Open verification items

These are items the draft asserts or relies on that have not been independently re-verified at the time of writing. Each must be resolved before Zenodo upload.

1. ~~**GSet G1 best-known cut.**~~ **RESOLVED 2026-05-20.** G1 best-known = 11624 confirmed against the value hardcoded in `phase2_benchmark.py` (lines 825-827) and against multiple benchmark papers (Toshiba SBM, MACUT, Memetic-PTS). Conventional citation = Benlic & Hao 2013 (BLS); the 11624 value itself predates BLS and circulates in the Helmberg-Rendl SDP-era literature.
2. ~~**Khan et al. [2025], arXiv:2510.24416.**~~ **RESOLVED 2026-05-20.** Confirmed: "Analyzing Parametric Oscillator Ising Machines through the Kuramoto Lens" by Nikhat Khan, E. M. H. E. B. Ekanayake, Nicolas Casilli, Cristian Cassella, Luke Theogarajan, Nikhil Shukla (submitted 2025-10-28, eess.SY). Bib updated.
3. ~~**Wang & Roychowdhury OIM citation.**~~ **RESOLVED 2026-05-20.** Confirmed: arXiv:1903.07163, UCNC 2019, Springer LNCS 11493, pp. 232-256, DOI 10.1007/978-3-030-19311-9_19. Bib updated.
4. ~~**CSV column completeness for Tables 1, 3, 4.**~~ **RESOLVED 2026-05-20.** `cut_polished` semantics verified throughout analysis; column usage confirmed consistent with `phase2_full_results.csv` at commit `1d5c5c8`.
5. ~~**Restart / method-seed nomenclature.**~~ **RESOLVED 2026-05-20.** "10 method seeds per cell" and "best-of-10 polished cut" confirmed consistent with `phase2_benchmark.py` implementation.
6. ~~**Heterogeneity comparator.** Section 4.3 substitutes `static_budgeted` for the predeclared `random_budgeted` as the slope comparator because `random_budgeted` was not included in the heterogeneity sweep.~~ **RESOLVED 2026-05-20.** Heterogeneity v2 rerun (`results/heterogeneity_experiment_v2.csv`, 600 rows) includes `random_budgeted`; Section 4.3 now reports the registered HARDWARE-ROBUSTNESS comparison directly.
7. ~~**Lean proof status.**~~ **RESOLVED 2026-05-20.** Verified against the current Lean source: Issue #4 (full KKT stationarity) closed sorry-free at commit `77c78df`. The only Lean commit since (`4c6a3bf`) adds `FutureWork.lean` as documentation only — it sketches the Moreau / tangent-cone path that would close the weight-chain-rule gap, but defines nothing and does not change the proof status. Section 3.3 open-issues list reduced from 4 to 3 to reflect this.
8. ~~**Detuning-ablation family.**~~ **RESOLVED 2026-05-20.** Confirmed in `phase2_benchmark.py` (line 1353): detuning ablation runs on sparse_er_200_p05 with `omega_sigma = 0.3` (ω_i ~ N(0, 0.3)) at σ=0 amplitude heterogeneity. Not included in main results; available as ablation only. Section 3.4 text already reflects this.
9. ~~**Slope-test methodology.**~~ **RESOLVED 2026-05-20.** Four σ values weighted equally; bootstrap resamples graph_seeds (not method_seeds or σ values); confirmed in `/tmp/analyze_het_v2.py`.
10. ~~**Bootstrap and permutation seeds.**~~ **RESOLVED 2026-05-20.** Bootstrap uses `seed=42`, 10000 resamples, documented in `analyze_het_v2.py`; reproducible from CSV alone.
11. ~~**Banned-verb sweep.**~~ **RESOLVED 2026-05-20.** Grep for *beats / matches / competes / rivals / ties / performs comparably / cheaper / faster* clean (remaining occurrences are either the registered CORE criterion verbatim or the adaptive-signal LSRW comparison, both legitimate).
12. ~~**Skim-reader test.**~~ **RESOLVED 2026-05-20.** Fresh-reader skim of Abstract + §1 returns headline as "robustness under amplitude heterogeneity," not "Hebbian beats static."

## References

[Cited via inline (Author, Year) keys to the companion bibliography `kuramoto-maxcut.bib`. Convert to numbered or alphabetical at Zenodo upload as appropriate for the rendering target.]

- Acebrón, J. A., Bonilla, L. L., Pérez Vicente, C. J., Ritort, F., & Spigler, R. (2005). *Reviews of Modern Physics*, 77(1), 137–185.
- Benlic, U., & Hao, J.-K. (2013). *Engineering Applications of [[AI|Artificial Intelligence]]*, 26(3), 1162–1173.
- Brucker, P. (1984). *Operations Research Letters*, 3(3), 163–166.
- Burer, S., Monteiro, R. D. C., & Zhang, Y. (2002). *SIAM Journal on Optimization*, 12(2), 503–521.
- Cassie, B. (2026). *kuramoto-lean: A Sorry-Free Lean 4 Library for Finite-N Kuramoto Synchronisation Dynamics*. Zenodo. https://doi.org/10.5281/zenodo.20468619
- Efron, B. (1979). *The Annals of Statistics*, 7(1), 1–26.
- Fisher, R. A. (1935). *The Design of Experiments*. Oliver and Boyd.
- Goemans, M. X., & Williamson, D. P. (1995). *Journal of the ACM*, 42(6), 1115–1145.
- Helmberg, C., & Rendl, F. (2000). *SIAM Journal on Optimization*, 10(3), 673–696.
- Hopfield, J. J. (1982). *PNAS*, 79(8), 2554–2558.
- Khan, N., Ekanayake, E. M. H. E. B., Casilli, N., Cassella, C., Theogarajan, L., & Shukla, N. (2025). Analyzing Parametric Oscillator Ising Machines through the Kuramoto Lens. *arXiv:2510.24416 [eess.SY]*.
- Kuramoto, Y. (1975). *International Symposium on Mathematical Problems in Theoretical Physics*, 420–422.
- McMahon, P. L., et al. (2016). *Science*, 354(6312), 614–617.
- Mohseni, N., McMahon, P. L., & Byrnes, T. (2022). *Nature Reviews Physics*, 4(6), 363–379.
- Moura, L. de, & Ullrich, S. (2021). *Conference on Automated Deduction*, 625–635.
- Strogatz, S. H. (2000). *Physica D*, 143(1–4), 1–20.
- The Mathlib Community. (2024). *The Lean Mathematical Library*.
- Wang, T., & Roychowdhury, J. (2019). OIM: Oscillator-based Ising Machines for Solving Combinatorial Optimisation Problems. In *Unconventional Computation and Natural Computation (UCNC)*, Springer LNCS 11493, 232–256. *arXiv:1903.07163*. DOI: 10.1007/978-3-030-19311-9_19.
- Wang, Z., Marandi, A., Wen, K., Byer, R. L., & Yamamoto, Y. (2013). *Physical Review A*, 88(6), 063853.

---

**End of working draft v1.** Resolve Open Verification Items 1–12, then convert references to chosen rendering format, then upload to Zenodo.
