# Hybrid Adaptive Calibration under Bounded Oscillator Coupling Budgets

[![Lyapunov Descent](https://img.shields.io/badge/theorem-Lean%204%20verified-brightgreen)](lean/RequestProject/LyapunovDescent.lean)

## 1. What this is

This repository studies budgeted Hebbian Kuramoto dynamics with a fixed sparsity support and symmetric-Frobenius projection. The model maintains a symmetric coupling matrix on a fixed edge mask, with per-node row-sum budgets enforced by exact projection at every step of an alternating phase / weight update. The contribution is a constrained descent flow on the joint state `(θ, W)` with an identified Lyapunov function, a hybrid freeze schedule that splits adaptive weight learning from fixed-weight phase settling, and a hardware-constraint framing in which the per-node budget represents the physical coupling-resource limit of an oscillator-based Ising machine. This is not a state-of-the-art Max-Cut solver, not a device-native learning law, and not a validated hardware primitive.

## 2. The theorem (informal)

Under zero detuning (`ω = 0`), a fixed support mask, symmetric weights, and exact symmetric-Frobenius projection onto the budget polytope at each step, the joint dynamics `(θ̇, Ẇ)` descend a constrained surrogate energy. Limit points satisfy KKT stationarity for that energy under the row-budget, non-negativity, and support constraints.

## 3. What is tested

- **Graph families.** Sparse Erdős–Rényi (`n=200`, `p=0.05`), dense Erdős–Rényi (`n=200`, `p=0.15`), random 10-regular (`n=200`), small instances with proven optima via Gray-code enumeration (`n ∈ [20, 30]`), and GSet G1 (`n=800`) as a literature calibration anchor.
- **Methods compared.** Static Kuramoto on the original adjacency; static plus randomised-hyperplane rounding; greedy single-flip local search; static-projected, static-budgeted (top-k symmetric support), random-budgeted (random non-negative weights then projected), learned-support-random-weights (support inherited from a converged Hebbian run, random weights on top); topology-scrambled control (double-edge-swap mask); Hebbian and hybrid variants under both Sinkhorn-style and symmetric-Frobenius projection; SDP relaxation with Goemans–Williamson rounding.
- **Budget levels.** Row-sum budget set to `mean_degree` and to `0.5 · mean_degree`.
- **Amplitude heterogeneity sweep.** Lognormal node gains `a_i ∼ LogNormal(−σ²/2, σ)` with `σ ∈ {0.0, 0.25, 0.5, 1.0}`, effective coupling `K_eff[i,j] = a_i a_j W_ij`, comparing `static_budgeted`, `hebbian_frobenius`, `hybrid_frobenius`, and an oracle compensation reference.

The benchmark tests whether the stationary coupling configuration produced by budgeted Hebbian adaptation is useful for graph optimisation under hardware-imposed coupling-resource constraints — specifically, whether the adaptive weight allocation provides signal beyond what the support sparsity and budget alone already provide.

## 4. What the results show

*Results from the full benchmark suite will be added here.*

## 5. What failed or is not claimed

- The Frobenius projection is implemented as a digital simulation proxy, not a hardware primitive. The dual root-finder satisfies KKT to machine epsilon on validation, but no physical realisation is claimed.
- The theorem requires zero detuning and a fixed mask. Detuning and adaptive topology are studied as ablations but lie outside the theorem.
- Classical Max-Cut solvers are not the target of this work. They appear as calibration baselines to anchor the cut-quality axis, not as competitive comparators.
- Biological plausibility of the quadratic weight decay term in the Lyapunov function is not claimed. The `λ‖W‖²_F / 4` regulariser is a control-theoretic ingredient, not a model of synaptic plasticity.
- This is a control / calibration algorithm. It is not yet a device-native learning law.

## 6. How to reproduce

```bash
git clone https://github.com/velvetmonkey/flywheel-universe
cd flywheel-universe
pip install -r requirements.txt
python validate_projection.py
python phase2_benchmark.py pilot
```

## 7. Citation / preprint

*Theory note and benchmark paper in preparation. Zenodo DOI will be added here.*

---

Earlier framing of this repository — the universe / boundary-rider / cosmology-web analogies, the "six primitives" rhetoric, and the pre-Phase-2 Max-Cut numbers — is preserved in [`EXPLORATORY_NOTES.md`](EXPLORATORY_NOTES.md) and is not part of the main technical claim.
