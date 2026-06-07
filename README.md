# Hybrid Adaptive Calibration under Bounded Oscillator Coupling Budgets

[![Live demo](https://img.shields.io/badge/demo-live-d4af37?style=flat-square&logo=github)](https://velvetmonkey.github.io/flywheel-universe/)

The algebraic core of the descent argument has been machine-verified in Lean 4 / Mathlib ([lean/RequestProject/LyapunovDescent.lean](lean/RequestProject/LyapunovDescent.lean)). The formalization establishes convexity and closedness of the constraint set, the key symmetry-based phase contribution identity, and the descent conclusion under explicit hypotheses on the dL/dt decomposition and the projection's variational inequality. Closing the remaining gaps (weight dynamics, chain-rule derivation, projection existence, full KKT stationarity) is the next formalisation step — see open issues.

The formal foundations have been extended and published as a separate sorry-free library: kuramoto-lean (Ben Cassie, 2026), available at https://github.com/velvetmonkey/kuramoto-lean, with a companion paper at https://doi.org/10.5281/zenodo.20468619. That library proves the unprojected algebraic joint descent core in `hebbian_joint_lyapunov_descent` with zero sorry, zero admit, and no new axioms.

## 1. What this is

This repository studies budgeted Hebbian Kuramoto dynamics with a fixed sparsity support and symmetric-Frobenius projection. The model maintains a symmetric coupling matrix on a fixed edge mask, with per-node row-sum budgets enforced by exact projection at every step of an alternating phase / weight update. The contribution is a constrained descent flow on the joint state `(θ, W)` with an identified Lyapunov function, a hybrid freeze schedule that splits adaptive weight learning from fixed-weight phase settling, and a hardware-constraint framing in which the per-node budget represents the physical coupling-resource limit of an oscillator-based Ising machine. This is not a state-of-the-art Max-Cut solver, not a device-native learning law, and not a validated hardware primitive.

## Paper

**Budgeted Hebbian Kuramoto Dynamics for Max-Cut under Amplitude Heterogeneity: Robustness, Not Cut Quality, Is the Signal**
Ben Cassie (2026). Zenodo.
https://zenodo.org/records/20303914

## 2. The theorem (informal)

Under zero detuning (`ω = 0`), a fixed support mask, symmetric weights, and exact symmetric-Frobenius projection onto the budget polytope at each step, the joint dynamics `(θ̇, Ẇ)` descend a constrained surrogate energy. Limit points satisfy KKT stationarity for that energy under the row-budget, non-negativity, and support constraints.

## 3. What is tested

- **Graph families.** Sparse Erdős–Rényi (`n=200`, `p=0.05`), dense Erdős–Rényi (`n=200`, `p=0.15`), random 10-regular (`n=200`), small instances with proven optima via Gray-code enumeration (`n ∈ [20, 30]`), and GSet G1 (`n=800`) as a literature calibration anchor.
- **Methods compared.** Static Kuramoto on the original adjacency; static plus randomised-hyperplane rounding; greedy single-flip local search; static-projected, static-budgeted (top-k symmetric support), random-budgeted (random non-negative weights then projected), learned-support-random-weights (support inherited from a converged Hebbian run, random weights on top); topology-scrambled control (double-edge-swap mask); Hebbian and hybrid variants under both Sinkhorn-style and symmetric-Frobenius projection; SDP relaxation with Goemans–Williamson rounding.
- **Budget levels.** Row-sum budget set to `mean_degree` and to `0.5 · mean_degree`.
- **Amplitude heterogeneity sweep.** Lognormal node gains `a_i ∼ LogNormal(−σ²/2, σ)` with `σ ∈ {0.0, 0.25, 0.5, 1.0}`, effective coupling `K_eff[i,j] = a_i a_j W_ij`. This targets the amplitude-heterogeneity regime identified in Khan et al. (arXiv:2510.24416), where quasi-steady amplitude variations in parametric-oscillator Ising machines rescale the effective spin couplings and degrade solution quality. Methods compared: `static_budgeted`, `hebbian_frobenius`, `hybrid_frobenius`, and an oracle-compensation reference.

All experiments use 10 method seeds per cell. The benchmark tests whether the stationary coupling configuration produced by budgeted Hebbian adaptation is useful for graph optimisation under hardware-imposed coupling-resource constraints — specifically, whether the adaptive weight allocation provides signal beyond what the support sparsity and budget alone already provide.

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

## Interactive Demo

**Live:** <https://velvetmonkey.github.io/flywheel-universe/> — no install, no download.

A network of Kuramoto oscillators synchronises via Hebbian learning. Five topologies (fully connected, ring, sparse random, hub-spoke, island chain), speed and node-count controls, senescence and learning toggles, and a perturbation slider to find the basin boundary.

Source: [`demos/hebbian-kuramoto.html`](demos/hebbian-kuramoto.html) (also runs standalone in any browser).

---

Earlier framing of this repository — the universe / boundary-rider / cosmology-web analogies, the "six primitives" rhetoric, and the pre-Phase-2 Max-Cut numbers — is preserved in [`EXPLORATORY_NOTES.md`](EXPLORATORY_NOTES.md) and is not part of the main technical claim.

## Related repositories

Part of the **Flywheel suite** — local-first knowledge infrastructure over a plain-markdown Obsidian vault:

- [vault-core](https://github.com/velvetmonkey/vault-core) — Shared infrastructure for the Flywheel ecosystem.
- [flywheel-memory](https://github.com/velvetmonkey/flywheel-memory) — Persistent knowledge-graph memory MCP server: semantic search, read, and write over your vault.
- [flywheel-crank](https://github.com/velvetmonkey/flywheel-crank) — Desktop window into your vault's Flywheel MCP server.
- [flywheel-gravity](https://github.com/velvetmonkey/flywheel-gravity) — A compressed, reality-filtered context field over a vault.
- [flywheel-ideas](https://github.com/velvetmonkey/flywheel-ideas) — Local-first decision ledger: falsifiable bets, accepted outcomes, reusable lessons.
- [mega-monkey](https://github.com/velvetmonkey/mega-monkey) — Telegram-native AI research cockpit over an Obsidian vault.
- [roundtable](https://github.com/velvetmonkey/roundtable) — Local MCP server for delegating tasks to multiple AI models.

Research and experiments:

- [flywheel-concept](https://github.com/velvetmonkey/flywheel-concept) — A falsifiable study of cross-model concept geometry.
- [flywheel-geometry](https://github.com/velvetmonkey/flywheel-geometry) — A pre-registered study of cross-domain knowledge retrieval.
- **flywheel-universe** (this repo) — Lean 4 / Mathlib-verified core of the descent argument.
- [flywheel-velvetgram](https://github.com/velvetmonkey/flywheel-velvetgram) — Local widescreen Telegram reader for long-form reading.

Verified-cognition demo: [mcp-seal](https://github.com/velvetmonkey/mcp-seal) (verified MCP approval gate) and [canary](https://github.com/velvetmonkey/canary) (the seal demo host).
