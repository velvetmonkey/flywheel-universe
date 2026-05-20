"""Heterogeneity experiment v2 — adds `random_budgeted` to the method list so the
predeclared HARDWARE-ROBUSTNESS comparator (hybrid_frobenius vs random_budgeted
slope difference) can be computed against the registered baseline.

Identical protocol to `run_heterogeneity_experiment()` in phase2_benchmark.py
otherwise: sparse_er_200_p05, σ ∈ {0.0, 0.25, 0.5, 1.0}, 3 graph seeds × 10 method seeds.

Writes results to `results/heterogeneity_experiment_v2.csv`.
Does NOT modify phase2_benchmark.py — preserves the registered v1 protocol.
"""
import os
import sys
import time

import numpy as np
import networkx as nx
import pandas as pd
from scipy.integrate import solve_ivp

# Import everything we need from phase2_benchmark
from phase2_benchmark import (
    make_frobenius_projector,
    project_frobenius_dual,
    kuramoto_rhs_fast,
    random_hyperplane_rounding,
    make_rngs,
    run_method,
    precompute_learned_support,
)


HET_METHODS_V2 = [
    'static_budgeted',
    'random_budgeted',          # <-- newly added for v2
    'hebbian_frobenius',
    'hybrid_frobenius',
    'oracle_compensation',
]
SIGMAS = [0.0, 0.25, 0.5, 1.0]
N_GRAPH_SEEDS = 3
N_METHOD_SEEDS = 10
N = 200


def run_heterogeneity_v2():
    print("\n--- HETEROGENEITY EXPERIMENT V2 (sparse_er, lognormal gains) ---")
    print(f"Methods: {HET_METHODS_V2}")
    print(f"Sigmas: {SIGMAS}")
    print(f"Graph seeds: 0..{N_GRAPH_SEEDS-1}, Method seeds: 0..{N_METHOD_SEEDS-1}")
    print(f"Expected total cells: {N_GRAPH_SEEDS * len(SIGMAS) * N_METHOD_SEEDS * len(HET_METHODS_V2)}")
    sys.stdout.flush()

    rows = []
    t_start = time.time()

    for g_seed in range(N_GRAPH_SEEDS):
        G = nx.erdos_renyi_graph(N, 0.05, seed=g_seed)
        A = nx.to_numpy_array(G)
        np.fill_diagonal(A, 0)
        A_mask = (A > 0)
        m = G.number_of_edges()
        mean_deg = 2.0 * m / N
        budget_value = mean_deg / 2.0
        prob, Y_param, X_var = make_frobenius_projector(N, A_mask, budget_value)

        learned_mask, _, _, _, _ = precompute_learned_support(
            G, A, mean_deg, budget_value, A_mask, N)

        print(f"\nGraph seed {g_seed}: n={N}, m={m}, mean_deg={mean_deg:.2f}, budget={budget_value:.2f}")
        sys.stdout.flush()

        for sigma in SIGMAS:
            for m_seed in range(N_METHOD_SEEDS):
                rngs = make_rngs(g_seed, m_seed)
                theta0 = rngs['phases'].uniform(0, 2 * np.pi, N)
                a_gains = rngs['amplitude'].lognormal(-0.5 * sigma**2, sigma, N)
                omega = np.zeros(N)

                for method in HET_METHODS_V2:
                    t0 = time.time()
                    if method == 'oracle_compensation':
                        a_outer = np.outer(a_gains, a_gains)
                        W_oracle = np.where(A_mask, A.astype(float) / (a_outer + 1e-12), 0.0)
                        W_oracle = 0.5 * (W_oracle + W_oracle.T)
                        W_final = project_frobenius_dual(W_oracle, A_mask, budget_value)
                        W_eff = W_final * a_outer
                        sol = solve_ivp(
                            kuramoto_rhs_fast, (0, 200), theta0,
                            args=(omega, -0.5, W_eff),
                            method='DOP853', rtol=1e-6, atol=1e-8)
                        if not sol.success:
                            raise RuntimeError(f"solve_ivp failed: {sol.message}")
                        cut_raw, cut_rounded, cut_polished = random_hyperplane_rounding(
                            sol.y[:, -1], A, seeds=20, rng=rngs['rounding'])
                        corr_W_inv_gain = 1.0
                        corr_W_gain = -1.0
                    else:
                        res = run_method(
                            method, A, theta0, omega, a_gains, budget_value,
                            prob, Y_param, X_var, rngs,
                            learned_support_mask=learned_mask)
                        cut_raw, cut_rounded, cut_polished = res[0], res[1], res[2]
                        W_final = res[-1]
                        if W_final is not None and method in ('hebbian_frobenius', 'hybrid_frobenius'):
                            active = (W_final > 1e-6) & np.triu(np.ones_like(W_final, dtype=bool), 1)
                            if active.any():
                                Wv = W_final[active]
                                a_outer = np.outer(a_gains, a_gains)
                                inv_gain = (1.0 / (a_outer + 1e-12))[active]
                                gain = a_outer[active]
                                if len(Wv) > 1 and np.std(Wv) > 1e-12:
                                    corr_W_inv_gain = float(np.corrcoef(Wv, inv_gain)[0, 1])
                                    corr_W_gain = float(np.corrcoef(Wv, gain)[0, 1])
                                else:
                                    corr_W_inv_gain = corr_W_gain = 0.0
                            else:
                                corr_W_inv_gain = corr_W_gain = 0.0
                        else:
                            corr_W_inv_gain = corr_W_gain = 0.0

                    rows.append({
                        'sigma': sigma, 'method': method,
                        'graph_seed': g_seed, 'method_seed': m_seed,
                        'cut_raw': cut_raw, 'cut_rounded': cut_rounded,
                        'cut_polished': cut_polished,
                        'corr_W_inverse_gain': corr_W_inv_gain,
                        'corr_W_gain': corr_W_gain,
                        'runtime_sec': time.time() - t0,
                    })
                    elapsed = time.time() - t_start
                    print(f"    g={g_seed} m={m_seed} sigma={sigma} {method:30s} "
                          f"cut={cut_polished:.0f} t={time.time()-t0:.1f}s "
                          f"[total elapsed {elapsed/60:.1f}min, n_rows={len(rows)}]")
                    sys.stdout.flush()

    df = pd.DataFrame(rows)
    baseline = (df[(df['method'] == 'static_budgeted') & (df['sigma'] == 0.0)]
                .groupby('graph_seed')['cut_polished'].mean().to_dict())
    df['cut_ratio_vs_baseline'] = df.apply(
        lambda r: r['cut_polished'] / baseline.get(r['graph_seed'], 1.0)
        if baseline.get(r['graph_seed'], 0.0) > 0 else 0.0,
        axis=1)

    out_path = 'results/heterogeneity_experiment_v2.csv'
    df.to_csv(out_path, index=False)
    total_min = (time.time() - t_start) / 60.0
    print(f"\nWrote {out_path} ({len(rows)} rows) in {total_min:.1f} min")
    return df


if __name__ == "__main__":
    os.makedirs('results', exist_ok=True)
    run_heterogeneity_v2()
