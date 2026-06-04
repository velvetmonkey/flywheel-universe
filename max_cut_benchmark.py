"""
Max-Cut benchmark for flywheel-universe Hebbian Kuramoto.
Tests budgeted Hebbian adaptive coupling vs static Kuramoto baseline.

Usage:
  python max_cut_benchmark.py          # uses random G(50,0.5) sanity check
  python max_cut_benchmark.py G1       # GSet G1 (download from Stanford GSet)

Defaults (T=50, n_runs=5) are tuned for fast sanity checks (~seconds).
For the real benchmark on GSet instances, call run_static / run_hebbian
with T=200, n_runs=_n_runs.
"""

import numpy as np
from scipy.integrate import solve_ivp
import networkx as nx
import time

def build_test_graph(n=50, p=0.5, seed=42):
    return nx.erdos_renyi_graph(n, p, seed=seed)

def load_gset(path):
    with open(path) as f:
        n, m = map(int, f.readline().split())
        G = nx.Graph()
        G.add_nodes_from(range(n))
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                u, v = int(parts[0])-1, int(parts[1])-1
                w = int(parts[2]) if len(parts) >= 3 else 1
                G.add_edge(u, v, weight=w)
    return G

def graph_to_matrices(G):
    """Dense weighted adjacency + boolean exist-mask + N."""
    nodes = sorted(G.nodes())
    W = nx.to_numpy_array(G, nodelist=nodes, weight='weight')
    return W, (W > 0), len(nodes)

def kuramoto_rhs_fast(t, theta, omega, K, W):
    """Vectorised Kuramoto RHS — no Python loops."""
    return omega + K * (W * np.sin(theta[np.newaxis, :] - theta[:, np.newaxis])).sum(axis=1)

def enforce_symmetric_support_budget(W, A_mask, budget):
    """Symmetrize, re-apply support, then enforce the per-node budget."""
    W = 0.5 * (W + W.T)
    W = np.where(A_mask, np.maximum(W, 0.0), 0.0)
    np.fill_diagonal(W, 0.0)
    row_totals = W.sum(axis=1)
    if float(np.max(row_totals - budget)) > 1e-12:
        scale = np.where(row_totals > budget,
                         budget / np.maximum(row_totals, 1e-12),
                         1.0)
        W = W * np.minimum(scale[:, np.newaxis], scale[np.newaxis, :])
        W = np.where(A_mask, np.maximum(W, 0.0), 0.0)
        np.fill_diagonal(W, 0.0)
    return W

def hebbian_update_fast(W, theta, A_mask, K_sign=1.0, dt=0.5, eta=0.05, lam=0.003, budget=1.0):
    """Vectorised budgeted Hebbian update. Sign-aware: reinforces pairs at the
    equilibrium preferred by the coupling regime (in-phase for K>0, anti-phase
    for K<0). exist_mask enforced by A_mask."""
    cos_diff = np.cos(theta[np.newaxis, :] - theta[:, np.newaxis])
    W = W + dt * eta * (K_sign * cos_diff - lam * W)
    return enforce_symmetric_support_budget(W, A_mask, budget)

def extract_cut(G, theta):
    nodes = sorted(G.nodes())
    assignment = {n: int(np.cos(theta[i]) <= 0) for i, n in enumerate(nodes)}
    return sum(G[u][v].get('weight', 1) for u, v in G.edges()
               if assignment[u] != assignment[v])

def extract_cut_best(G, theta, n_hyperplanes=100, rng=None):
    """Best cut over n_hyperplanes random hyperplane roundings."""
    if rng is None:
        rng = np.random.default_rng()
    nodes = sorted(G.nodes())
    best = 0
    for _ in range(n_hyperplanes):
        angle = rng.uniform(0, 2 * np.pi)
        assignment = {n: int(np.cos(theta[i] - angle) <= 0)
                      for i, n in enumerate(nodes)}
        cut = sum(G[u][v].get('weight', 1) for u, v in G.edges()
                  if assignment[u] != assignment[v])
        best = max(best, cut)
    return best

def greedy_max_cut(W_adj, seed=0):
    """Greedy local-search Max-Cut on a symmetric weighted adjacency matrix.

    Starts at a random ±1 spin assignment; flips the single node with the
    largest cut-improving gain each step until no flip improves. For symmetric
    W with no self-loops, cut = 0.25 * (sum(W) - s^T W s); flipping s_k changes
    s^T W s by -4 s_k (Ws)_k, so the gain of flipping node k is +s_k (Ws)_k.
    Returns (cut_value, spins).
    """
    rng = np.random.default_rng(seed)
    N = W_adj.shape[0]
    s = rng.choice([-1, 1], size=N).astype(np.float64)
    while True:
        Ws = W_adj @ s
        gains = s * Ws
        k = int(np.argmax(gains))
        if gains[k] <= 0:
            break
        s[k] = -s[k]
    cut = 0.25 * (W_adj.sum() - s @ W_adj @ s)
    return float(cut), s.astype(int)

def run_static(G, K=2.0, T=50, n_runs=5, seed=0):
    W, _, N = graph_to_matrices(G)
    rng = np.random.default_rng(seed)
    best = 0
    for _ in range(n_runs):
        omega = rng.normal(0, 0.3, N)
        theta0 = rng.uniform(0, 2*np.pi, N)
        sol = solve_ivp(kuramoto_rhs_fast, (0, T), theta0,
                        args=(omega, K, W),
                        method='DOP853', rtol=1e-6, atol=1e-8)
        best = max(best, extract_cut(G, sol.y[:, -1]))
    return best

def run_hebbian(G, K=2.0, T=50, dt=0.5, n_runs=5,
                eta=0.05, lam=0.003, budget=None, seed=0):
    _, A_mask, N = graph_to_matrices(G)
    deg = A_mask.sum(axis=1, keepdims=True).clip(min=1)
    if budget is None:
        budget = float(deg.mean())  # match Static row-sum so coupling scales agree
    rng = np.random.default_rng(seed)
    best = 0
    for _ in range(n_runs):
        omega = rng.normal(0, 0.3, N)
        theta = rng.uniform(0, 2*np.pi, N)
        W = A_mask.astype(float) * (budget / deg)
        t = 0.0
        while t < T:
            sol = solve_ivp(kuramoto_rhs_fast, (t, t+dt), theta,
                            args=(omega, K, W),
                            method='DOP853', rtol=1e-5, atol=1e-7)
            theta = sol.y[:, -1]
            t += dt
            W = hebbian_update_fast(W, theta, A_mask, K_sign=np.sign(K), dt=dt, eta=eta, lam=lam, budget=budget)
        best = max(best, extract_cut_best(G, theta, rng=rng))
    return best

def run_hybrid(G, K=-0.5, T=200, T_prune_frac=0.3, dt=0.5, n_runs=20,
               eta=0.05, lam=0.003, budget=None, seed=0):
    """
    Phase 1 (T_prune_frac * T): Hebbian adaptive coupling learns sparse W.
    Phase 2 (remaining T): freeze W, run static Kuramoto on the learned graph.
    """
    _, A_mask, N = graph_to_matrices(G)
    if budget is None:
        budget = float(A_mask.sum(axis=1).mean())
    deg = A_mask.sum(axis=1, keepdims=True).clip(min=1)
    rng = np.random.default_rng(seed)
    best = 0
    T_prune = T * T_prune_frac
    T_static = T * (1 - T_prune_frac)
    for _ in range(n_runs):
        omega = rng.normal(0, 0.3, N)
        theta = rng.uniform(0, 2 * np.pi, N)
        W = A_mask.astype(float) * (budget / deg)
        # Phase 1: adaptive pruning
        t = 0.0
        while t < T_prune:
            sol = solve_ivp(kuramoto_rhs_fast, (t, t+dt), theta,
                            args=(omega, K, W),
                            method='DOP853', rtol=1e-5, atol=1e-7)
            theta = sol.y[:, -1]
            t += dt
            W = hebbian_update_fast(W, theta, A_mask, K_sign=np.sign(K),
                                    dt=dt, eta=eta, lam=lam, budget=budget)
        # Phase 2: freeze W, run static on the learned sparse graph
        sol = solve_ivp(kuramoto_rhs_fast, (0, T_static), theta,
                        args=(omega, K, W),
                        method='DOP853', rtol=1e-6, atol=1e-8)
        theta = sol.y[:, -1]
        best = max(best, extract_cut_best(G, theta, rng=rng))
    return best

# ── RUN ──
import os, csv

def run_method(method, G, seed, K=-0.5, T=200):
    """One method × one seed on G. Returns (cut, runtime_seconds)."""
    t0 = time.time()
    if method == 'static':
        cut = run_static(G, K=K, T=T, n_runs=1, seed=seed)
    elif method == 'random_rounding':
        W, _, N = graph_to_matrices(G)
        rng = np.random.default_rng(seed)
        omega = rng.normal(0, 0.3, N)
        theta0 = rng.uniform(0, 2*np.pi, N)
        sol = solve_ivp(kuramoto_rhs_fast, (0, T), theta0,
                        args=(omega, K, W),
                        method='DOP853', rtol=1e-6, atol=1e-8)
        cut = extract_cut_best(G, sol.y[:, -1], rng=rng)
    elif method == 'hebbian':
        cut = run_hebbian(G, K=K, T=T, n_runs=1, seed=seed)
    elif method == 'hybrid':
        cut = run_hybrid(G, K=K, T=T, n_runs=1, seed=seed)
    elif method == 'greedy':
        W, _, _ = graph_to_matrices(G)
        cut, _ = greedy_max_cut(W, seed=seed)
    else:
        raise ValueError(method)
    return float(cut), time.time() - t0

graphs = {
    'sparse_er_200_p05': build_test_graph(n=200, p=0.05),
    'dense_er_200_p15':  build_test_graph(n=200, p=0.15),
}
if os.path.exists('G1'):
    graphs['gset_g1'] = load_gset('G1')
else:
    print("G1 not found in repo root — skipping GSet benchmark")

methods = ['static', 'random_rounding', 'hebbian', 'hybrid', 'greedy']
seeds = list(range(10))

rows = []
for graph_name, G in graphs.items():
    print(f"\n=== {graph_name}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges ===")
    for method in methods:
        cuts, times = [], []
        for seed in seeds:
            cut, rt = run_method(method, G, seed)
            cuts.append(cut)
            times.append(rt)
            print(f"  {method:18s} seed={seed} cut={cut:.0f} ({rt:.1f}s)")
        rows.append({
            'graph': graph_name,
            'method': method,
            'mean_cut': float(np.mean(cuts)),
            'max_cut':  float(np.max(cuts)),
            'std_cut':  float(np.std(cuts)),
            'mean_runtime': float(np.mean(times)),
        })

os.makedirs('results', exist_ok=True)
csv_path = 'results/benchmark_results.csv'
with open(csv_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['graph', 'method', 'mean_cut', 'max_cut', 'std_cut', 'mean_runtime'])
    w.writeheader()
    w.writerows(rows)

print(f"\nWrote {csv_path}\n")
print(f"{'graph':22s} {'method':18s} {'mean':>10s} {'max':>8s} {'std':>8s} {'time(s)':>10s}")
for r in rows:
    print(f"{r['graph']:22s} {r['method']:18s} {r['mean_cut']:10.1f} {r['max_cut']:8.0f} {r['std_cut']:8.2f} {r['mean_runtime']:10.1f}")
