"""
Max-Cut benchmark for flywheel-universe Hebbian Kuramoto.
Tests budgeted Hebbian adaptive coupling vs static Kuramoto baseline.

Usage:
  python max_cut_benchmark.py          # uses random G(50,0.5) sanity check
  python max_cut_benchmark.py G1       # GSet G1 (download from Stanford GSet)

Defaults (T=50, n_runs=5) are tuned for fast sanity checks (~seconds).
For the real benchmark on GSet instances, call run_static / run_hebbian
with T=200, n_runs=20.
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

def hebbian_update_fast(W, theta, A_mask, K_sign=1.0, dt=0.5, eta=0.05, lam=0.003, budget=1.0):
    """Vectorised budgeted Hebbian update. Sign-aware: reinforces pairs at the
    equilibrium preferred by the coupling regime (in-phase for K>0, anti-phase
    for K<0). exist_mask enforced by A_mask."""
    cos_diff = np.cos(theta[np.newaxis, :] - theta[:, np.newaxis])
    W = W + dt * eta * (K_sign * cos_diff - lam * W)
    W = np.where(A_mask, np.maximum(W, 0.0), 0.0)
    row_totals = W.sum(axis=1) + 1e-12
    scale = np.minimum(1.0, budget / row_totals)
    return W * scale[:, np.newaxis]

def extract_cut(G, theta):
    nodes = sorted(G.nodes())
    assignment = {n: int(np.cos(theta[i]) <= 0) for i, n in enumerate(nodes)}
    return sum(G[u][v].get('weight', 1) for u, v in G.edges()
               if assignment[u] != assignment[v])

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
        best = max(best, extract_cut(G, theta))
    return best

# ── RUN ──
G = build_test_graph(n=200, p=0.05)   # sparse: 200 nodes, mean degree ~10
# G = build_test_graph(n=50, p=0.5)    # dense baseline (static AFM is hard to beat)
# G = load_gset('G1')                  # GSet G1 — 800 nodes, best known 11624

N, M = G.number_of_nodes(), G.number_of_edges()
mean_deg = 2 * M / N
print(f"Graph: {N} nodes, {M} edges, mean degree {mean_deg:.1f}\n")

print("Static Kuramoto (20 runs, K=-0.5, T=200)...")
t0 = time.time()
s = run_static(G, K=-0.5, T=200, n_runs=20)
print(f"  Best cut: {s}  ({time.time()-t0:.1f}s)\n")

print(f"Hebbian Kuramoto (20 runs, K=-0.5, T=200, budget=mean_deg={mean_deg:.1f})...")
t0 = time.time()
h = run_hebbian(G, K=-0.5, T=200, n_runs=20)
print(f"  Best cut: {h}  ({time.time()-t0:.1f}s)\n")

print(f"Static:  {s}")
print(f"Hebbian: {h}  ({'+' if h >= s else ''}{h-s} vs static)")
