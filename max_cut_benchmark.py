"""
Max-Cut benchmark for flywheel-universe Hebbian Kuramoto.
Tests budgeted Hebbian adaptive coupling vs static Kuramoto baseline.

Usage:
  python max_cut_benchmark.py          # uses random G(50,0.5) sanity check
  python max_cut_benchmark.py G1       # GSet G1 (download from Stanford GSet)
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

def build_adj(G):
    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    adj = [[] for _ in range(N)]
    wts = [[] for _ in range(N)]
    for u, v, d in G.edges(data=True):
        w = d.get('weight', 1)
        adj[idx[u]].append(idx[v]); wts[idx[u]].append(w)
        adj[idx[v]].append(idx[u]); wts[idx[v]].append(w)
    return adj, wts, N

def kuramoto_rhs(t, theta, omega, K, adj, weights):
    dtheta = omega.copy()
    for i in range(len(theta)):
        for j, w in zip(adj[i], weights[i]):
            dtheta[i] += K * w * np.sin(theta[j] - theta[i])
    return dtheta

def hebbian_update(W, theta, adj, dt=0.5, eta=0.05, lam=0.003, budget=1.0):
    """Budgeted Hebbian update — the novel piece. Exist-mask enforced by structure."""
    W_new = [list(row) for row in W]
    for i in range(len(theta)):
        for k, j in enumerate(adj[i]):
            W_new[i][k] += dt * eta * (np.cos(theta[j] - theta[i]) - lam * W[i][k])
            W_new[i][k] = max(W_new[i][k], 0.0)
        total = sum(W_new[i]) + 1e-12
        if total > budget:
            W_new[i] = [w * budget / total for w in W_new[i]]
    return W_new

def extract_cut(G, theta):
    nodes = sorted(G.nodes())
    assignment = {n: int(np.cos(theta[i]) <= 0) for i, n in enumerate(nodes)}
    cut = sum(G[u][v].get('weight', 1) for u, v in G.edges()
              if assignment[u] != assignment[v])
    return cut

def run_static(G, K=2.0, T=200, n_runs=20, seed=0):
    adj, wts, N = build_adj(G)
    rng = np.random.default_rng(seed)
    best = 0
    for _ in range(n_runs):
        omega = rng.normal(0, 0.3, N)
        theta0 = rng.uniform(0, 2*np.pi, N)
        sol = solve_ivp(kuramoto_rhs, (0, T), theta0,
                        args=(omega, K, adj, wts),
                        method='DOP853', rtol=1e-6, atol=1e-8)
        best = max(best, extract_cut(G, sol.y[:, -1]))
    return best

def run_hebbian(G, K=2.0, T=200, dt=0.5, n_runs=20,
                eta=0.05, lam=0.003, budget=1.0, seed=0):
    adj, wts, N = build_adj(G)
    rng = np.random.default_rng(seed)
    best = 0
    for _ in range(n_runs):
        omega = rng.normal(0, 0.3, N)
        theta = rng.uniform(0, 2*np.pi, N)
        W = [[budget / max(len(adj[i]), 1)] * len(adj[i]) for i in range(N)]
        t = 0.0
        while t < T:
            sol = solve_ivp(kuramoto_rhs, (t, t+dt), theta,
                            args=(omega, K, adj, W),
                            method='DOP853', rtol=1e-5, atol=1e-7)
            theta = sol.y[:, -1]
            t += dt
            W = hebbian_update(W, theta, adj, dt=dt, eta=eta, lam=lam, budget=budget)
        best = max(best, extract_cut(G, theta))
    return best

# ── RUN ──
G = build_test_graph(n=50, p=0.5)
# G = load_gset('G1')  # uncomment for real benchmark

N, M = G.number_of_nodes(), G.number_of_edges()
print(f"Graph: {N} nodes, {M} edges\n")

print("Static Kuramoto (20 runs)...")
t0 = time.time()
s = run_static(G)
print(f"  Best cut: {s}  ({time.time()-t0:.1f}s)\n")

print("Hebbian Kuramoto (20 runs)...")
t0 = time.time()
h = run_hebbian(G)
print(f"  Best cut: {h}  ({time.time()-t0:.1f}s)\n")

print(f"Static:  {s}")
print(f"Hebbian: {h}  ({'+' if h >= s else ''}{h-s} vs static)")
