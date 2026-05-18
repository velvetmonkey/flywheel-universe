"""
Phase 2 Max-Cut benchmark — symmetric-Frobenius projection, Lyapunov certification,
multi-baseline comparison, amplitude heterogeneity ablation.

Predeclared success criteria (amendment 2):

  CORE
    Setup: sigma = 0, budget = half_mean_degree.
    Claim: hybrid_frobenius beats random_budgeted on polished cut, aggregated
           best-of-restarts per graph_seed, paired Wilcoxon across graph_seeds.

  HARDWARE-ROBUSTNESS
    Setup: amplitude sigma >= 0.5 (heterogeneity block).
    Claim: hybrid_frobenius degrades less than random_budgeted, and tracks
           oracle_compensation more closely as sigma rises.

Statistical claim unit (amendment 1): the graph instance, not the restart.
Aggregate restarts within a (family, graph_seed, method) cell to best-of-r
polished cut, then paired Wilcoxon across graph_seeds within a family.
G1 is calibration only, not a primary proof point.

Self-contained; does NOT import from Phase 1 (max_cut_benchmark.py). Run:

    python phase2_benchmark.py pilot   # restricted pilot — Amendment A pass gates
    python phase2_benchmark.py         # full Phase 2 (long, ~6h)
"""

import os
import sys
import time
import itertools
import numpy as np
import networkx as nx
import pandas as pd
import cvxpy as cp
from scipy.integrate import solve_ivp

# Module-level CVXPY-projection timing accumulator. Reset between probe windows.
_PROJ_STATS = {'time_sec': 0.0, 'count': 0}

def _proj_reset():
    _PROJ_STATS['time_sec'] = 0.0
    _PROJ_STATS['count'] = 0

def _proj_snapshot():
    return _PROJ_STATS['time_sec'], _PROJ_STATS['count']

# =============================================================================
# 0. PHASE 1 LOADER (G1 GSet parser, kept locally so this file is standalone)
# =============================================================================

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

# =============================================================================
# 1. PROJECTIONS & UNIT TESTS
# =============================================================================

def make_frobenius_projector(n, A_mask, budget):
    """Compile the CVXPY symmetric-Frobenius projection (validator only) and stash
    (A_mask, budget) on the prob so the dual solver can pick them up at runtime.

    Fix 1 — uses cp.multiply(1 - A_mask.astype(float), X) == 0 instead of indexing
    X[~A_mask], which is fragile for cvxpy expressions.
    """
    Y_param = cp.Parameter((n, n), symmetric=True)
    X = cp.Variable((n, n), symmetric=True)
    A_off = (1.0 - A_mask.astype(float))
    obj = cp.Minimize(cp.sum_squares(X - Y_param))
    constraints = [
        X >= 0,
        cp.multiply(A_off, X) == 0,
        cp.diag(X) == 0,
        cp.sum(X, axis=1) <= budget,
    ]
    prob = cp.Problem(obj, constraints)
    prob._A_mask = A_mask.copy()
    prob._budget = float(budget)
    return prob, Y_param, X

# -----------------------------------------------------------------------------
# Dual root-finder projection (runtime path)
#
# KKT form: X_ij = max(0, Y_ij - lambda_i - lambda_j) on the active mask, with
# lambda_i >= 0 and complementary slackness on the row-budget constraint
# (sum_j X_ij <= budget). Single-row sub-problem is the classical projection
# onto the budget-simplex (Brucker, O(k log k) via sort). Outer Gauss-Seidel
# sweep over rows until ||Δλ||_∞ < tol.
# -----------------------------------------------------------------------------

def _solve_lambda_row(offs, budget):
    """Find lambda >= 0 such that sum_j max(0, offs[j] - lambda) <= budget,
    with complementary slackness. offs is 1D float array of (Y_ij - lam_j) over
    j active. Returns the unique lambda.
    """
    if offs.size == 0:
        return 0.0
    pos = offs[offs > 0]
    if pos.size == 0 or pos.sum() <= budget + 1e-14:
        return 0.0
    s = np.sort(pos)[::-1]
    cum = np.cumsum(s)
    k = np.arange(1, len(s) + 1, dtype=np.float64)
    # Brucker: rho = largest k with k * s[k-1] > cum[k-1] - budget
    cond = k * s > (cum - budget)
    idxs = np.where(cond)[0]
    if idxs.size == 0:
        return 0.0
    rho = int(idxs[-1])  # 0-indexed
    lam = (cum[rho] - budget) / (rho + 1)
    return float(max(0.0, lam))

def project_frobenius_dual(Y, A_mask, budget, max_outer=300, tol=1e-12,
                           init_lam=None, return_lam=False):
    """Project Y onto {symmetric X : X >= 0, X[~A_mask]=0, diag(X)=0, row_sum <= budget}.

    Dual root-finder via Gauss-Seidel block coordinate descent on lambda_i.
    init_lam: warm-start vector (huge speedup when Y changes slowly between
              calls — e.g. inside the Hebbian staggered loop).
    return_lam: if True, also returns the converged lambda for warm-starting
                the next call.
    """
    n = Y.shape[0]
    Y = 0.5 * (Y + Y.T).astype(np.float64)
    mask = (A_mask & ~np.eye(n, dtype=bool))
    Y = np.where(mask, Y, 0.0)

    if init_lam is not None and init_lam.shape == (n,):
        lam = init_lam.astype(np.float64).copy()
    else:
        lam = np.zeros(n)

    for outer in range(max_outer):
        lam_prev = lam.copy()
        for i in range(n):
            active = mask[i]
            if not active.any():
                lam[i] = 0.0
                continue
            offs = Y[i, active] - lam[active]
            lam[i] = _solve_lambda_row(offs, budget)
        if np.max(np.abs(lam - lam_prev)) < tol:
            break

    D = Y - lam[:, None] - lam[None, :]
    X = np.where(mask, np.maximum(D, 0.0), 0.0)
    # Symmetric by construction (Y, mask, lam_i+lam_j all symmetric); guard for fp noise
    X = 0.5 * (X + X.T)
    X = np.where(mask, np.maximum(X, 0.0), 0.0)

    # Strict budget enforcement: GS converges to budget within finite tol; after
    # 100s of staggered calls residual overshoot can reach ~1e-5. Symmetric
    # min(row, col) factor scales any overshoot back to feasibility while
    # preserving symmetry.
    row_sums = X.sum(axis=1)
    overshoot = row_sums - budget
    if float(np.max(overshoot)) > 1e-12:
        factor = np.where(row_sums > budget,
                          budget / np.maximum(row_sums, 1e-12),
                          1.0)
        # Symmetric scaling: scale entry (i,j) by min(f_i, f_j) preserves symmetry
        sym_factor = np.minimum(factor[:, None], factor[None, :])
        X = X * sym_factor
        X = np.where(mask, np.maximum(X, 0.0), 0.0)

    if return_lam:
        return X, lam
    return X

def _project_frobenius_cvxpy(W, prob, Y_param, X_var, first_call=False):
    """CVXPY-based projection — used only inside unit_test_projection as the
    reference solver for validating project_frobenius_dual."""
    Y_param.value = 0.5 * (W + W.T)
    try:
        prob.solve(solver=cp.OSQP, warm_start=(not first_call),
                   eps_abs=1e-9, eps_rel=1e-9, max_iter=50000, verbose=False)
    except Exception:
        prob.solve(solver=cp.SCS, warm_start=False, verbose=False)
    if X_var.value is None or prob.status not in ('optimal', 'optimal_inaccurate'):
        prob.solve(solver=cp.SCS, warm_start=False, verbose=False,
                   eps=1e-9, max_iters=50000)
    return X_var.value

def project_frobenius(W, prob, Y_param, X_var, first_call=False):
    """Runtime projection: dispatches to the dual root-finder.

    Signature unchanged for backward compatibility — prob, Y_param, X_var are
    ignored at runtime; A_mask and budget are pulled from attributes stashed on
    prob by make_frobenius_projector. Warm-starts lambda from prob._lam_cache
    across calls; first_call=True clears the cache (fresh-method start).
    Wall-time accumulated in _PROJ_STATS.
    """
    t0 = time.time()
    if first_call:
        prob._lam_cache = None
    init_lam = getattr(prob, '_lam_cache', None)
    X, prob._lam_cache = project_frobenius_dual(
        W, prob._A_mask, prob._budget,
        init_lam=init_lam, return_lam=True)
    _PROJ_STATS['time_sec'] += time.time() - t0
    _PROJ_STATS['count'] += 1
    return X

def project_sinkhorn(W, budget, A_mask, max_iters=100, tol=1e-9):
    """Sinkhorn-style alternating projection (named ablation).

    Tighter defaults (was max_iters=50, tol=1e-6) plus a final symmetric
    budget-enforcement clip so residual row-sum overshoot lands at machine
    epsilon rather than the alternating-projection tolerance.
    """
    W = 0.5 * (W + W.T)
    W = np.clip(W, 0, None)
    W[~A_mask] = 0
    np.fill_diagonal(W, 0)
    for _ in range(max_iters):
        row_sums = W.sum(axis=1).copy()
        row_sums[row_sums < 1e-12] = 1.0
        scale = np.minimum(1.0, budget / row_sums)
        W_new = W * scale[:, None]
        W_new = 0.5 * (W_new + W_new.T)
        W_new[~A_mask] = 0
        np.fill_diagonal(W_new, 0)
        if np.max(np.abs(W - W_new)) < tol:
            break
        W = W_new

    # Strict budget enforcement (same as project_frobenius_dual)
    row_sums = W.sum(axis=1)
    overshoot = row_sums - budget
    if float(np.max(overshoot)) > 1e-12:
        factor = np.where(row_sums > budget,
                          budget / np.maximum(row_sums, 1e-12),
                          1.0)
        sym_factor = np.minimum(factor[:, None], factor[None, :])
        W = W * sym_factor
        W[~A_mask] = 0
        np.fill_diagonal(W, 0)
    return W

def unit_test_projection(n_problems=20, n=10):
    """Validate project_frobenius_dual against CVXPY across N random problems at n=10.

    Fix 2 — seven explicit KKT checks on the dual output: symmetry, budget,
    mask, max-off-mask weight, non-negativity, idempotence, zero diagonal.
    Plus the cross-solver agreement check ||X_dual - X_cvxpy||_F.
    """
    print(f"[unit test] dual vs CVXPY at n={n}, {n_problems} random problems...")
    rng = np.random.default_rng(42)
    max_diff = 0.0
    worst_checks = {}
    for trial in range(n_problems):
        budget = float(rng.uniform(1.0, 5.0))
        density = float(rng.uniform(0.3, 0.8))
        A_mask = rng.random((n, n)) < density
        A_mask = A_mask | A_mask.T
        np.fill_diagonal(A_mask, 0)
        prob, Y_param, X_var = make_frobenius_projector(n, A_mask, budget)
        W_raw = rng.normal(0, 1, (n, n))

        X_dual = project_frobenius_dual(W_raw, A_mask, budget)
        X_cvxpy = _project_frobenius_cvxpy(W_raw, prob, Y_param, X_var, first_call=True)

        diff = float(np.linalg.norm(X_dual - X_cvxpy, 'fro'))
        sym_err = float(np.linalg.norm(X_dual - X_dual.T, 'fro') / (np.linalg.norm(X_dual, 'fro') + 1e-12))
        bud_viol = float(np.max(X_dual.sum(axis=1) - budget))
        mask_viol = float(np.max(np.abs(X_dual[~A_mask])))
        max_off_mask_weight = mask_viol
        nonneg_min = float(np.min(X_dual))
        diag_max = float(np.max(np.abs(np.diag(X_dual))))
        X_dual2 = project_frobenius_dual(X_dual, A_mask, budget)
        idem = float(np.linalg.norm(X_dual2 - X_dual, 'fro'))

        if diff > max_diff:
            max_diff = diff
            worst_checks = dict(trial=trial, diff=diff, sym_err=sym_err, bud_viol=bud_viol,
                                mask_viol=mask_viol, nonneg_min=nonneg_min, idem=idem, diag_max=diag_max)

        fails = []
        if diff > 1e-5: fails.append(f'cvxpy_diff={diff:.3e}')
        if sym_err > 1e-6: fails.append(f'sym_err={sym_err:.3e}')
        if bud_viol > 1e-6: fails.append(f'bud_viol={bud_viol:.3e}')
        if mask_viol > 1e-6: fails.append(f'mask_viol={mask_viol:.3e}')
        if max_off_mask_weight > 1e-6: fails.append(f'max_off_mask_weight={max_off_mask_weight:.3e}')
        if nonneg_min < -1e-6: fails.append(f'nonneg_min={nonneg_min:.3e}')
        if idem > 1e-6: fails.append(f'idem={idem:.3e}')
        if diag_max > 1e-6: fails.append(f'diag_max={diag_max:.3e}')
        if fails:
            raise RuntimeError(f"Projection unit test FAILED trial {trial}: " + "; ".join(fails))
    print(f"✓ Projection unit test passed: {n_problems} problems, "
          f"max ||X_dual-X_cvxpy||_F = {max_diff:.3e}; worst-case = {worst_checks}")
    return {'max_diff': max_diff, 'worst': worst_checks}

def profile_dual_at_n200(n_calls=50, target_ms=50.0):
    """Median wall time of project_frobenius_dual at n=200 over `n_calls` random
    problems (mix of sparse and dense ER). Raises if median > target_ms."""
    print(f"[profile] project_frobenius_dual at n=200, {n_calls} calls...")
    rng = np.random.default_rng(0)
    n = 200
    times = []
    for c in range(n_calls):
        p = float(rng.uniform(0.03, 0.20))
        G = nx.erdos_renyi_graph(n, p, seed=int(rng.integers(0, 10_000_000)))
        A = nx.to_numpy_array(G)
        A_mask = (A > 0) & ~np.eye(n, dtype=bool)
        m = G.number_of_edges()
        if m == 0:
            continue
        mean_deg = 2.0 * m / n
        budget = mean_deg / 2.0
        Y = rng.normal(0, 1, (n, n))
        t0 = time.time()
        _ = project_frobenius_dual(Y, A_mask, budget)
        times.append(time.time() - t0)
    times_ms = np.array(times) * 1000.0
    median_ms = float(np.median(times_ms))
    mean_ms = float(np.mean(times_ms))
    p95_ms = float(np.percentile(times_ms, 95))
    max_ms = float(np.max(times_ms))
    print(f"  median: {median_ms:.1f}ms   mean: {mean_ms:.1f}ms   p95: {p95_ms:.1f}ms   max: {max_ms:.1f}ms")
    if median_ms > target_ms:
        raise RuntimeError(f"Dual solver too slow: median {median_ms:.1f}ms > target {target_ms:.1f}ms")
    print(f"✓ Profile passed: median {median_ms:.1f}ms < {target_ms}ms")
    return dict(median_ms=median_ms, mean_ms=mean_ms, p95_ms=p95_ms, max_ms=max_ms)

# =============================================================================
# 2. METRICS & BASELINES
# =============================================================================

def get_cut_value(A, spins):
    return 0.25 * np.sum(A * (1 - np.outer(spins, spins)))

def greedy_max_cut(A, spins):
    """1-flip local search starting from given spins."""
    spins = spins.copy().astype(np.float64)
    n = len(spins)
    improved = True
    while improved:
        improved = False
        As = A @ spins
        gains = spins * As
        k = int(np.argmax(gains))
        if gains[k] > 0:
            spins[k] = -spins[k]
            improved = True
    return spins.astype(int)

def random_hyperplane_rounding(theta, A, seeds=20, rng=None):
    """Returns (raw, best rounded, polished) cut values.

    rng is an np.random.Generator (typically a spawned sub-stream).
    """
    if rng is None:
        rng = np.random.default_rng()
    spins_raw = np.sign(np.cos(theta))
    spins_raw[spins_raw == 0] = 1
    cut_raw = get_cut_value(A, spins_raw)

    best_rounded = cut_raw
    best_spins = spins_raw.copy()
    vecs = np.vstack([np.cos(theta), np.sin(theta)]).T
    for _ in range(seeds):
        r = rng.normal(0, 1, 2)
        spins = np.sign(vecs @ r)
        spins[spins == 0] = 1
        cut = get_cut_value(A, spins)
        if cut > best_rounded:
            best_rounded = cut
            best_spins = spins

    polished_spins = greedy_max_cut(A, best_spins)
    cut_polished = get_cut_value(A, polished_spins)
    return float(cut_raw), float(best_rounded), float(cut_polished)

def run_sdp_gw(A, rng):
    """SDP relaxation + Goemans-Williamson rounding. CLARABEL primary, SCS fallback."""
    n = A.shape[0]
    X = cp.Variable((n, n), PSD=True)
    obj = cp.Maximize(0.25 * cp.sum(cp.multiply(A, 1 - X)))
    constraints = [cp.diag(X) == 1]
    prob = cp.Problem(obj, constraints)
    try:
        prob.solve(solver=cp.CLARABEL, verbose=False)
    except Exception:
        prob.solve(solver=cp.SCS, verbose=False)
    if X.value is None or prob.status not in ('optimal', 'optimal_inaccurate'):
        try:
            prob.solve(solver=cp.SCS, verbose=False, eps=1e-8, max_iters=100000)
        except Exception:
            return 0.0, 0.0
    if X.value is None:
        return 0.0, 0.0
    sdp_bound = prob.value
    X_val = X.value
    try:
        L = np.linalg.cholesky(X_val + 1e-9 * np.eye(n))
    except np.linalg.LinAlgError:
        e, v = np.linalg.eigh(X_val)
        e[e < 0] = 0
        L = v @ np.diag(np.sqrt(e))
    gw_cut = 0
    for _ in range(20):
        r = rng.normal(0, 1, n)
        spins = np.sign(L @ r)
        spins[spins == 0] = 1
        cut = get_cut_value(A, spins)
        if cut > gw_cut:
            gw_cut = cut
    return float(sdp_bound), float(gw_cut)

def brute_force_max_cut(A):
    """Exact Max-Cut via Gray-code 2^(n-1) enumeration with incremental updates.

    Fix 8 — used only when exact_small picks n in [20, 30]. Spin 0 fixed to +1
    by anti-symmetry of the problem. (As)_k is maintained incrementally so each
    bit-flip costs O(n) rather than O(n^2).
    """
    n = A.shape[0]
    if n > 30:
        return float('nan')  # Don't enumerate beyond 30 — caller must set best_known
    spins = np.ones(n, dtype=np.float64)
    As = A @ spins
    sAs = float(spins @ As)
    total_A = float(A.sum())
    best_cut = 0.25 * (total_A - sAs)
    half = 1 << (n - 1)
    for i in range(1, half):
        low = i & -i
        bit = int(low).bit_length() - 1
        k = bit + 1
        sAs += -4.0 * spins[k] * As[k]
        As -= 2.0 * A[:, k] * spins[k]
        spins[k] = -spins[k]
        cut = 0.25 * (total_A - sAs)
        if cut > best_cut:
            best_cut = cut
    return float(best_cut)

def calc_lyapunov(theta, W, K_sign, lam):
    """Lyapunov function for the budgeted Kuramoto + Hebbian system.

    L = -K_sign * sum_{i<j} W_ij cos(theta_i - theta_j) + (lambda/4) * ||W||^2_F

    Fix 3 — coefficient is lam/4 (not lam/2). For symmetric, zero-diagonal W,
    np.sum(W**2) counts each undirected edge twice, so the proper regulariser
    is lam * 0.5 * (sum over edges) = lam * 0.5 * 0.5 * np.sum(W**2) = (lam/4) * sum(W**2).
    """
    theta_diff = theta[:, None] - theta[None, :]
    energy = -K_sign * np.sum(np.triu(W * np.cos(theta_diff), 1))
    reg = (lam / 4.0) * np.sum(W**2)
    return float(energy + reg)

# =============================================================================
# 3. TOP-K SYMMETRIC SUPPORT (fix 6)
# =============================================================================

def top_k_symmetric_support(W, k, rng):
    """Symmetric union top-k support with seeded tie-breaking.

    Edge (i,j) is in support iff i in j's top-k OR j in i's top-k. Tiny seeded
    noise breaks ties deterministically per graph_seed.
    """
    N = W.shape[0]
    noise = rng.uniform(0, 1e-9, W.shape)
    W_noisy = W + noise
    in_top_k = np.zeros_like(W, dtype=bool)
    for i in range(N):
        idx = np.argsort(W_noisy[i])[::-1][:k]
        in_top_k[i, idx] = True
    support = in_top_k | in_top_k.T
    np.fill_diagonal(support, False)
    return support

# =============================================================================
# 4. ODE RHS & STAGGERED INTEGRATION
# =============================================================================

def kuramoto_rhs_fast(t, theta, omega, K, W_eff):
    theta_diff = theta - theta[:, None]
    return omega + K * np.sum(W_eff * np.sin(theta_diff), axis=1)

def run_static_ode(theta0, omega, K, W_eff, T_static):
    sol = solve_ivp(kuramoto_rhs_fast, (0, T_static), theta0, args=(omega, K, W_eff),
                    method='DOP853', rtol=1e-6, atol=1e-8)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    return sol.y[:, -1], sol.nfev

def run_hebbian_staggered(theta0, W0, omega, K_sign, a_gains, prob, Y_param, X_var,
                          proj_type, T_adapt, dt, eta, lam, A_mask, budget):
    """Discrete staggered ODE-then-project loop. Fix 9: L_lower assertion per log step."""
    N = len(theta0)
    steps = int(T_adapt / dt)
    theta = theta0.copy()
    W = W0.copy()
    a_outer = np.outer(a_gains, a_gains)

    nfev_total = 0
    proj_calls = 0
    adaptive_updates = 0
    energy_violations = 0

    L_lower = -abs(K_sign) * budget * (N / 2.0)
    prev_L = calc_lyapunov(theta, W, K_sign, lam)
    final_L = prev_L

    first = True
    for step in range(steps):
        # 1. Fast ODE step on the current (frozen) effective coupling
        W_eff = W * a_outer
        sol = solve_ivp(kuramoto_rhs_fast, (0, dt), theta, args=(omega, K_sign, W_eff),
                        method='DOP853', rtol=1e-5, atol=1e-7)
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed: {sol.message}")
        theta = sol.y[:, -1]
        nfev_total += sol.nfev

        # 2. Hebbian raw update
        theta_diff = theta[:, None] - theta[None, :]
        # Sign-aware: anti-phase reinforced when K_sign < 0
        dW = eta * (K_sign * np.cos(theta_diff)) - lam * W
        W_raw = W + dt * dW

        # 3. Projection (Frobenius or Sinkhorn)
        if proj_type == 'frobenius':
            W = project_frobenius(W_raw, prob, Y_param, X_var, first_call=first)
            first = False
        else:
            W = project_sinkhorn(W_raw, budget, A_mask)
        proj_calls += 1
        adaptive_updates += 1

        # 4. Lyapunov check (fix 9: L_lower assertion; amendment 4: relative threshold)
        current_L = calc_lyapunov(theta, W, K_sign, lam)
        delta_L = current_L - prev_L
        viol_threshold = max(1e-8, 1e-6 * abs(prev_L))
        if delta_L > viol_threshold:
            energy_violations += 1
        if current_L < L_lower - 1e-6:
            print(f"WARNING: L < L_lower at step={step} L={current_L:.4f} L_lower={L_lower:.4f} delta={L_lower-current_L:.6f}")
        prev_L = current_L
        final_L = current_L

    return theta, W, nfev_total, proj_calls, adaptive_updates, energy_violations, final_L

# =============================================================================
# 5. CORE METHOD DISPATCHER
# =============================================================================

# Canonical return tuple from run_method — 20 items, see comment below.
RET_FIELDS = [
    'cut_raw', 'cut_rounded', 'cut_polished',
    'sdp_bound', 'gw_cut', 'exact_opt',
    'active_edges', 'mean_active_degree', 'max_active_degree',
    'row_sum_mean', 'row_sum_max', 'budget_violation_max', 'symmetry_error',
    'energy_descent_violations', 'final_L',
    'runtime_sec', 'projection_calls', 'ode_nfev', 'adaptive_updates',
    'W_final',
]

def _zero_metrics():
    """Placeholder structural metrics for methods that don't produce a W."""
    return dict(active_edges=0, mean_active_degree=0.0, max_active_degree=0,
                row_sum_mean=0.0, row_sum_max=0.0, budget_violation_max=0.0,
                symmetry_error=0.0)

def _structural_metrics(W_final, budget):
    """Extract diagnostic metrics from a final weight matrix."""
    active_mask = W_final > 1e-6
    active_edges = int(np.sum(active_mask) / 2)
    degrees = active_mask.sum(axis=1)
    mean_deg = float(np.mean(degrees))
    max_deg = int(np.max(degrees))
    row_sums = W_final.sum(axis=1)
    row_mean = float(np.mean(row_sums))
    row_max = float(np.max(row_sums))
    bud_max = float(np.max(row_sums - budget))
    fro_W = np.linalg.norm(W_final, 'fro') + 1e-12
    sym_err = float(np.linalg.norm(W_final - W_final.T, 'fro') / fro_W)
    return dict(active_edges=active_edges, mean_active_degree=mean_deg, max_active_degree=max_deg,
                row_sum_mean=row_mean, row_sum_max=row_max, budget_violation_max=bud_max,
                symmetry_error=sym_err)

def run_method(method, A, theta0, omega, a_gains, budget_value,
               prob, Y_param, X_var, rngs,
               learned_support_mask=None, fixed_support_mask=None):
    """One method × one seed. Returns 20-tuple per RET_FIELDS.

    rngs: dict of named np.random.Generator sub-streams (Amendment C).
    """
    N = A.shape[0]
    A_mask = (A > 0)
    K_sign = -1.0  # AFM for Max-Cut

    T_adapt = 50.0
    T_static = 50.0
    dt = 0.5
    eta = 1.0
    lam = 0.1

    start_time = time.time()

    theta_final = theta0.copy()
    W_final = A.astype(float)
    nfev = proj_calls = adaptive_updates = energy_viol = 0
    final_L = 0.0
    sdp_bound = gw_cut = exact_opt = 0.0
    cut_raw = cut_rounded = cut_polished = 0.0
    metrics = _zero_metrics()
    a_outer = np.outer(a_gains, a_gains)

    # Effective support mask: override-able for fixed_mask budget modes
    eff_mask = fixed_support_mask if fixed_support_mask is not None else A_mask

    # -------------------------------------------------------------------------
    # Static family
    # -------------------------------------------------------------------------
    if method in ('static', 'random_rounding'):
        W_final = A.astype(float)
        W_eff = W_final * a_outer
        theta_final, nfev = run_static_ode(theta0, omega, K_sign, W_eff, T_static)

    elif method == 'static_projected':
        W_final = project_frobenius(A.astype(float), prob, Y_param, X_var, first_call=True)
        proj_calls += 1
        W_eff = W_final * a_outer
        theta_final, nfev = run_static_ode(theta0, omega, K_sign, W_eff, T_static)

    elif method == 'static_budgeted':
        k = max(1, int(round(budget_value)))
        support = top_k_symmetric_support(A.astype(float), k, rngs['support'])
        W_crop = np.where(support, A.astype(float), 0.0)
        W_crop = 0.5 * (W_crop + W_crop.T)
        W_final = project_frobenius(W_crop, prob, Y_param, X_var, first_call=True)
        proj_calls += 1
        W_eff = W_final * a_outer
        theta_final, nfev = run_static_ode(theta0, omega, K_sign, W_eff, T_static)

    elif method == 'random_budgeted':
        W_rand = rngs['random_init'].uniform(0, 1, (N, N))
        W_rand = 0.5 * (W_rand + W_rand.T)
        W_rand[~eff_mask] = 0
        W_final = project_frobenius(W_rand, prob, Y_param, X_var, first_call=True)
        proj_calls += 1
        W_eff = W_final * a_outer
        theta_final, nfev = run_static_ode(theta0, omega, K_sign, W_eff, T_static)

    elif method == 'learned_support_random_weights':
        # Fix 11: use cached support from hebbian_frobenius m_seed=0
        if learned_support_mask is None:
            raise RuntimeError("learned_support_random_weights called without learned_support_mask")
        W_rand = rngs['random_init'].uniform(0, 1, (N, N))
        W_rand = 0.5 * (W_rand + W_rand.T)
        W_rand[~learned_support_mask] = 0
        W_final = project_frobenius(W_rand, prob, Y_param, X_var, first_call=True)
        proj_calls += 1
        W_eff = W_final * a_outer
        theta_final, nfev = run_static_ode(theta0, omega, K_sign, W_eff, T_static)

    elif method == 'greedy':
        # Greedy starts from cos(theta0) sign, polishes — no Kuramoto involved
        spins_raw = np.sign(np.cos(theta0))
        spins_raw[spins_raw == 0] = 1
        polished = greedy_max_cut(A, spins_raw)
        c = get_cut_value(A, polished)
        cut_raw = cut_rounded = cut_polished = float(c)
        runtime = time.time() - start_time
        return (cut_raw, cut_rounded, cut_polished, 0.0, 0.0, 0.0,
                metrics['active_edges'], metrics['mean_active_degree'], metrics['max_active_degree'],
                metrics['row_sum_mean'], metrics['row_sum_max'], metrics['budget_violation_max'], metrics['symmetry_error'],
                0, 0.0, runtime, 0, 0, 0, None)

    elif method == 'sdp_gw':
        sdp_bound, gw_cut = run_sdp_gw(A, rngs['rounding'])
        runtime = time.time() - start_time
        # Echo gw_cut into cut_polished so SDP results appear in standard
        # cross-method aggregations; sdp_bound is the relaxation upper bound.
        return (gw_cut, gw_cut, gw_cut, sdp_bound, gw_cut, 0.0,
                metrics['active_edges'], metrics['mean_active_degree'], metrics['max_active_degree'],
                metrics['row_sum_mean'], metrics['row_sum_max'], metrics['budget_violation_max'], metrics['symmetry_error'],
                0, 0.0, runtime, 0, 0, 0, None)

    # -------------------------------------------------------------------------
    # Hebbian family
    # -------------------------------------------------------------------------
    elif method in ('hebbian_frobenius', 'hebbian_sinkhorn',
                    'topology_scrambled'):
        proj_type = 'frobenius' if 'frobenius' in method else (
            'sinkhorn' if 'sinkhorn' in method else 'frobenius')
        mask = A_mask.copy()
        local_prob, local_Yp, local_Xv = prob, Y_param, X_var

        if method == 'topology_scrambled':
            # Fix 12: nswap = 10 * |E|
            Gp = nx.from_numpy_array(A)
            nswap = 10 * len(Gp.edges)
            try:
                nx.double_edge_swap(Gp, nswap=nswap, max_tries=10 * nswap)
            except nx.NetworkXAlgorithmError:
                pass  # graph too small to swap; fall through with original mask
            mask = nx.to_numpy_array(Gp) > 0
            # Recompile projector for scrambled mask
            local_prob, local_Yp, local_Xv = make_frobenius_projector(N, mask, budget_value)
            proj_type = 'frobenius'

        # Init W0 by projecting A onto the (possibly scrambled) mask
        if proj_type == 'frobenius':
            W0 = project_frobenius(A.astype(float) * mask.astype(float), local_prob, local_Yp, local_Xv, first_call=True)
            proj_calls += 1
        else:
            W0 = project_sinkhorn(A.astype(float) * mask.astype(float), budget_value, mask)

        theta_final, W_final, nfev, pc, au, energy_viol, final_L = run_hebbian_staggered(
            theta0, W0, omega, K_sign, a_gains, local_prob, local_Yp, local_Xv,
            proj_type, T_adapt, dt, eta, lam, mask, budget_value)
        proj_calls += pc
        adaptive_updates += au

    elif method in ('hybrid_frobenius', 'hybrid_sinkhorn'):
        proj_type = 'frobenius' if 'frobenius' in method else 'sinkhorn'
        if proj_type == 'frobenius':
            W0 = project_frobenius(A.astype(float), prob, Y_param, X_var, first_call=True)
            proj_calls += 1
        else:
            W0 = project_sinkhorn(A.astype(float), budget_value, A_mask)

        # Phase 1: adaptive
        theta_mid, W_frozen, nfev1, pc, au, energy_viol, final_L = run_hebbian_staggered(
            theta0, W0, omega, K_sign, a_gains, prob, Y_param, X_var,
            proj_type, T_adapt, dt, eta, lam, A_mask, budget_value)
        # Phase 2: frozen static
        W_eff = W_frozen * a_outer
        theta_final, nfev2 = run_static_ode(theta_mid, omega, K_sign, W_eff, T_static)
        W_final = W_frozen
        nfev = nfev1 + nfev2
        proj_calls += pc
        adaptive_updates += au

    else:
        raise ValueError(f"Unknown method: {method}")

    # Hyperplane rounding + greedy polish
    cut_raw, cut_rounded, cut_polished = random_hyperplane_rounding(
        theta_final, A, seeds=20, rng=rngs['rounding'])

    metrics = _structural_metrics(W_final, budget_value)
    runtime = time.time() - start_time

    return (cut_raw, cut_rounded, cut_polished, sdp_bound, gw_cut, exact_opt,
            metrics['active_edges'], metrics['mean_active_degree'], metrics['max_active_degree'],
            metrics['row_sum_mean'], metrics['row_sum_max'], metrics['budget_violation_max'], metrics['symmetry_error'],
            energy_viol, final_L,
            runtime, proj_calls, nfev, adaptive_updates,
            W_final)

# =============================================================================
# 6. RNG MANAGEMENT (Amendment C)
# =============================================================================

def make_rngs(graph_seed, method_seed):
    """Spawn deterministic sub-streams from a single root rng.

    Sub-streams returned: 'phases', 'omega', 'amplitude', 'rounding', 'random_init',
    'support', 'aux'. Each is a fresh np.random.Generator derived via spawn().
    """
    root = np.random.default_rng((int(graph_seed) + 1) * 1_000_003 + int(method_seed))
    streams = root.spawn(7)
    return {
        'phases': streams[0],
        'omega': streams[1],
        'amplitude': streams[2],
        'rounding': streams[3],
        'random_init': streams[4],
        'support': streams[5],
        'aux': streams[6],
    }

# =============================================================================
# 7. BENCHMARK ORCHESTRATION
# =============================================================================

ALL_METHODS = [
    'static', 'random_rounding', 'greedy',
    'static_projected', 'static_budgeted', 'random_budgeted',
    'learned_support_random_weights', 'topology_scrambled',
    'hebbian_sinkhorn', 'hybrid_sinkhorn',
    'hebbian_frobenius', 'hybrid_frobenius',
    'sdp_gw',
]

def build_graph(family, graph_seed, rng, fixed_n=None):
    """Returns (G, n, m, mean_degree, best_known_or_nan).

    fixed_n: optional override for families with random n (exact_small).
    """
    if family == 'sparse_er_200_p05':
        G = nx.erdos_renyi_graph(200, 0.05, seed=graph_seed)
        best = float('nan')
    elif family == 'dense_er_200_p15':
        G = nx.erdos_renyi_graph(200, 0.15, seed=graph_seed)
        best = float('nan')
    elif family == 'random_regular_200_d10':
        G = nx.random_regular_graph(10, 200, seed=graph_seed)
        best = float('nan')
    elif family == 'sbm_200':
        G = nx.planted_partition_graph(2, 100, 0.3, 0.05, seed=graph_seed)
        best = float('nan')
    elif family == 'exact_small':
        # Fix 8: n in [20, 30]; pilot pins n=20 via fixed_n
        n_target = int(fixed_n) if fixed_n is not None else int(rng.integers(20, 31))
        G = nx.erdos_renyi_graph(n_target, 0.5, seed=graph_seed)
        best = float('nan')  # filled by brute_force later
    elif family == 'gset_g1':
        # Fix 7: real G1, skip if missing
        if not os.path.exists('G1'):
            return None, 0, 0, 0.0, 11624.0
        G = load_gset('G1')
        best = 11624.0
    else:
        raise ValueError(family)
    n = G.number_of_nodes()
    m = G.number_of_edges()
    mean_deg = 2.0 * m / n if n > 0 else 0.0
    return G, n, m, mean_deg, best

def save_W_snapshot(W_final, out_dir, graph_family, graph_seed, method, method_seed):
    """Amendment B: persist final W for hebbian-class methods."""
    if W_final is None:
        return
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{graph_family}_{graph_seed}_{method}_{method_seed}.npz")
    np.savez_compressed(path, W=W_final)

def compute_sdp_scope_ok(family, graph_seed, n):
    """Fix 15: SDP only on exact_small (always) and sparse_er_200_p05 g_seed=0."""
    if family == 'exact_small':
        return True
    if family == 'sparse_er_200_p05' and graph_seed == 0:
        return True
    return False

def precompute_learned_support(G, A, mean_deg, budget_value, A_mask, n):
    """Fix 11: run hebbian_frobenius m_seed=0 once to obtain canonical support mask."""
    rngs = make_rngs(graph_seed=0, method_seed=0)
    theta0 = rngs['phases'].uniform(0, 2 * np.pi, n)
    omega = np.zeros(n)
    a_gains = np.ones(n)
    prob, Y_param, X_var = make_frobenius_projector(n, A_mask, budget_value)
    res = run_method('hebbian_frobenius', A, theta0, omega, a_gains, budget_value,
                     prob, Y_param, X_var, rngs)
    W_final = res[-1]
    mask = W_final > 1e-6
    return mask, W_final, prob, Y_param, X_var

# =============================================================================
# 8. PILOT (Amendment A)
# =============================================================================

# Pilot scope: 3 ER-class families plus one exact-small case (amendment 3)
PILOT_FAMILIES = ['sparse_er_200_p05', 'dense_er_200_p15', 'random_regular_200_d10', 'exact_small']
PILOT_GRAPH_SEEDS = [0]
PILOT_METHOD_SEEDS = list(range(3))
PILOT_EXACT_SMALL_N = 20  # pinned by amendment 3
W_SNAPSHOT_DIR = 'results/W_snapshots'
HEBBIAN_LIKE = ('hebbian_frobenius', 'hebbian_sinkhorn', 'hybrid_frobenius', 'hybrid_sinkhorn',
                'topology_scrambled')

def run_pilot():
    """Restricted run. Amendment A pass gates evaluated after."""
    print("\n--- PILOT MODE (Amendment A) ---")
    print(f"Families: {PILOT_FAMILIES}")
    print(f"Graph seeds: {PILOT_GRAPH_SEEDS}, Method seeds: {PILOT_METHOD_SEEDS}")

    os.makedirs('results', exist_ok=True)
    os.makedirs(W_SNAPSHOT_DIR, exist_ok=True)

    print("\n[unit test]")
    unit_test_projection()
    print("[unit test] OK\n")

    print("[profile] dual solver at n=200")
    profile_dual_at_n200(n_calls=50, target_ms=50.0)
    print("[profile] OK\n")

    rows = []
    projector_times = []
    energy_violation_total = 0

    total_runtime_accumulator = 0.0
    _proj_reset()
    proj_time_start_pilot, proj_count_start_pilot = _proj_snapshot()

    for family in PILOT_FAMILIES:
        for g_seed in PILOT_GRAPH_SEEDS:
            graph_rng = np.random.default_rng(g_seed)
            fixed_n = PILOT_EXACT_SMALL_N if family == 'exact_small' else None
            G, n, m, mean_deg, best_known = build_graph(family, g_seed, graph_rng, fixed_n=fixed_n)
            if G is None:
                print(f"Skipping {family} g_seed={g_seed} (graph unavailable)")
                continue
            A = nx.to_numpy_array(G)
            np.fill_diagonal(A, 0)
            A_mask = (A > 0)
            budget_value = mean_deg / 2.0
            budget_mode = 'half_mean_deg'

            # Exact-small: compute brute-force optimum and tag
            if family == 'exact_small':
                bf_t0 = time.time()
                best_known = brute_force_max_cut(A)
                exact_status_for_family = 'proven_optimal'
                print(f"  [brute force] n={n} optimum = {best_known:.0f} ({time.time()-bf_t0:.1f}s)")
            elif family == 'gset_g1':
                exact_status_for_family = 'best_known'
            else:
                exact_status_for_family = 'timeout'

            print(f"\n=== {family} (n={n}, m={m}, mean_deg={mean_deg:.1f}, budget={budget_value:.2f}) ===")

            # Precompute learned support
            print("  [precompute] hebbian_frobenius m_seed=0 for canonical support...")
            t0 = time.time()
            learned_mask, _W_learned, prob, Y_param, X_var = precompute_learned_support(
                G, A, mean_deg, budget_value, A_mask, n)
            print(f"  [precompute] support density: {learned_mask.sum() / (n * (n - 1)):.4f}  ({time.time()-t0:.1f}s)")

            for method in ALL_METHODS:
                # Fix 15: SDP scope guard
                if method == 'sdp_gw' and not compute_sdp_scope_ok(family, g_seed, n):
                    continue
                for m_seed in PILOT_METHOD_SEEDS:
                    if method == 'sdp_gw' and m_seed > 0:
                        continue
                    rngs = make_rngs(g_seed, m_seed)
                    theta0 = rngs['phases'].uniform(0, 2 * np.pi, n)
                    omega = np.zeros(n)
                    a_gains = np.ones(n)

                    # Per-call CVXPY timing window
                    pt_before, pc_before = _proj_snapshot()
                    res = run_method(method, A, theta0, omega, a_gains, budget_value,
                                     prob, Y_param, X_var, rngs,
                                     learned_support_mask=learned_mask)
                    pt_after, pc_after = _proj_snapshot()
                    cvxpy_time_this_run = pt_after - pt_before
                    cvxpy_calls_this_run = pc_after - pc_before

                    runtime = res[15]
                    proj_calls = res[16]
                    W_final = res[-1]
                    total_runtime_accumulator += runtime

                    if cvxpy_calls_this_run > 0:
                        projector_times.append((family, method, m_seed, n,
                                                cvxpy_time_this_run, cvxpy_calls_this_run,
                                                runtime))

                    energy_violation_total += res[13]

                    # Amendment B: snapshot for hebbian-like methods
                    if method in HEBBIAN_LIKE and W_final is not None:
                        save_W_snapshot(W_final, W_SNAPSHOT_DIR, family, g_seed, method, m_seed)

                    L_valid = (method in ('hebbian_frobenius', 'hybrid_frobenius')
                               and 'frobenius' in method and budget_mode == 'half_mean_deg')

                    row = {
                        'graph_family': family, 'graph_seed': g_seed, 'n': n, 'm': m,
                        'mean_degree': mean_deg,
                        'budget_mode': budget_mode, 'budget_value': budget_value,
                        'projection_type': ('frobenius' if 'frobenius' in method else
                                            ('sinkhorn' if 'sinkhorn' in method else 'none')),
                        'method': method, 'amplitude_sigma': 0.0, 'amplitude_label': 'ideal',
                        'omega_mode': 'zero', 'restart_seed': m_seed,
                        'support_source': 'hebbian_seed0' if method == 'learned_support_random_weights' else 'graph',
                        'T_adapt': 50.0, 'T_static': 50.0,
                        'adaptive_updates': res[18], 'projection_calls': res[16], 'ode_nfev': res[17],
                        'runtime_sec': res[15],
                        'cvxpy_time_sec': cvxpy_time_this_run,
                        'cvxpy_calls': cvxpy_calls_this_run,
                        'cut_raw': res[0], 'cut_rounded': res[1], 'cut_polished': res[2],
                        'sdp_bound': res[3], 'gw_cut': res[4],
                        'exact_optimum': best_known if not np.isnan(best_known) else None,
                        'exact_status': exact_status_for_family,
                        'cut_ratio_exact': (res[2] / best_known) if best_known and not np.isnan(best_known) and best_known > 0 else 0.0,
                        'active_edges': res[6], 'mean_active_degree': res[7], 'max_active_degree': res[8],
                        'row_sum_mean': res[9], 'row_sum_max': res[10], 'budget_violation_max': res[11],
                        'symmetry_error': res[12],
                        'energy_descent_violations': res[13] if L_valid else 0,
                        'energy_violations_raw': res[13],
                        'final_L': res[14],
                        'L_valid': L_valid,
                    }
                    rows.append(row)
                    print(f"    {method:32s} m_seed={m_seed} cut_pol={res[2]:.0f} t={runtime:.2f}s pc={proj_calls} viol={res[13]}")

    df = pd.DataFrame(rows)
    df.to_csv('results/pilot_results.csv', index=False)
    print(f"\nWrote results/pilot_results.csv  ({len(rows)} rows)")

    # ---------------- Amendment 1: best-of-restarts summary per (family, method) ----------------
    print("\n--- BEST-OF-RESTARTS SUMMARY (statistical unit = graph instance) ---")
    bor = (df.groupby(['graph_family', 'method'])
             .agg(best_cut=('cut_polished', 'max'),
                  mean_cut=('cut_polished', 'mean'),
                  cut_ratio=('cut_ratio_exact', 'max'),
                  exact=('exact_optimum', 'first'),
                  status=('exact_status', 'first'))
             .reset_index())
    print(bor.to_string(index=False))

    # ---------------- Amendment A pass-gate evaluation ----------------
    print("\n--- PILOT PASS-GATE EVALUATION ---")
    # 1. Unit tests passed (we got here)
    gate_unit = True
    print(f"[gate 1] unit tests: PASS")

    # 2. symmetry_error < 1e-6 on every row
    sym_max = float(df['symmetry_error'].max())
    gate_sym = sym_max < 1e-6
    print(f"[gate 2] symmetry_error max {sym_max:.3e}: {'PASS' if gate_sym else 'FAIL'}")

    # 3. budget_violation_max < 1e-6 on rows that actually project (filter out
    # 'static', 'random_rounding', 'greedy' whose W_final = A by design).
    projecting = df[df['projection_calls'] > 0]
    bud_max = float(projecting['budget_violation_max'].max()) if len(projecting) else 0.0
    gate_bud = bud_max < 1e-6
    print(f"[gate 3] budget_violation_max on projecting methods {bud_max:.3e}: "
          f"{'PASS' if gate_bud else 'FAIL'}")
    if not gate_bud:
        worst_bud = projecting.nlargest(5, 'budget_violation_max')[
            ['graph_family', 'method', 'restart_seed', 'budget_violation_max']]
        print(worst_bud.to_string(index=False))

    # 4. L_valid rows: zero meaningful energy_descent_violations
    hebbian_rows = df[df['L_valid'] == True]
    total_viol = int(hebbian_rows['energy_descent_violations'].sum())
    gate_L = (total_viol == 0)
    print(f"[gate 4] energy_descent_violations on L_valid rows ({len(hebbian_rows)}): {total_viol} -> {'PASS' if gate_L else 'FAIL'}")

    # 5. Cut sanity: 0 <= cut_polished <= m AND cut_polished >= cut_raw (rounding never worsens)
    bad_cuts = df[(df['cut_polished'] < 0) | (df['cut_polished'] > df['m'])]
    worse_round = df[df['cut_polished'] < df['cut_raw'] - 1e-6]
    gate_cuts = len(bad_cuts) == 0 and len(worse_round) == 0
    print(f"[gate 5] cut sanity ({len(bad_cuts)} out-of-range, {len(worse_round)} where rounding worsens): "
          f"{'PASS' if gate_cuts else 'FAIL'}")
    if len(worse_round) > 0:
        print(worse_round[['graph_family', 'method', 'restart_seed', 'cut_raw', 'cut_polished']].to_string(index=False))

    # 6. CVXPY projector wall time < 500ms per call at n=200
    proj_n200 = df[(df['n'] == 200) & (df['cvxpy_calls'] > 0)]
    if len(proj_n200) > 0:
        per_call_ms = (proj_n200['cvxpy_time_sec'] / proj_n200['cvxpy_calls']) * 1000.0
        mean_proj_ms = float(per_call_ms.mean())
        max_proj_ms = float(per_call_ms.max())
        gate_proj_speed = mean_proj_ms < 500.0
        print(f"[gate 6] projector wall mean {mean_proj_ms:.1f}ms (max {max_proj_ms:.1f}ms) at n=200: "
              f"{'PASS' if gate_proj_speed else 'FAIL'}")
    else:
        mean_proj_ms = max_proj_ms = 0.0
        gate_proj_speed = False
        print(f"[gate 6] no n=200 projector runs - FAIL")

    # 7. Projection time fraction (informational only now — dual is the runtime path).
    # Original CVXPY-80%-trigger has already been remediated by switching to dual.
    total_proj = float(df['cvxpy_time_sec'].sum())
    total_runtime = float(df['runtime_sec'].sum()) or 1.0
    proj_frac = total_proj / total_runtime
    print(f"[gate 7] dual-projection time fraction = {proj_frac:.2%}: INFO "
          f"(CVXPY-80%-remediation already done; dual is the runtime path)")

    overall = (gate_unit and gate_sym and gate_bud and gate_L and gate_cuts
               and gate_proj_speed)
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}  "
          f"(mean projector {mean_proj_ms:.1f}ms, projection frac {proj_frac:.2%})")
    return overall, mean_proj_ms

# =============================================================================
# 9. FULL PHASE 2 (not run by default)
# =============================================================================

FULL_FAMILIES = ['sparse_er_200_p05', 'dense_er_200_p15', 'random_regular_200_d10', 'exact_small']
FULL_ER_GRAPH_SEEDS = [0, 1, 2]
FULL_GSET_INCLUDED = True
FULL_METHOD_SEEDS = list(range(10))
FULL_BUDGET_MODES = [('half_mean_deg', 0.5), ('mean_deg', 1.0)]  # multiplier on mean_degree

def _run_cell(family, G, g_seed, n, m, mean_deg, best_known, exact_status,
              budget_value, budget_mode, omega_mode, method_seeds,
              learned_mask, prob, Y_param, X_var,
              save_snapshots=True):
    """Run all methods × method_seeds on a single (family, g_seed, budget_mode, omega_mode) cell.
    Returns a list of row dicts (one per method × method_seed)."""
    A = nx.to_numpy_array(G)
    np.fill_diagonal(A, 0)

    rows = []
    for method in ALL_METHODS:
        if method == 'sdp_gw' and not compute_sdp_scope_ok(family, g_seed, n):
            continue
        for m_seed in method_seeds:
            if method == 'sdp_gw' and m_seed > 0:
                continue
            rngs = make_rngs(g_seed, m_seed)
            theta0 = rngs['phases'].uniform(0, 2 * np.pi, n)
            if omega_mode == 'zero':
                omega = np.zeros(n)
            else:
                omega = rngs['omega'].normal(0, 0.3, n)
            a_gains = np.ones(n)

            pt_before, pc_before = _proj_snapshot()
            res = run_method(method, A, theta0, omega, a_gains, budget_value,
                             prob, Y_param, X_var, rngs,
                             learned_support_mask=learned_mask)
            pt_after, pc_after = _proj_snapshot()
            cvxpy_time_this_run = pt_after - pt_before
            cvxpy_calls_this_run = pc_after - pc_before

            W_final = res[-1]
            if save_snapshots and method in HEBBIAN_LIKE and W_final is not None:
                tag = f"{family}_{g_seed}_{budget_mode}_{omega_mode}"
                save_W_snapshot(W_final, W_SNAPSHOT_DIR, tag, '', method, m_seed)

            L_valid = (method in ('hebbian_frobenius', 'hybrid_frobenius')
                       and omega_mode == 'zero'
                       and budget_mode == 'half_mean_deg')

            row = {
                'graph_family': family, 'graph_seed': g_seed, 'n': n, 'm': m,
                'mean_degree': mean_deg,
                'budget_mode': budget_mode, 'budget_value': budget_value,
                'projection_type': ('frobenius' if 'frobenius' in method else
                                    ('sinkhorn' if 'sinkhorn' in method else 'none')),
                'method': method, 'amplitude_sigma': 0.0, 'amplitude_label': 'ideal',
                'omega_mode': omega_mode, 'restart_seed': m_seed,
                'support_source': 'hebbian_seed0' if method == 'learned_support_random_weights' else 'graph',
                'T_adapt': 50.0, 'T_static': 50.0,
                'adaptive_updates': res[18], 'projection_calls': res[16], 'ode_nfev': res[17],
                'runtime_sec': res[15],
                'cvxpy_time_sec': cvxpy_time_this_run,
                'cvxpy_calls': cvxpy_calls_this_run,
                'cut_raw': res[0], 'cut_rounded': res[1], 'cut_polished': res[2],
                'sdp_bound': res[3], 'gw_cut': res[4],
                'exact_optimum': best_known if not np.isnan(best_known) else None,
                'exact_status': exact_status,
                'cut_ratio_exact': (res[2] / best_known) if best_known and not np.isnan(best_known) and best_known > 0 else 0.0,
                'active_edges': res[6], 'mean_active_degree': res[7], 'max_active_degree': res[8],
                'row_sum_mean': res[9], 'row_sum_max': res[10], 'budget_violation_max': res[11],
                'symmetry_error': res[12],
                'energy_descent_violations': res[13] if L_valid else 0,
                'energy_violations_raw': res[13],
                'final_L': res[14],
                'L_valid': L_valid,
            }
            rows.append(row)
            print(f"    {method:32s} m_seed={m_seed} cut_pol={res[2]:.0f} "
                  f"t={res[15]:.1f}s pc={res[16]} viol={res[13]}")
    return rows

def run_heterogeneity_experiment():
    """Amplitude heterogeneity ablation on sparse_er_200_p05."""
    print("\n--- HETEROGENEITY EXPERIMENT (sparse_er, lognormal gains) ---")
    sigmas = [0.0, 0.25, 0.5, 1.0]
    het_methods = ['static_budgeted', 'hebbian_frobenius', 'hybrid_frobenius', 'oracle_compensation']
    rows = []
    n = 200
    for g_seed in range(3):
        graph_rng = np.random.default_rng(g_seed)
        G = nx.erdos_renyi_graph(n, 0.05, seed=g_seed)
        A = nx.to_numpy_array(G)
        np.fill_diagonal(A, 0)
        A_mask = (A > 0)
        m = G.number_of_edges()
        mean_deg = 2.0 * m / n
        budget_value = mean_deg / 2.0
        prob, Y_param, X_var = make_frobenius_projector(n, A_mask, budget_value)

        # Precompute learned mask
        learned_mask, _, _, _, _ = precompute_learned_support(G, A, mean_deg, budget_value, A_mask, n)

        for sigma in sigmas:
            for m_seed in range(10):  # match main suite seed count
                rngs = make_rngs(g_seed, m_seed)
                theta0 = rngs['phases'].uniform(0, 2 * np.pi, n)
                # Lognormal with unit-mean correction (fix 4)
                a_gains = rngs['amplitude'].lognormal(-0.5 * sigma**2, sigma, n)
                omega = np.zeros(n)

                for method in het_methods:
                    t0 = time.time()
                    if method == 'oracle_compensation':
                        # Oracle: W chosen to cancel heterogeneity at unit static budget
                        a_outer = np.outer(a_gains, a_gains)
                        W_oracle = np.where(A_mask, A.astype(float) / (a_outer + 1e-12), 0.0)
                        W_oracle = 0.5 * (W_oracle + W_oracle.T)
                        W_final = project_frobenius_dual(W_oracle, A_mask, budget_value)
                        W_eff = W_final * a_outer
                        sol = solve_ivp(kuramoto_rhs_fast, (0, 200), theta0,
                                        args=(omega, -0.5, W_eff),
                                        method='DOP853', rtol=1e-6, atol=1e-8)
                        if not sol.success:
                            raise RuntimeError(f"solve_ivp failed: {sol.message}")
                        cut_raw, cut_rounded, cut_polished = random_hyperplane_rounding(
                            sol.y[:, -1], A, seeds=20, rng=rngs['rounding'])
                        corr_W_inv_gain = 1.0
                        corr_W_gain = -1.0
                    else:
                        res = run_method(method, A, theta0, omega, a_gains, budget_value,
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
                        'cut_raw': cut_raw, 'cut_rounded': cut_rounded, 'cut_polished': cut_polished,
                        'corr_W_inverse_gain': corr_W_inv_gain,
                        'corr_W_gain': corr_W_gain,
                        'runtime_sec': time.time() - t0,
                    })
                    print(f"    g={g_seed} m={m_seed} sigma={sigma} {method:25s} "
                          f"cut={cut_polished:.0f} t={time.time()-t0:.1f}s")
    df = pd.DataFrame(rows)
    # Compute cut_ratio_vs_baseline (vs static_budgeted at sigma=0 for each graph_seed)
    baseline = (df[(df['method'] == 'static_budgeted') & (df['sigma'] == 0.0)]
                  .groupby('graph_seed')['cut_polished'].mean().to_dict())
    df['cut_ratio_vs_baseline'] = df.apply(
        lambda r: r['cut_polished'] / baseline.get(r['graph_seed'], 1.0)
                  if baseline.get(r['graph_seed'], 0.0) > 0 else 0.0,
        axis=1)
    df.to_csv('results/heterogeneity_experiment.csv', index=False)
    print(f"\nWrote results/heterogeneity_experiment.csv ({len(rows)} rows)")
    return df

def run_full():
    """Full Phase 2 sweep — graph_seed=0..2 for ER, 2 budget modes, all methods,
    detuning ablation on sparse_er, heterogeneity experiment. Long (~12-16 h)."""
    print("\n--- FULL PHASE 2 SUITE ---")
    print(f"Families: {FULL_FAMILIES} (+ gset_g1 if present)")
    print(f"ER graph seeds: {FULL_ER_GRAPH_SEEDS}, method seeds: {FULL_METHOD_SEEDS}")
    print(f"Budget modes: {[m for m, _ in FULL_BUDGET_MODES]}")

    os.makedirs('results', exist_ok=True)
    os.makedirs(W_SNAPSHOT_DIR, exist_ok=True)

    print("\n[unit test]")
    unit_test_projection()
    print("[unit test] OK\n")

    print("[profile] dual solver at n=200")
    profile_dual_at_n200(n_calls=20, target_ms=50.0)
    print("[profile] OK\n")

    all_rows = []
    csv_path = 'results/phase2_full_results.csv'

    # === Core loop: families × graph_seeds × budget_modes × omega_mode=zero ===
    for family in FULL_FAMILIES:
        graph_seeds = FULL_ER_GRAPH_SEEDS if family != 'exact_small' else [0]
        for g_seed in graph_seeds:
            graph_rng = np.random.default_rng(g_seed)
            fixed_n = (20 + g_seed) if family == 'exact_small' else None  # vary n in [20, 22]
            G, n, m, mean_deg, best_known = build_graph(family, g_seed, graph_rng, fixed_n=fixed_n)
            if G is None:
                print(f"Skipping {family} g_seed={g_seed}")
                continue
            A = nx.to_numpy_array(G)
            np.fill_diagonal(A, 0)
            A_mask = (A > 0)

            if family == 'exact_small':
                bf_t0 = time.time()
                best_known = brute_force_max_cut(A)
                exact_status = 'proven_optimal'
                print(f"[brute force] {family} g_seed={g_seed} n={n} opt={best_known:.0f} ({time.time()-bf_t0:.1f}s)")
            elif family == 'gset_g1':
                exact_status = 'best_known'
            else:
                exact_status = 'timeout'

            for budget_mode, mult in FULL_BUDGET_MODES:
                budget_value = mult * mean_deg
                print(f"\n=== {family} g_seed={g_seed} budget={budget_mode}({budget_value:.2f}) "
                      f"n={n} m={m} ===")

                # Precompute canonical support per (g_seed, budget_mode)
                print(f"  [precompute] hebbian_frobenius m_seed=0 for canonical support...")
                t0 = time.time()
                learned_mask, _W_learned, prob, Y_param, X_var = precompute_learned_support(
                    G, A, mean_deg, budget_value, A_mask, n)
                print(f"  [precompute] support density "
                      f"{learned_mask.sum() / (n * (n - 1)):.4f}  ({time.time()-t0:.1f}s)")

                cell_rows = _run_cell(family, G, g_seed, n, m, mean_deg, best_known, exact_status,
                                      budget_value, budget_mode, 'zero', FULL_METHOD_SEEDS,
                                      learned_mask, prob, Y_param, X_var, save_snapshots=True)
                all_rows.extend(cell_rows)
                # Incremental save so we don't lose progress on crash
                pd.DataFrame(all_rows).to_csv(csv_path, index=False)

    # === Add GSet G1 if present ===
    if FULL_GSET_INCLUDED and os.path.exists('G1'):
        family = 'gset_g1'
        g_seed = 0
        G, n, m, mean_deg, best_known = build_graph(family, g_seed, np.random.default_rng(g_seed))
        if G is not None:
            A = nx.to_numpy_array(G)
            np.fill_diagonal(A, 0)
            A_mask = (A > 0)
            exact_status = 'best_known'
            for budget_mode, mult in FULL_BUDGET_MODES:
                budget_value = mult * mean_deg
                print(f"\n=== {family} g_seed=0 budget={budget_mode}({budget_value:.2f}) n={n} m={m} ===")
                print(f"  [precompute] hebbian_frobenius m_seed=0 for canonical support...")
                t0 = time.time()
                learned_mask, _W, prob, Y_param, X_var = precompute_learned_support(
                    G, A, mean_deg, budget_value, A_mask, n)
                print(f"  [precompute] support density "
                      f"{learned_mask.sum() / (n * (n - 1)):.4f}  ({time.time()-t0:.1f}s)")
                cell_rows = _run_cell(family, G, g_seed, n, m, mean_deg, best_known, exact_status,
                                      budget_value, budget_mode, 'zero', FULL_METHOD_SEEDS,
                                      learned_mask, prob, Y_param, X_var, save_snapshots=True)
                all_rows.extend(cell_rows)
                pd.DataFrame(all_rows).to_csv(csv_path, index=False)

    # === Detuning ablation: sparse_er with omega_sigma=0.3 ===
    print("\n--- DETUNING ABLATION (sparse_er, omega_sigma=0.3) ---")
    for g_seed in FULL_ER_GRAPH_SEEDS:
        graph_rng = np.random.default_rng(g_seed)
        G, n, m, mean_deg, best_known = build_graph('sparse_er_200_p05', g_seed, graph_rng)
        A = nx.to_numpy_array(G)
        np.fill_diagonal(A, 0)
        A_mask = (A > 0)
        for budget_mode, mult in FULL_BUDGET_MODES:
            budget_value = mult * mean_deg
            print(f"\n=== sparse_er_200_p05_ablation g_seed={g_seed} "
                  f"budget={budget_mode}({budget_value:.2f}) ===")
            t0 = time.time()
            learned_mask, _W, prob, Y_param, X_var = precompute_learned_support(
                G, A, mean_deg, budget_value, A_mask, n)
            print(f"  [precompute] support density "
                  f"{learned_mask.sum() / (n * (n - 1)):.4f}  ({time.time()-t0:.1f}s)")
            cell_rows = _run_cell('sparse_er_200_p05_ablation', G, g_seed, n, m, mean_deg,
                                  best_known, 'timeout',
                                  budget_value, budget_mode, 'detuning_ablation', FULL_METHOD_SEEDS,
                                  learned_mask, prob, Y_param, X_var, save_snapshots=False)
            all_rows.extend(cell_rows)
            pd.DataFrame(all_rows).to_csv(csv_path, index=False)

    # === Final core CSV ===
    df = pd.DataFrame(all_rows)
    df.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path} ({len(df)} rows)")

    # === Heterogeneity experiment ===
    run_heterogeneity_experiment()

    print("\n--- FULL PHASE 2 COMPLETE ---")

# =============================================================================
# 10. ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else 'pilot'
    if mode == 'pilot':
        ok, mean_ms = run_pilot()
        sys.exit(0 if ok else 1)
    elif mode == 'full':
        run_full()
    else:
        print(f"Usage: {sys.argv[0]} [pilot|full]")
        sys.exit(2)
