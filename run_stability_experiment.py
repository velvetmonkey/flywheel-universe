"""Stability experiment for budgeted Hebbian-Kuramoto Max-Cut.

The headline claim is degradation under amplitude heterogeneity, not absolute
Max-Cut performance. The script writes incremental CSVs and can resume safely.

Examples:
  python run_stability_experiment.py smoke
  python run_stability_experiment.py full --max-hours 10
  python run_stability_experiment.py analyze
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from phase2_benchmark import (
    build_graph,
    enforce_symmetric_support_budget,
    get_cut_value,
    greedy_max_cut,
    kuramoto_rhs_fast,
    load_gset,
    make_frobenius_projector,
    make_rngs,
    random_hyperplane_rounding,
    run_hebbian_staggered,
    run_method,
    run_sdp_gw,
)


RESULTS_DIR = Path("results")
HET_CSV = RESULTS_DIR / "stability_heterogeneity.csv"
BENCH_CSV = RESULTS_DIR / "stability_supporting_benchmarks.csv"
SLOPE_CSV = RESULTS_DIR / "stability_degradation_slopes.csv"
LOG_PATH = RESULTS_DIR / "stability_run.log"
WRITEUP_PATH = RESULTS_DIR / "stability_writeup.md"

SIGMAS = [0.0, 0.25, 0.5, 1.0]
K_SIGN = -1.0
G1_BEST_KNOWN = 11624.0


def log(msg: str) -> None:
    print(msg, flush=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as f:
        f.write(msg + "\n")


def append_row(path: Path, row: dict) -> None:
    path.parent.mkdir(exist_ok=True)
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def load_done(path: Path, keys: list[str]) -> set[tuple]:
    if not path.exists():
        return set()
    df = pd.read_csv(path)
    if not set(keys).issubset(df.columns):
        return set()
    return set(map(tuple, df[keys].itertuples(index=False, name=None)))


def graph_to_A(G: nx.Graph) -> np.ndarray:
    A = nx.to_numpy_array(G)
    np.fill_diagonal(A, 0.0)
    return A


def final_dynamics_metrics(theta: np.ndarray, omega: np.ndarray, W_final: np.ndarray | None,
                           a_gains: np.ndarray | None) -> dict:
    if W_final is None or a_gains is None:
        return {"rhs_rms": 0.0, "rhs_max": 0.0}
    W_eff = W_final * np.outer(a_gains, a_gains)
    rhs = kuramoto_rhs_fast(0.0, theta, omega, K_SIGN, W_eff)
    return {
        "rhs_rms": float(np.sqrt(np.mean(rhs * rhs))),
        "rhs_max": float(np.max(np.abs(rhs))),
    }


def structural_metrics(W_final: np.ndarray | None, A_mask: np.ndarray, budget: float) -> dict:
    if W_final is None:
        return {
            "active_edges": 0,
            "mean_active_degree": 0.0,
            "row_sum_max": 0.0,
            "budget_violation_max": 0.0,
            "symmetry_error": 0.0,
        }
    active = W_final > 1e-6
    row_sums = W_final.sum(axis=1)
    fro = np.linalg.norm(W_final, "fro") + 1e-12
    return {
        "active_edges": int(active.sum() // 2),
        "mean_active_degree": float(active.sum(axis=1).mean()),
        "row_sum_max": float(row_sums.max()),
        "budget_violation_max": float(np.max(row_sums - budget)),
        "symmetry_error": float(np.linalg.norm(W_final - W_final.T, "fro") / fro),
    }


def run_hybrid_with_theta(A: np.ndarray, theta0: np.ndarray, omega: np.ndarray,
                          a_gains: np.ndarray, budget: float, prob, Y_param, X_var,
                          rngs: dict) -> dict:
    A_mask = A > 0
    start = time.time()
    from phase2_benchmark import project_frobenius, run_static_ode
    W0 = project_frobenius(A.astype(float), prob, Y_param, X_var, first_call=True)
    W0 = enforce_symmetric_support_budget(W0, A_mask, budget)
    theta_mid, W_frozen, nfev1, pc, au, energy_viol, final_L = run_hebbian_staggered(
        theta0, W0, omega, K_SIGN, a_gains, prob, Y_param, X_var,
        "frobenius", 50.0, 0.5, 1.0, 0.1, A_mask, budget)
    W_eff = W_frozen * np.outer(a_gains, a_gains)
    theta_final, nfev2 = run_static_ode(theta_mid, omega, K_SIGN, W_eff, 50.0)
    cut_raw, cut_rounded, cut_polished = random_hyperplane_rounding(
        theta_final, A, seeds=20, rng=rngs["rounding"])
    metrics = structural_metrics(W_frozen, A_mask, budget)
    metrics.update(final_dynamics_metrics(theta_final, omega, W_frozen, a_gains))
    return {
        "cut_raw": cut_raw,
        "cut_rounded": cut_rounded,
        "cut_polished": cut_polished,
        "runtime_sec": time.time() - start,
        "projection_calls": pc + 1,
        "ode_nfev": nfev1 + nfev2,
        "adaptive_updates": au,
        "energy_violations_raw": energy_viol,
        "final_L": final_L,
        "theta_final": theta_final,
        "W_final": W_frozen,
        **metrics,
    }


def run_phase2_method(method: str, A: np.ndarray, theta0: np.ndarray, omega: np.ndarray,
                      a_gains: np.ndarray, budget: float, prob, Y_param, X_var,
                      rngs: dict, fixed_support_mask: np.ndarray | None = None) -> dict:
    res = run_method(method, A, theta0, omega, a_gains, budget, prob, Y_param, X_var,
                     rngs, fixed_support_mask=fixed_support_mask)
    W_final = res[-1]
    metrics = structural_metrics(W_final, A > 0, budget)
    return {
        "cut_raw": res[0],
        "cut_rounded": res[1],
        "cut_polished": res[2],
        "runtime_sec": res[15],
        "projection_calls": res[16],
        "ode_nfev": res[17],
        "adaptive_updates": res[18],
        "energy_violations_raw": res[13],
        "final_L": res[14],
        "W_final": W_final,
        **metrics,
    }


def run_support_random(A: np.ndarray, theta0: np.ndarray, omega: np.ndarray,
                       a_gains: np.ndarray, budget: float, prob, Y_param, X_var,
                       rngs: dict, support: np.ndarray) -> dict:
    return run_phase2_method("random_budgeted", A, theta0, omega, a_gains, budget,
                             prob, Y_param, X_var, rngs, fixed_support_mask=support)


def maybe_stop(start_time: float, max_hours: float | None) -> bool:
    return max_hours is not None and (time.time() - start_time) >= max_hours * 3600.0


def run_heterogeneity(scope: str, max_hours: float | None) -> None:
    n_graph = 1 if scope == "smoke" else (3 if scope == "staged" else 10)
    n_method = 1 if scope == "smoke" else (5 if scope == "staged" else 10)
    sigmas = [0.0, 1.0] if scope == "smoke" else SIGMAS
    methods = ["static", "static_budgeted", "random_budgeted", "hybrid_frobenius",
               "hybrid_support_random_weights"]
    key_cols = ["graph_family", "graph_seed", "method_seed", "sigma", "method"]
    done = load_done(HET_CSV, key_cols)
    start = time.time()

    log(f"[heterogeneity] scope={scope} graph_seeds={n_graph} method_seeds={n_method} sigmas={sigmas}")
    for g_seed in range(n_graph):
        G = nx.erdos_renyi_graph(200, 0.05, seed=g_seed)
        A = graph_to_A(G)
        A_mask = A > 0
        n = A.shape[0]
        m = G.number_of_edges()
        mean_deg = 2.0 * m / n
        budget = mean_deg / 2.0
        prob, Y_param, X_var = make_frobenius_projector(n, A_mask, budget)

        for method_seed in range(n_method):
            for sigma in sigmas:
                rngs = make_rngs(g_seed, method_seed)
                theta0 = rngs["phases"].uniform(0, 2 * np.pi, n)
                omega = np.zeros(n)
                a_gains = rngs["amplitude"].lognormal(-0.5 * sigma ** 2, sigma, n)
                hybrid_result = None

                for method in methods:
                    key = ("sparse_er_200_p05", g_seed, method_seed, sigma, method)
                    if key in done:
                        continue
                    if maybe_stop(start, max_hours):
                        log("[heterogeneity] max-hours reached; stopping cleanly")
                        return

                    rngs_method = make_rngs(g_seed, method_seed)
                    if method == "hybrid_frobenius":
                        result = run_hybrid_with_theta(A, theta0, omega, a_gains, budget,
                                                       prob, Y_param, X_var, rngs_method)
                        hybrid_result = result
                    elif method == "hybrid_support_random_weights":
                        if hybrid_result is None:
                            hybrid_result = run_hybrid_with_theta(A, theta0, omega, a_gains,
                                                                   budget, prob, Y_param, X_var,
                                                                   rngs_method)
                        support = hybrid_result["W_final"] > 1e-6
                        result = run_support_random(A, theta0, omega, a_gains, budget,
                                                    prob, Y_param, X_var, rngs_method, support)
                    else:
                        result = run_phase2_method(method, A, theta0, omega, a_gains, budget,
                                                   prob, Y_param, X_var, rngs_method)

                    row = {
                        "graph_family": "sparse_er_200_p05",
                        "graph_seed": g_seed,
                        "method_seed": method_seed,
                        "sigma": sigma,
                        "method": method,
                        "n": n,
                        "m": m,
                        "mean_degree": mean_deg,
                        "budget_mode": "half_mean_deg",
                        "budget_value": budget,
                        "cut_raw": result["cut_raw"],
                        "cut_rounded": result["cut_rounded"],
                        "cut_polished": result["cut_polished"],
                        "runtime_sec": result["runtime_sec"],
                        "projection_calls": result["projection_calls"],
                        "ode_nfev": result["ode_nfev"],
                        "adaptive_updates": result["adaptive_updates"],
                        "rhs_rms": result.get("rhs_rms", np.nan),
                        "rhs_max": result.get("rhs_max", np.nan),
                        "active_edges": result["active_edges"],
                        "mean_active_degree": result["mean_active_degree"],
                        "row_sum_max": result["row_sum_max"],
                        "budget_violation_max": result["budget_violation_max"],
                        "symmetry_error": result["symmetry_error"],
                        "energy_violations_raw": result["energy_violations_raw"],
                    }
                    append_row(HET_CSV, row)
                    done.add(key)
                    log(f"[heterogeneity] g={g_seed} seed={method_seed} sigma={sigma} "
                        f"{method} cut={row['cut_polished']:.0f} t={row['runtime_sec']:.1f}s")


def run_greedy_row(A: np.ndarray, theta0: np.ndarray) -> dict:
    start = time.time()
    spins = np.sign(np.cos(theta0))
    spins[spins == 0] = 1
    spins = greedy_max_cut(A, spins)
    cut = get_cut_value(A, spins)
    return {
        "cut_raw": cut,
        "cut_rounded": cut,
        "cut_polished": cut,
        "runtime_sec": time.time() - start,
        "projection_calls": 0,
        "ode_nfev": 0,
        "adaptive_updates": 0,
        "budget_violation_max": 0.0,
        "symmetry_error": 0.0,
    }


def sdp_worker(A: np.ndarray) -> tuple[float, float]:
    return run_sdp_gw(A, np.random.default_rng(987654321))


def run_sdp_with_timeout(A: np.ndarray, timeout_sec: float) -> tuple[float, float, str]:
    with ProcessPoolExecutor(max_workers=1) as ex:
        fut = ex.submit(sdp_worker, A)
        try:
            bound, gw = fut.result(timeout=timeout_sec)
            status = "ok" if bound > 0 else "failed"
            return float(bound), float(gw), status
        except TimeoutError:
            fut.cancel()
            return 0.0, 0.0, "timeout"


def run_supporting(scope: str, max_hours: float | None) -> None:
    n_graph = 1 if scope == "smoke" else (3 if scope == "staged" else 10)
    n_method = 1 if scope == "smoke" else (5 if scope == "staged" else 10)
    if scope == "staged":
        families = ["gset_g1"]
        methods = ["static_budgeted", "random_budgeted", "hebbian_frobenius",
                   "hybrid_frobenius", "greedy"]
        budgets = [("half_mean_deg", 0.5)]
    else:
        families = ["dense_er_200_p15", "gset_g1"] if scope != "smoke" else ["dense_er_200_p15"]
        methods = ["static", "static_budgeted", "random_budgeted", "hebbian_frobenius",
                   "hybrid_frobenius", "topology_scrambled", "greedy"]
        budgets = [("half_mean_deg", 0.5), ("mean_deg", 1.0)]
    key_cols = ["graph_family", "graph_seed", "method_seed", "budget_mode", "method"]
    done = load_done(BENCH_CSV, key_cols)
    start = time.time()

    log(f"[supporting] scope={scope} families={families}")
    for family in families:
        graph_seeds = [0] if family == "gset_g1" else list(range(n_graph))
        if family == "gset_g1" and not Path("G1").exists():
            log("[supporting] G1 missing; skipping")
            continue
        for g_seed in graph_seeds:
            if family == "gset_g1":
                G = load_gset("G1")
                best_known = G1_BEST_KNOWN
                exact_status = "best_known"
            else:
                G, _, _, _, best_known = build_graph(family, g_seed, np.random.default_rng(g_seed))
                exact_status = "sdp_or_timeout"
            A = graph_to_A(G)
            A_mask = A > 0
            n = A.shape[0]
            m = G.number_of_edges()
            mean_deg = 2.0 * m / n
            sdp_bound = gw_cut = 0.0
            sdp_status = "not_run"
            if family != "gset_g1" and scope != "smoke":
                sdp_bound, gw_cut, sdp_status = run_sdp_with_timeout(A, timeout_sec=180.0)

            for budget_mode, mult in budgets:
                budget = mean_deg * mult
                prob, Y_param, X_var = make_frobenius_projector(n, A_mask, budget)
                for method_seed in range(n_method):
                    rngs = make_rngs(g_seed, method_seed)
                    theta0 = rngs["phases"].uniform(0, 2 * np.pi, n)
                    omega = np.zeros(n)
                    a_gains = np.ones(n)
                    for method in methods:
                        key = (family, g_seed, method_seed, budget_mode, method)
                        if key in done:
                            continue
                        if maybe_stop(start, max_hours):
                            log("[supporting] max-hours reached; stopping cleanly")
                            return
                        if method == "greedy":
                            result = run_greedy_row(A, theta0)
                        else:
                            result = run_phase2_method(method, A, theta0, omega, a_gains,
                                                       budget, prob, Y_param, X_var, rngs)
                        ceiling = best_known if family == "gset_g1" else sdp_bound
                        ceiling_label = "known_best_11624" if family == "gset_g1" else sdp_status
                        row = {
                            "graph_family": family,
                            "graph_seed": g_seed,
                            "method_seed": method_seed,
                            "budget_mode": budget_mode,
                            "method": method,
                            "n": n,
                            "m": m,
                            "mean_degree": mean_deg,
                            "budget_value": budget,
                            "cut_raw": result["cut_raw"],
                            "cut_rounded": result["cut_rounded"],
                            "cut_polished": result["cut_polished"],
                            "runtime_sec": result["runtime_sec"],
                            "projection_calls": result["projection_calls"],
                            "ode_nfev": result["ode_nfev"],
                            "adaptive_updates": result["adaptive_updates"],
                            "budget_violation_max": result["budget_violation_max"],
                            "symmetry_error": result["symmetry_error"],
                            "sdp_bound": sdp_bound,
                            "gw_cut": gw_cut,
                            "ceiling_value": ceiling,
                            "ceiling_label": ceiling_label,
                            "ratio_to_ceiling": result["cut_polished"] / ceiling
                            if ceiling and not math.isnan(ceiling) else 0.0,
                            "exact_status": exact_status,
                        }
                        append_row(BENCH_CSV, row)
                        done.add(key)
                        log(f"[supporting] {family} g={g_seed} seed={method_seed} "
                            f"{budget_mode} {method} cut={row['cut_polished']:.0f} "
                            f"t={row['runtime_sec']:.1f}s")

    if not Path("G14").exists():
        log("[supporting] G14 not present in repo root; skipped to avoid burning time.")


def summarize_mean_std(df: pd.DataFrame, group_cols: list[str], value: str) -> pd.DataFrame:
    return df.groupby(group_cols)[value].agg(["mean", "std", "count"]).reset_index()


def analyze() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    if not HET_CSV.exists():
        raise SystemExit("No heterogeneity CSV to analyze")
    het = pd.read_csv(HET_CSV)
    if len(het) == 0:
        raise SystemExit("Empty heterogeneity CSV")

    base = (het[het["sigma"] == 0.0]
            .groupby(["graph_seed", "method_seed", "method"])["cut_polished"]
            .mean().rename("method_sigma0"))
    het = het.join(base, on=["graph_seed", "method_seed", "method"])
    het["cut_ratio_vs_method_sigma0"] = het["cut_polished"] / het["method_sigma0"]

    slope_rows = []
    complete_keys = []
    for (method, g_seed, method_seed), g in het.groupby(["method", "graph_seed", "method_seed"]):
        observed_sigmas = sorted(g["sigma"].unique())
        if 0.0 in observed_sigmas and len(observed_sigmas) >= 2:
            gm = g.groupby("sigma")["cut_ratio_vs_method_sigma0"].mean().reindex(observed_sigmas)
            slope = np.polyfit(np.array(observed_sigmas), gm.to_numpy(), 1)[0]
            slope_rows.append({
                "method": method,
                "graph_seed": g_seed,
                "method_seed": method_seed,
                "n_sigmas": len(observed_sigmas),
                "slope_ratio_per_sigma": float(slope),
            })
        if set(SIGMAS).issubset(set(observed_sigmas)):
            complete_keys.append((method, g_seed, method_seed))
    slopes = pd.DataFrame(slope_rows)
    slopes.to_csv(SLOPE_CSV, index=False)

    complete_index = pd.MultiIndex.from_tuples(
        complete_keys, names=["method", "graph_seed", "method_seed"])
    row_index = pd.MultiIndex.from_frame(het[["method", "graph_seed", "method_seed"]])
    het_complete = het[row_index.isin(complete_index)].copy()
    if het_complete.empty:
        het_complete = het

    plot_degradation(het_complete)
    plot_convergence(het_complete)
    write_report(het_complete, slopes, total_het_rows=len(het))
    log(f"[analyze] wrote {SLOPE_CSV}, plots, and {WRITEUP_PATH}")


def plot_degradation(het: pd.DataFrame) -> None:
    methods = ["static", "static_budgeted", "random_budgeted", "hybrid_frobenius",
               "hybrid_support_random_weights"]
    plt.figure(figsize=(8, 5))
    for method in methods:
        sub = het[het["method"] == method]
        if sub.empty:
            continue
        stats = sub.groupby("sigma")["cut_ratio_vs_method_sigma0"].agg(["mean", "std"]).reset_index()
        x = stats["sigma"].to_numpy()
        y = stats["mean"].to_numpy()
        e = stats["std"].fillna(0.0).to_numpy()
        plt.plot(x, y, marker="o", label=method)
        plt.fill_between(x, y - e, y + e, alpha=0.15)
    plt.axhline(1.0, color="black", linewidth=0.8)
    plt.xlabel("Amplitude heterogeneity sigma")
    plt.ylabel("Cut ratio vs same-method sigma=0")
    plt.title("Amplitude-Heterogeneity Degradation")
    plt.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "stability_degradation_cut.png", dpi=160)
    plt.close()


def plot_convergence(het: pd.DataFrame) -> None:
    plt.figure(figsize=(8, 5))
    for method in ["static", "static_budgeted", "hybrid_frobenius"]:
        sub = het[het["method"] == method]
        if sub.empty or sub["rhs_rms"].isna().all():
            continue
        stats = sub.groupby("sigma")["rhs_rms"].agg(["mean", "std"]).reset_index()
        x = stats["sigma"].to_numpy()
        y = stats["mean"].to_numpy()
        e = stats["std"].fillna(0.0).to_numpy()
        plt.plot(x, y, marker="o", label=method)
        plt.fill_between(x, y - e, y + e, alpha=0.15)
    plt.xlabel("Amplitude heterogeneity sigma")
    plt.ylabel("Final RHS RMS")
    plt.title("Convergence Under Heterogeneity")
    plt.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "stability_convergence.png", dpi=160)
    plt.close()


def write_report(het: pd.DataFrame, slopes: pd.DataFrame, total_het_rows: int | None = None) -> None:
    cut_stats = summarize_mean_std(
        het, ["sigma", "method"], "cut_polished").sort_values(["sigma", "method"])
    ratio_stats = summarize_mean_std(
        het, ["sigma", "method"], "cut_ratio_vs_method_sigma0").sort_values(["sigma", "method"])
    if slopes.empty:
        slope_stats = pd.DataFrame(columns=["method", "mean", "std", "count"])
    else:
        slope_stats = summarize_mean_std(
            slopes, ["method"], "slope_ratio_per_sigma").sort_values("method")

    static_slope = slope_stats[slope_stats["method"] == "static_budgeted"]["mean"]
    hybrid_slope = slope_stats[slope_stats["method"] == "hybrid_frobenius"]["mean"]
    if len(static_slope) and len(hybrid_slope):
        delta = float(hybrid_slope.iloc[0] - static_slope.iloc[0])
        if delta > 0.002:
            verdict = "holds preliminarily"
        elif delta > -0.002:
            verdict = "holds partially / inconclusive"
        else:
            verdict = "does not hold"
    else:
        verdict = "incomplete"

    g1_stats = pd.DataFrame()
    if BENCH_CSV.exists():
        bench = pd.read_csv(BENCH_CSV)
        g1 = bench[(bench["graph_family"] == "gset_g1") &
                   (bench["budget_mode"] == "half_mean_deg")]
        if not g1.empty:
            g1_stats = (g1.groupby("method")["cut_polished"]
                          .agg(["mean", "std", "count"]).reset_index()
                          .sort_values("method"))
            g1_means = g1_stats.set_index("method")["mean"]
            if "hybrid_frobenius" in g1_means and "random_budgeted" in g1_means:
                if g1_means["hybrid_frobenius"] < g1_means["random_budgeted"] and verdict != "incomplete":
                    verdict = "does not hold"

    caps = []
    staged_expected = 3 * 5 * len(SIGMAS) * 5
    full_expected = 10 * 10 * len(SIGMAS) * 5
    observed_rows = total_het_rows if total_het_rows is not None else len(het)
    if observed_rows < staged_expected:
        caps.append(f"heterogeneity rows {observed_rows} of staged diagnostic target {staged_expected}; "
                    f"headline tables use {len(het)} rows from complete four-sigma paired cells")
    caps.append(f"full robust 10x10 heterogeneity matrix ({full_expected} rows) was not run")
    if not Path("G14").exists():
        caps.append("G14 skipped because it was not present locally")

    def table(df: pd.DataFrame) -> str:
        if df.empty:
            return "(no complete rows)"
        return df.to_string(index=False, float_format=lambda x: f"{x:.5f}")

    text = [
        "# Stability Experiment Writeup",
        "",
        f"Verdict: **{verdict}**.",
        "",
        "This is a preliminary, competitive, hardware-aware stability run. It is not a validation or proof.",
        "",
        "## Degradation Slopes",
        "",
        "```text",
        table(slope_stats),
        "```",
        "",
        "## Cut Quality Mean +/- Std",
        "",
        "```text",
        table(cut_stats),
        "```",
        "",
        "## Normalized Degradation Mean +/- Std",
        "",
        "```text",
        table(ratio_stats),
        "```",
        "",
        "## G1 Half-Mean-Degree Falsification",
        "",
        "```text",
        table(g1_stats),
        "```",
        "",
        "## Discipline Notes",
        "",
        "- Headline comparison is `hybrid_frobenius` vs `static_budgeted`; raw `static` is context.",
        "- The plotted degradation ratio normalizes each method to its own sigma=0 paired baseline.",
        "- G1 SDP was not run; the ceiling reference is known best cut 11624.",
        "- Results are reported as mean +/- std, not best-of.",
    ]
    if caps:
        text.extend(["", "## Runtime Caps", ""])
        text.extend(f"- {cap}" for cap in caps)
    WRITEUP_PATH.write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["smoke", "staged", "full", "analyze"])
    parser.add_argument("--max-hours", type=float, default=None)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(exist_ok=True)
    if args.mode == "analyze":
        analyze()
        return
    log(f"[start] mode={args.mode} max_hours={args.max_hours}")
    run_heterogeneity(args.mode, args.max_hours)
    run_supporting(args.mode, args.max_hours)
    analyze()


if __name__ == "__main__":
    main()
