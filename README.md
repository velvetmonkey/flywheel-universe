# Flywheel Universe

> embrace the cycles.
> become the chain.

**Notebooks** (run any in-browser, no install):

| | Notebook | Open |
|---|---|---|
| 1 | `kuramoto_ring.ipynb` — 1D nearest-neighbour ring | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/velvetmonkey/flywheel-universe/blob/main/kuramoto_ring.ipynb) |
| 2 | `meteorology.ipynb` — structured frequency coarsening | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/velvetmonkey/flywheel-universe/blob/main/meteorology.ipynb) |
| 3 | `kuramoto_universe.ipynb` — hierarchical plastic oscillator fluid | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/velvetmonkey/flywheel-universe/blob/main/kuramoto_universe.ipynb) |

> Inspired by this image, which suggested the topology — the math wasn't in the picture, the picture suggested which math.

![Banner](readme_assets/banner.jpg)

## What this is

120 oscillators arranged in a ring. Each one only talks to its immediate left and right neighbours. That single constraint changes the physics completely.

The standard Kuramoto model assumes all-to-all coupling — every oscillator feels every other. This produces a clean mean-field phase transition with a well-defined critical K_c. Nearest-neighbour coupling on a ring is different. There is no such transition. The system selects a **winding number** based on initial conditions and locks there permanently.

## What actually happens

**Winding number selection.** At any K, the ring settles into a twisted state where phases advance uniformly around the circle q times: θ_i ≈ 2πqi/N. These states have winding number q = 0, ±1, ±2... and are stable fixed points in a rotating frame. Because the phases are spread evenly around the circle by construction, the global order parameter r stays low — not because the system failed to synchronise, but because it found a different kind of order entirely.

**r is the wrong metric — for a topological reason.** r measures global phase alignment. Twisted states have zero global alignment by design. The right question is which winding number was selected, and why.

**Chimera states.** At intermediate K, a local order parameter (averaged over a neighbourhood of oscillators) reveals coherent islands — r_local up to ~0.85 in a cluster — floating in an otherwise incoherent ring. The global r (~0.1) hides this completely. The phase heatmap shows it as flat locked bands drifting in a sea of diagonal incoherence.

**Mean-field intuitions fail here.** K_c ≈ 1.2 is the mean-field result and does not apply. Each oscillator talks to 2 out of 119 others. Global sync on a 1D ring would require K ~ O(N), and even then the system may prefer a twisted state.

## Visualisations

- **Bifurcation diagram** — r vs K. No transition. Noise in the 0.05–0.30 band. The story is in the winding number, not r.
- **Phase heatmap** — chimera states visible as locked horizontal bands in drifting incoherence.
- **Complex-plane portrait** — z(t) traces a quasi-periodic orbit near the origin. Not a spiral to the unit circle.
- **The crystal chain** — oscillators on the unit circle connected to neighbours. Twisted states show as smooth colour gradients winding around the ring.

![Phase portrait](readme_assets/03_phase_portrait.png)

![Phase heatmap](readme_assets/04_phase_heatmap.png)

![K-sweep rings — twisted states at increasing K](readme_assets/06_K_sweep_rings.png)

![Local order parameter — chimera islands hidden by global r](readme_assets/07_local_order.png)

![r(t) traces across K](readme_assets/08_rg_traces.png)

![Output 1](readme_assets/output1.png)

![Output](readme_assets/output.png)

![Attractor geometry comparison panel](readme_assets/attractor_geometry_comparison_panel.png)

![LPA heteroclinic cycle chain 3D](readme_assets/lpa_heteroclinic_cycle_chain_3d.png)

![McKay affine E8 rep dimensions](readme_assets/mckay_affine_e8_rep_dimensions.png)

![McKay affine E8 rep dimensions (SVG)](readme_assets/mckay_affine_e8_rep_dimensions.svg)

![Stuart-Landau normal form (animated)](readme_assets/stuart_landau_normal_form_animated.gif)

![Stuart-Landau normal form (static)](readme_assets/stuart_landau_normal_form_static.png)

![Oph](readme_assets/oph.png)

![Meteorology — Gaussian omega control. Locked/drifting heatmap at K=0.8, sigma=1.5, N=200, T=200, testing whether the coarsening pattern survives without structured frequencies](readme_assets/dl_01.png)

*Meteorology — Gaussian ω control. Tests whether the lock/drift coarsening pattern is generic or structure-dependent; even with unstructured natural frequencies, persistent locked and drifting bands form.*

![Meteorology — unconstrained Hebbian: network at t=0 with lock/drift boundary, then weight histograms at t=50 and t=100 showing graph fully densified to 100% pairs active](readme_assets/dl_02.png)

*Meteorology — unconstrained Hebbian (cell 4). Without a conservation law the Hebbian rule densifies the graph completely: by t=50, 100% of pairs are active. No preferential attachment, no scale-free structure. The negative result is the result.*

![Meteorology — conserved Hebbian: network snapshots at t=0, 50, 100 with cyan lock/drift boundary nodes and degree distributions](readme_assets/dl_03.png)

*Meteorology — conserved Hebbian (cell 5). A per-node connection budget slows densification but does not produce scale-free topology. No power-law tail emerges. Three ingredients are required: kinematics, feedback, and dissipation.*

![Kuramoto Universe — filament evolution across t=0, 999, 1999, 2999, 3999, 4999; global R rises from 0.015 to 0.772](readme_assets/dl_04.png)

*Universe — filament evolution across the full run. Cluster colour tracks the dominant phase basin (rainbow → blue → red → orange → blue) as global R climbs from 0.015 to 0.772.*

![Global order parameter R(t) over t=0 to 5000 — saturates near 0.77 with bounded fluctuation](readme_assets/dl_05.png)

*Universe — global R(t) saturates near 0.77 by t≈100 and holds for the full t=5000 run, with one brief excursion near t≈3500.*

![Final-state network with 216 edges and mean degree 4.5, plus degree distribution histogram peaked at 4–6](readme_assets/dl_06.png)

*Universe — final network state from run 10. 216 undirected edges, mean degree 4.5, degree distribution peaked near 5.*

![Anomaly cell R(t) traces — cyan lines are four persistent boundary riders stuck below R=0.5; grey lines are cells absorbed into the main cluster](readme_assets/dl_07.png)

*Universe — four cells (6, 70, 86, 92) stayed below R=0.5 for the entire t=5000 run while the rest locked into the bulk. 100% fraction_below_0.5 throughout — these are the persistent boundary riders flagged in the commit log.*

## Why it matters

Chimera states — partial coherence, structured incoherence, winding-number selection — are well-studied (Kuramoto-Battogtokh 2002, Abrams-Strogatz 2004) and have speculative links to neural dynamics. This repo makes them executable and visible. The connection to consciousness is a direction, not yet a result.

## Run it

```bash
pip install numpy scipy matplotlib
python kuramoto_ring.py
```

## Related

**Own work:**

- Universe feature/bug task: [zenodo.org/records/20179566](https://zenodo.org/records/20179566)

**Potential collaborators — independent work, listed for context, not endorsement:**

- FloatingPragma — observer-patch holography (theory): [github.com/FloatingPragma/observer-patch-holography](https://github.com/FloatingPragma/observer-patch-holography)
- Dula (DULA2025) — OPH toroidal emulator (visualization): [codepen.io/DULA2025/pen/LERMgxv](https://codepen.io/DULA2025/pen/LERMgxv)
- Dula (DULA2025) — Prime Inertia Engine: [codepen.io/DULA2025/pen/yyVejbr](https://codepen.io/DULA2025/pen/yyVejbr) · [github.com/DULA2025/prime-inertia-engine](https://github.com/DULA2025/prime-inertia-engine)
- dmytronic — preprint *Logic Instantiated in Time is Ontically Experienced*: [zenodo.org/records/20233225](https://zenodo.org/records/20233225)

*These are researchers whose work intersects thematically. Inclusion is for context, not endorsement; connections are speculative and the work is unverified by this author.*

## Source

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.cm as cm

# ================================================
# Kuramoto Model on a Ring (1D periodic chain)
# Nearest-neighbor coupling — exactly as requested
# ================================================

def kuramoto(t, theta, omega, K, N=None):
    """
    Right-hand side of Kuramoto ODE on a ring (nearest-neighbor coupling).
    Vectorized: np.roll handles the periodic ±1 neighbour shifts in one shot,
    ~10x faster than the python for-loop on the bifurcation sweep.
    N is unused (kept for backward-compatible signature).
    """
    return omega + K * (np.sin(np.roll(theta, -1) - theta) + np.sin(np.roll(theta, 1) - theta))


def order_parameter(theta):
    """Compute magnitude of complex order parameter r"""
    z = np.mean(np.exp(1j * theta))
    return np.abs(z)


# Parameters (tuned for clear visuals + reasonable runtime)
N = 120                    # Number of oscillators on the ring
np.random.seed(42)
omega = np.random.normal(0.0, 0.65, N)   # Natural frequencies (Gaussian spread)

t_span = (0, 150)
t_eval_long = np.linspace(0, 150, 3000)

# ------------------------------------------------
# 1. Bifurcation diagram: order parameter r vs K
# ------------------------------------------------
print("Computing bifurcation diagram...")
K_values = np.linspace(0.0, 4.0, 35)
r_final = []

for K in K_values:
    theta0 = np.random.uniform(0, 2*np.pi, N)
    sol = solve_ivp(kuramoto, t_span, theta0, args=(omega, K, N),
                    method='RK45', rtol=1e-6, atol=1e-8)
    theta_final = sol.y[:, -1]
    r = order_parameter(theta_final)
    r_final.append(r)
    print(f"K = {K:.2f}  →  r = {r:.3f}")

# Plot bifurcation
plt.figure(figsize=(9, 5))
plt.plot(K_values, r_final, 'o-', linewidth=2, markersize=4)
# Pedagogical contrast: overlay where mean-field theory *would* predict K_c.
# Nothing actually happens at K=1.2 on a nearest-neighbour ring — the line is
# kept to make the absence of a transition visible against the textbook claim.
plt.axvline(x=1.2, color='red', linestyle='--', alpha=0.3, label='Mean-field K_c (does not apply here)')
plt.xlabel('Coupling strength K')
plt.ylabel('Order parameter r')
plt.title('Bifurcation Diagram — Kuramoto Oscillators on a Ring')
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------------------------------
# 2. Simulation at K = 2.5 — chimera / twisted regime
# ------------------------------------------------
K_chimera = 2.5   # Intermediate coupling — chimera states and twisted attractors live here

theta0 = np.random.uniform(0, 2*np.pi, N)
sol = solve_ivp(kuramoto, t_span, theta0, args=(omega, K_chimera, N),
                method='RK45', dense_output=True, rtol=1e-6, atol=1e-8)

t_plot = np.linspace(0, 80, 1200)

# Compute order parameter trajectory
R_real = []
R_imag = []
r_time = []

for t in t_plot:
    theta_t = sol.sol(t)
    z = np.mean(np.exp(1j * theta_t))
    R_real.append(z.real)
    R_imag.append(z.imag)
    r_time.append(np.abs(z))

R_real = np.array(R_real)
R_imag = np.array(R_imag)

# ------------------------------------------------
# 3. Plots (exactly what you asked for)
# ------------------------------------------------

# A. Time evolution of order parameter r(t)
plt.figure(figsize=(8, 4))
plt.plot(t_plot, r_time, linewidth=2.5)
plt.xlabel('Time')
plt.ylabel('Order parameter r(t)')
plt.title(f'Time evolution of order parameter r(t)  (K = {K_chimera})')
plt.grid(True, alpha=0.5)
plt.ylim(0, 1.05)
plt.tight_layout()
plt.show()

# B. Phase portrait: Order parameter in complex plane
#     → This is the spiraling/winding-path "portal" visualization
plt.figure(figsize=(7, 7))
plt.plot(R_real, R_imag, 'b-', linewidth=1.8, alpha=0.85, label='Trajectory')
plt.plot(R_real[0], R_imag[0], 'go', markersize=8, label='Initial (incoherent)')
plt.plot(R_real[-1], R_imag[-1], 'ro', markersize=8, label='Final (chimera/twisted state)')
plt.xlabel('Re(R)')
plt.ylabel('Im(R)')
plt.title('Phase Portrait of the Order Parameter\nQuasi-periodic orbit — not a spiral to sync')
plt.axis('equal')
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()

# C. Phase evolution heatmap (shows "become the chain")
theta_sol = np.array([sol.sol(t) for t in t_plot]).T   # shape (N, len(t_plot))
phases_wrapped = theta_sol % (2 * np.pi)

plt.figure(figsize=(11, 6))
plt.pcolormesh(t_plot, np.arange(N), phases_wrapped, cmap='twilight_shifted', shading='auto')
plt.colorbar(label='Phase θ_i (mod 2π)')
plt.xlabel('Time')
plt.ylabel('Oscillator index i (along the ring)')
plt.title('Phase Evolution on the Ring\n(Random → chimera bands floating in drift)')
plt.tight_layout()
plt.show()

# D. Ring visualization at beginning and end (the "crystal chain")
def plot_ring(theta, title, ax):
    x = np.cos(theta)
    y = np.sin(theta)
    # Scatter oscillators colored by phase
    sc = ax.scatter(x, y, c=theta, cmap='twilight_shifted', s=60, edgecolor='black', linewidth=0.5)
    # Connect nearest neighbors to show the "chain"
    for i in range(N):
        ax.plot([x[i], x[(i+1)%N]], [y[i], y[(i+1)%N]], 'k-', alpha=0.25, linewidth=1)
    ax.set_title(title)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    return sc

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
theta_init = theta0
theta_final = sol.sol(t_plot[-1])

plot_ring(theta_init, 'Initial: Incoherent phases', ax1)
plot_ring(theta_final, f'Final ring state (K={K_chimera}) — twisted/chimera, not global sync', ax2)
plt.suptitle('Oscillators on a Ring — winding number selection and chimera states')
plt.tight_layout()
plt.show()

print("\nSimulation complete. Phase heatmap + ring plot show the chimera/twisted geometry; r alone does not.")
```
