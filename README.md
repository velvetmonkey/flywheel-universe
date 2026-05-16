# Flywheel Universe

> embrace the cycles.
> become the chain.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/velvetmonkey/flywheel-universe/blob/main/kuramoto_ring.ipynb)

> All code below was reverse-engineered from this image.

![Banner](readme_assets/banner.jpg)

## What this is

120 oscillators arranged in a ring. Each one only talks to its immediate left and right neighbours. That constraint changes everything.

The standard Kuramoto model assumes all-to-all coupling — every oscillator feels every other. This produces a clean mean-field phase transition: below critical coupling K_c, incoherence. Above it, global synchrony. Textbook.

Nearest-neighbour coupling on a ring is different. The system finds **traveling wave solutions** and **chimera states** — islands of coherent oscillators floating in an incoherent sea. The global order parameter r stays low even when the system is highly organised. r is the wrong metric. The geometry is the point.

## What you'll see

- **Bifurcation diagram** — r vs K across the full range. Noisy, non-monotonic. That's the topology, not noise.
- **Phase portrait** — the order parameter R traced in the complex plane. Not a clean spiral to the edge of the unit circle. A wander with structure.
- **Phase heatmap** — where chimera states live. Flat bands (locked) floating in diagonal drift (incoherent).
- **The crystal chain** — oscillators mapped to the unit circle, connected to neighbours. Tangled web becomes coherent ring.

![Output 1](readme_assets/output1.png)

![Output](readme_assets/output.png)

![Attractor geometry comparison panel](readme_assets/attractor_geometry_comparison_panel.png)

![LPA heteroclinic cycle chain 3D](readme_assets/lpa_heteroclinic_cycle_chain_3d.png)

![McKay affine E8 rep dimensions](readme_assets/mckay_affine_e8_rep_dimensions.png)

![McKay affine E8 rep dimensions (SVG)](readme_assets/mckay_affine_e8_rep_dimensions.svg)

![Stuart-Landau normal form (animated)](readme_assets/stuart_landau_normal_form_animated.gif)

![Stuart-Landau normal form (static)](readme_assets/stuart_landau_normal_form_static.png)

![Oph](readme_assets/oph.png)

Observer Patch Holography — by [FloatingPragma](https://github.com/FloatingPragma/observer-patch-holography)

[Interactive demo: Prime Inertia Engine →](https://codepen.io/DULA2025/pen/yyVejbr) — by [DULA2025](https://github.com/DULA2025/prime-inertia-engine)

## Why it matters

The hard problem of consciousness asks what collective order feels like from inside one oscillator. This simulation is the geometry of that question, made executable.

Chimera states — partial coherence, structured incoherence — may be the most honest mathematical model of a mind anyone has written down.

## Run it

```bash
pip install numpy scipy matplotlib
python kuramoto_ring.py
```

## Related

- Witness theory: `github.com/velvetmonkey/witness-descent`
- Dula's OPH toroidal emulator: `codepen.io/DULA2025/pen/LERMgxv`

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

def kuramoto(t, theta, omega, K, N):
    """
    Right-hand side of Kuramoto ODE on a ring.
    theta: array of phases [N]
    """
    dtheta = np.zeros(N)
    for i in range(N):
        ip = (i + 1) % N
        im = (i - 1) % N
        # Coupling from left and right neighbors only
        dtheta[i] = omega[i] + K * (np.sin(theta[ip] - theta[i]) + np.sin(theta[im] - theta[i]))
    return dtheta


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
plt.axvline(x=1.2, color='red', linestyle='--', alpha=0.6, label='Approximate critical K_c')
plt.xlabel('Coupling strength K')
plt.ylabel('Order parameter r')
plt.title('Bifurcation Diagram — Kuramoto Oscillators on a Ring')
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------------------------------
# 2. Detailed simulation at supercritical K (harmony!)
# ------------------------------------------------
K_super = 2.5   # Well above critical coupling → strong synchronization

theta0 = np.random.uniform(0, 2*np.pi, N)
sol = solve_ivp(kuramoto, t_span, theta0, args=(omega, K_super, N),
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
plt.title(f'Time evolution of synchronization (K = {K_super})')
plt.grid(True, alpha=0.5)
plt.ylim(0, 1.05)
plt.tight_layout()
plt.show()

# B. Phase portrait: Order parameter in complex plane
#     → This is the spiraling/winding-path "portal" visualization
plt.figure(figsize=(7, 7))
plt.plot(R_real, R_imag, 'b-', linewidth=1.8, alpha=0.85, label='Trajectory')
plt.plot(R_real[0], R_imag[0], 'go', markersize=8, label='Initial (incoherent)')
plt.plot(R_real[-1], R_imag[-1], 'ro', markersize=8, label='Final (synchronized attractor)')
plt.xlabel('Re(R)')
plt.ylabel('Im(R)')
plt.title('Phase Portrait of the Order Parameter\nSpiraling into the synchronized attractor')
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
plt.title('Phase Evolution on the Ring\n(Random → Coherent synchronization)')
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
plot_ring(theta_final, f'Final: Synchronized (K={K_super}) — "Become the Chain"', ax2)
plt.suptitle('Oscillators on a Ring — "Embrace the Cycles, Become the Chain"')
plt.tight_layout()
plt.show()

print("\nSimulation complete. The complex-plane spiral and the phase heatmap/ring plots best evoke the 'portal / vortex / harmony threshold' geometry from the image.")
```
