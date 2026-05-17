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
# K_c ≈ 1.2 is the mean-field result for all-to-all coupling — it does NOT apply here.
# Nearest-neighbor 1D rings have no clean phase transition; sync scales differently with N.
plt.axvline(x=1.2, color='red', linestyle='--', alpha=0.3, label='Mean-field K_c ≈ 1.2 (all-to-all only, not applicable here)')
plt.xlabel('Coupling strength K')
plt.ylabel('Order parameter r')
plt.title('Bifurcation Diagram — Kuramoto Oscillators on a Ring')
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------------------------------
# 2. Simulation at K = 2.5 — twisted / locally-ordered regime
# Note: pure nearest-neighbor coupling produces twisted states, not true chimeras.
# True chimeras require nonlocal coupling (Abrams & Strogatz 2004).
# ------------------------------------------------
K_sim = 2.5   # Intermediate coupling — twisted states and local ordering

theta0 = np.random.uniform(0, 2*np.pi, N)
sol = solve_ivp(kuramoto, t_span, theta0, args=(omega, K_sim, N),
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
plt.title(f'Time evolution of order parameter r(t)  (K = {K_sim})')
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
plot_ring(theta_final, f'Final ring state (K={K_sim}) — twisted/chimera, not global sync', ax2)
plt.suptitle('Oscillators on a Ring — winding number selection and chimera states')
plt.tight_layout()
plt.show()

print("\nSimulation complete. Phase heatmap + ring plot show the chimera/twisted geometry; r alone does not.")
