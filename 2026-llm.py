exec(r'''
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")
rng = np.random.default_rng(7)

def unit_vectors(n, d):
    x = rng.normal(size=(n, d))
    return x / np.linalg.norm(x, axis=1, keepdims=True)

# ------------------------------------------------------------
# Figure 1: quasi-orthogonal coexistence
# ------------------------------------------------------------

fig1, ax1 = plt.subplots(figsize=(5.3, 4.2))

eps = 0.25
dims_capacity = [128, 192, 256, 320]
alphas = np.linspace(0.05, 0.33, 8)
trials_capacity = 20

for d in dims_capacity:
    probs = []

    for alpha in alphas:
        N = max(3, int(round(np.exp(alpha * d * eps**2))))
        good = 0

        for trial in range(trials_capacity):
            V = unit_vectors(N, d)
            G = np.abs(V @ V.T)
            np.fill_diagonal(G, 0.0)
            good += int(G.max() <= eps)

        probs.append(good / trials_capacity)

    ax1.plot(alphas, probs, "o-", label=fr"$d_{{eff}}={d}$")

ax1.set_xlabel(r"$\alpha$ in $N=\exp(\alpha d_{eff}\epsilon^2)$")
ax1.set_ylabel(r"Pr(all pairwise overlaps $\leq \epsilon$)")
ax1.set_title(r"Quasi-orthogonal coexistence, $\epsilon=0.25$")
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, alpha=0.25)
ax1.legend(frameon=False, fontsize=8)

fig1.tight_layout()
fig1.savefig("quasi_orthogonal_coexistence.pdf", bbox_inches="tight")
fig1.savefig("quasi_orthogonal_coexistence.png", dpi=240, bbox_inches="tight")


# ------------------------------------------------------------
# Figure 2: readout interference
# ------------------------------------------------------------

fig2, ax2 = plt.subplots(figsize=(5.3, 4.2))

ks = np.array([1, 2, 4, 8, 16, 32, 64, 128])
dims_readout = [128, 256, 512]
trials_readout = 5000

for d in dims_readout:
    measured = []

    for k in ks:
        overlaps = rng.normal(0.0, 1 / np.sqrt(d), size=(trials_readout, k))
        signs = rng.choice([-1.0, 1.0], size=(trials_readout, k))
        interference = np.sum(signs * overlaps, axis=1)
        measured.append(np.std(interference, ddof=1))

    measured = np.array(measured)
    theory = np.sqrt(ks / d)

    line, = ax2.plot(ks, measured, "o-", label=fr"simulation $d_{{eff}}={d}$")
    ax2.plot(ks, theory, "--", color=line.get_color(), alpha=0.75)

ax2.set_xscale("log", base=2)
ax2.set_yscale("log")
ax2.set_xlabel(r"active distractors $k$")
ax2.set_ylabel(r"interference std.")
ax2.set_title(r"Readout noise $\sim \sqrt{k/d_{eff}}$")
ax2.set_xlim(1, 128)
ax2.set_ylim(0.03, 1.2)
ax2.grid(True, which="both", alpha=0.25)
ax2.legend(frameon=False, fontsize=8)

fig2.tight_layout()
fig2.savefig("readout_interference.pdf", bbox_inches="tight")
fig2.savefig("readout_interference.png", dpi=240, bbox_inches="tight")

print("saved:")
print("  quasi_orthogonal_coexistence.pdf")
print("  quasi_orthogonal_coexistence.png")
print("  readout_interference.pdf")
print("  readout_interference.png")

plt.show()
''')
