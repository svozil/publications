# -*- coding: utf-8 -*-
import sympy as sp
from sympy import symbols, exp, sin, cos, diff, simplify, diag, Matrix, zeros, eye, pi

# Define symbols
r, theta, phi, Nv, G = symbols('r theta phi N_v G', real=True, positive=True)
a = (1 - Nv)**(sp.Rational(2, 3))  # a = (1 - Nv)^(2/3)
f = a * exp(-r**2)  # f(r) = a * exp(-r^2)

# Define coordinates: [t, r, theta, phi]
coords = [symbols('t'), r, theta, phi]

# Define metric tensor g_{mu nu} (diagonal)
g = Matrix([
    [-1, 0, 0, 0],
    [0, f, 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sin(theta)**2]
])

# Compute inverse metric g^{mu nu}
g_inv = g.inv()

# Compute Christoffel symbols Gamma^{lambda}_{mu nu}
n = len(coords)
christoffel = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]

for lam in range(n):
    for mu in range(n):
        for nu in range(n):
            # Gamma^{lambda}_{mu nu} = (1/2) g^{lambda sigma} (d_mu g_{nu sigma} + d_nu g_{mu sigma} - d_sigma g_{mu nu})
            for sigma in range(n):
                term1 = diff(g[nu, sigma], coords[mu])
                term2 = diff(g[mu, sigma], coords[nu])
                term3 = diff(g[mu, nu], coords[sigma])
                christoffel[lam][mu][nu] += (sp.Rational(1, 2) * g_inv[lam, sigma] *
                                             (term1 + term2 - term3))

# Compute Riemann tensor R^{rho}_{sigma mu nu}
riemann = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]

for rho in range(n):
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                # R^{rho}_{sigma mu nu} = d_mu Gamma^{rho}_{nu sigma} - d_nu Gamma^{rho}_{mu sigma}
                # + Gamma^{rho}_{mu lambda} Gamma^{lambda}_{nu sigma} - Gamma^{rho}_{nu lambda} Gamma^{lambda}_{mu sigma}
                term1 = diff(christoffel[rho][nu][sigma], coords[mu])
                term2 = diff(christoffel[rho][mu][sigma], coords[nu])

                term3 = 0
                term4 = 0
                for lam in range(n):
                    term3 += christoffel[rho][mu][lam] * christoffel[lam][nu][sigma]
                    term4 += christoffel[rho][nu][lam] * christoffel[lam][mu][sigma]

                riemann[rho][sigma][mu][nu] = term1 - term2 + term3 - term4

# Compute Ricci tensor R_{mu nu} = R^{lambda}_{mu lambda nu}
ricci = zeros(n, n)

for mu in range(n):
    for nu in range(n):
        for lam in range(n):
            ricci[mu, nu] += riemann[lam][mu][lam][nu]

# Compute Ricci scalar R = g^{mu nu} R_{mu nu}
ricci_scalar = 0
for mu in range(n):
    for nu in range(n):
        ricci_scalar += g_inv[mu, nu] * ricci[mu, nu]

# Compute Einstein tensor G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R
einstein = ricci - sp.Rational(1, 2) * g * ricci_scalar

# Convert to mixed Einstein tensor G^{mu}_{nu} = g^{mu alpha} G_{alpha nu}
einstein_mixed = zeros(n, n)
for mu in range(n):
    for nu in range(n):
        for alpha in range(n):
            einstein_mixed[mu, nu] += g_inv[mu, alpha] * einstein[alpha, nu]

# Extract components (t=0, r=1, theta=2, phi=3)
G_tt_computed = einstein_mixed[0, 0]
G_rr_computed = einstein_mixed[1, 1]
G_theta_theta_computed = einstein_mixed[2, 2]
G_phi_phi_computed = einstein_mixed[3, 3]

# Paper's expressions for comparison
f_prime = diff(f, r)  # f' = -2*r*a*exp(-r^2)
G_tt_paper = (1 - f + 2*r**2*f) / r**2
G_rr_paper = (1 - f) / r**2
G_theta_theta_paper = -f
G_phi_phi_paper = -f

# Verify expressions match
verification = {
    "G^t_t matches paper": simplify(G_tt_computed - G_tt_paper) == 0,
    "G^r_r matches paper": simplify(G_rr_computed - G_rr_paper) == 0,
    "G^theta_theta matches paper": simplify(G_theta_theta_computed - G_theta_theta_paper) == 0,
    "G^phi_phi matches paper": simplify(G_phi_phi_computed - G_phi_phi_paper) == 0
}

# Compute effective energy density and pressure
rho_computed = G_tt_computed / (8 * pi * G)
p_computed = -G_rr_computed / (8 * pi * G)

rho_paper = G_tt_paper / (8 * pi * G)
p_paper = -G_rr_paper / (8 * pi * G)

# Compute repulsion criterion
repulsion_computed = rho_computed + 3 * p_computed
repulsion_paper = rho_paper + 3 * p_paper

# Simplified repulsion criterion (Eq. 18)
simplified_repulsion = (a * (1 + r**2) * exp(-r**2) - 1) / (4 * pi * G * r**2)

# Verify repulsion criterion matches
repulsion_match = simplify(repulsion_computed - repulsion_paper) == 0
simplified_match = simplify(repulsion_computed - simplified_repulsion) == 0

# Analyze repulsion criterion sign
factor = a * exp(-r**2) * (1 + r**2) - 1

# Check if factor < 0 for r>0 and Nv>0
h_r = (1 + r**2) * exp(-r**2)
h_max = sp.limit(h_r, r, 0)  # = 1
h_derivative = diff(h_r, r)
critical_points = sp.solve(h_derivative, r)
h_at_critical = [h_r.subs(r, cp) for cp in critical_points if cp.is_real and cp > 0]

print("=== Verification of Computations ===")
print(f"Ricci scalar R = {simplify(ricci_scalar)}")
print(f"Repulsion criterion rho+3p = {simplify(repulsion_computed)}")
print(f"Factor in repulsion: {simplify(factor)}")
print(f"Maximum of h(r) = (1+r^2)exp(-r^2) at r=0: {h_max}")
print(f"h(r) at critical points (r>0): {h_at_critical}")
print("\n=== Verification Against Paper ===")
for key, value in verification.items():
    print(f"{key}: {value}")
print(f"Repulsion criterion matches paper: {repulsion_match}")
print(f"Simplified repulsion matches: {simplified_match}")

# Check sign for specific values
print("\n=== Sign Check for Specific Values ===")
test_values = [
    (0.1, 0.1),   # Small r, small Nv
    (1.0, 0.2),   # Medium r, medium Nv
    (2.0, 0.3),   # Larger r, larger Nv
    (0.5, 0.5)    # Medium r, large Nv
]

for r_val, nv_val in test_values:
    a_val = (1 - nv_val)**(2/3)
    factor_val = a_val * (1 + r_val**2) * exp(-r_val**2) - 1
    repulsion_val = factor_val / (4 * 3.1416 * 1 * r_val**2)  # G=1
    print(f"r={r_val}, Nv={nv_val}: rho+3p = {repulsion_val:.6f} (should be <0)")

# Check behavior at r->0 and r->inf
print("\n=== Asymptotic Behavior ===")
print("As r->0+:")
print(f"  rho ~ {simplify(sp.series(rho_computed, r, 0, 3).removeO())}")
print(f"  p ~ {simplify(sp.series(p_computed, r, 0, 3).removeO())}")
print(f"  rho+3p ~ {simplify(sp.series(repulsion_computed, r, 0, 3).removeO())}")

print("\nAs r->inf:")
print(f"  rho ~ {simplify(sp.series(rho_computed, r, sp.oo, 3).removeO())}")
print(f"  p ~ {simplify(sp.series(p_computed, r, sp.oo, 3).removeO())}")
print(f"  rho+3p ~ {simplify(sp.series(repulsion_computed, r, sp.oo, 3).removeO())}")

# Verify the paper's key result: rho+3p < 0 for all r>0 when Nv>0
print("\n=== Verification of Paper's Key Result ===")
print("Key result: rho+3p < 0 for all r>0 when Nv>0")
print("Proof:")
print("1. a = (1-Nv)^(2/3) < 1 when Nv>0")
print("2. h(r) = (1+r^2)exp(-r^2) <= 1 for all r, with equality only at r=0")
print("3. Therefore a*h(r) < 1 for all r>0")
print("4. Factor = a*h(r) - 1 < 0 for all r>0")
print("5. Denominator 4*pi*G*r^2 > 0 for all r>0")
print("6. Therefore rho+3p = factor/(4*pi*G*r^2) < 0 for all r>0")
