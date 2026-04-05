"""
ex10_topological_defect.py
============================
Scenariusz C: cialo jako defekt topologiczny pola TGP.

WYNIKAJACE PYTANIE:
  Czy sila zrodlowa C moze byc wyprowadzona z nonlinearnego rownania pola,
  zamiast byc zadana ad hoc?

STRUKTURA MATEMATYCZNA:
  Pole TGP: g(r) = Phi(r)/Phi_0, gdzie g -> 1 daleko od ciala (proznia S1)
             i g -> g_0 < 1 w centrum "ciala" (lokalny defekt w kierunku S0).

  Statyczne sferyczne ODE (f=1 dla uproszczenia, beta=gamma=1):
     g'' + (2/r)g' = g^2*(g-1)

  Warunek brzegowy:
     g(inf) = 1  (daleko: proznia S1)
     g'(0)  = 0  (gladkosc w centrum)
     g(0)   = g0 (glebokosc defektu - jedyny parametr)

  Zachowanie asymptotyczne dla duzych r:
     g(r) ~ 1 - A*exp(-r)/r   (ogon Yukawa, m_sp=1 dla beta=1)

  => C_eff = A    (amplituda ogona = sila zrodlowa defektu)

  Energia defektu:
     E = 4*pi * integral_0^inf r^2 * [1/2*(g')^2 + V_eff(g)] dr
  gdzie V_eff(g) = -(1/3)g^3 + (1/4)g^4  (potencjal z minimum w g=1)

  KLUCZOWE PYTANIE: czy C_eff jest proporcjonalne do E_defect?
  Jesli tak: zasada rownowaznosci (m_inercyjne = m_ciezkie) wynika
  z teorii pola, nie jest aksjomatem.
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# ===========================================================================
# Krok 1: analiza rownania pola
# ===========================================================================
print("=" * 65)
print("KROK 1: Rownanie defektu topologicznego TGP")
print("=" * 65)
print()
print("ODE (sferyczna symetria, beta=gamma=1):")
print("  g'' + (2/r)*g' = g^2*(g-1)")
print()
print("Zachowanie liniowe przy g~1 (delta = 1-g << 1):")
print("  delta'' + (2/r)*delta' = delta")
print("  => (Laplacian - 1)*delta = 0  => delta ~ A*exp(-r)/r  (Yukawa!)")
print()
print("Zachowanie przy r->0 (Taylor):")
print("  g(r) ~ g0 + [g0^2*(g0-1)/6]*r^2")
print("  g''(0) = g0^2*(g0-1)/3")
print()
print("  Dla g0 < 1: g''(0) < 0 => g maleje od g0 w gore (r wzrasta)")
print("  ALE: to jdzie w strone 0, nie w strone 1!")
print()
print("  => Defekt nie istnieje dla g0 in (0,1) z g'(0)=0, g(inf)=1.")
print()
print("WNIOSEK: Cialo jako 'dziura w prozni' (g0<1) NIE daje rozwiazan")
print("         laczacych g0 z g=1 przy tej postaci ODE.")
print()

# ===========================================================================
# Krok 2: Shooting od duzego r wstecz
# ===========================================================================
print("=" * 65)
print("KROK 2: Shooting wsteczne od r=R_max do r=0")
print("=" * 65)
print()
print("Strategia: zamiast strzelac od r=0 (trudne), strzelamy od r=R_max:")
print("  g(R_max) = 1 - A*exp(-R_max)/R_max")
print("  g'(R_max) = -A*(-1/R_max - 1)*exp(-R_max)/R_max  (approx)")
print()
print("Dla roznych amplitud A > 0, integrujemy do r=0 i patrzymy na g(0).")
print()

def rhs(r, y):
    """ODE: g'' + (2/r)*g' = g^2*(g-1). y = [g, g']."""
    g, gp = y
    if r < 1e-10:
        return [gp, g**2 * (g - 1) / 3.0]
    gpp = g**2 * (g - 1) - 2.0 * gp / r
    return [gp, gpp]

R_max = 30.0
results = []

A_vals = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]
print(f"  {'A':>6}  {'g(0)':>10}  {'g_min':>10}  {'monoton.?':>12}")
print("-" * 45)

for A in A_vals:
    # IC at R_max from Yukawa tail: g = 1 - A*exp(-r)/r
    r0 = R_max
    exp_r0 = np.exp(-r0)
    g0_bc  = 1.0 - A * exp_r0 / r0
    gp0_bc = A * exp_r0 * (1.0/r0 + 1.0/r0**2)  # d/dr[-A*exp(-r)/r] = A*exp(-r)*(1/r+1/r^2)

    sol = solve_ivp(rhs, [R_max, 1e-3], [g0_bc, gp0_bc],
                    method='DOP853', max_step=0.05,
                    rtol=1e-10, atol=1e-12, dense_output=True)

    r_arr = sol.t[::-1]  # reverse: from r=0 outward
    g_arr = sol.y[0][::-1]

    g_center = g_arr[0]  # g at r~0
    g_min = g_arr.min()
    monotone = "si" if np.all(np.diff(g_arr) >= -1e-8) else "nie"

    results.append({'A': A, 'g0': g_center, 'g_min': g_min,
                    'r': r_arr, 'g': g_arr})
    print(f"  {A:>6.2f}  {g_center:>10.5f}  {g_min:>10.5f}  {monotone:>12}")

print()

# ===========================================================================
# Krok 3: Interpretacja
# ===========================================================================
print("=" * 65)
print("KROK 3: Interpretacja fizyczna wynikow")
print("=" * 65)
print()
print("Obserwacja: dla kazdego A > 0, istnieje rozwiazanie laczace")
print("ogon Yukawa g~1-A*exp(-r)/r (r->inf) z dobrze okreslona wartoscia g(0).")
print()

# Wyznacz C_eff i g(0) dla kazdego A
print(f"  {'A=C_eff':>8}  {'g(0)':>10}  {'1-g(0)':>10}  {'sens fiz.':>12}")
print("-" * 48)
for res in results:
    A = res['A']
    g0 = res['g0']
    depth = 1.0 - g0
    phys = "OK (g>0)" if g0 > 0 else "PROBLEM (g<0)"
    if g0 > 0 and g0 < 1:
        print(f"  {A:>8.2f}  {g0:>10.5f}  {depth:>10.5f}  {phys:>12}")
    elif g0 >= 1:
        print(f"  {A:>8.2f}  {g0:>10.5f}  {depth:>10.5f}  nadmierne pole")
    else:
        print(f"  {A:>8.2f}  {g0:>10.5f}  {depth:>10.5f}  g<0 (poza S1)")

print()

# ===========================================================================
# Krok 4: Energia defektu i stosunek C/E
# ===========================================================================
print("=" * 65)
print("KROK 4: Energia defektu E i stosunek C/E")
print("=" * 65)
print()
print("V_eff(g) = -(1/3)*g^3 + (1/4)*g^4  (V(0)=0, V(1)=-1/12)")
print("E = 4*pi * integral r^2 * [1/2*(g')^2 + V_eff(g) - V_eff(1)] dr")
print()

def V_eff(g):
    return -(1.0/3.0)*g**3 + (1.0/4.0)*g**4

V_vac = V_eff(1.0)  # -1/12

print(f"  {'A=C':>6}  {'E_defect':>12}  {'C/E':>10}  {'E/C':>10}")
print("-" * 45)

CE_ratios = []
for res in results:
    A = res['A']
    r_arr = res['r']
    g_arr = res['g']

    if len(r_arr) < 5:
        continue

    # Numeryczne g'
    gp_arr = np.gradient(g_arr, r_arr)

    # Calk energetyczna
    integrand = r_arr**2 * (0.5 * gp_arr**2 + V_eff(g_arr) - V_vac)
    # Usun r=0 (singularnosc)
    mask = r_arr > 0.01
    E = 4.0 * np.pi * np.trapezoid(integrand[mask], r_arr[mask])

    g0 = res['g0']
    if g0 > 0 and g0 < 1 and E > 0:
        ratio = A / E
        CE_ratios.append({'A': A, 'E': E, 'ratio': ratio})
        print(f"  {A:>6.2f}  {E:>12.4f}  {ratio:>10.4f}  {1/ratio:>10.4f}")

print()
if CE_ratios:
    ratios = [r['ratio'] for r in CE_ratios]
    print(f"Ratio C/E: min={min(ratios):.4f}, max={max(ratios):.4f}")
    span = (max(ratios) - min(ratios)) / np.mean(ratios) * 100
    print(f"Zmiennosc C/E: {span:.1f}%")
    if span < 5.0:
        print()
        print("=> C/E ~ stala! Zasada rownowaznosci wynika z teorii pola.")
        print("   C = (const) * E_defect  dla wszystkich defektow TGP.")
    else:
        print()
        print("=> C/E NIE jest stala. Zasada rownowaznosci NIE wynika automatycznie.")
        print("   Rozne defekty maja rozny stosunek ladunku do energii.")

# ===========================================================================
# Wykresy
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
for res in results:
    A = res['A']
    r = res['r']
    g = res['g']
    if r[0] < 0.05:
        mask = r < 15
        ax.plot(r[mask], g[mask], lw=1.5, label=f"A={A:.1f}")
ax.axhline(1.0, color='k', ls='--', lw=1, alpha=0.5, label="g=1 (proznia)")
ax.set_xlabel("r", fontsize=12)
ax.set_ylabel("g(r) = Phi(r)/Phi_0", fontsize=12)
ax.set_title("Profile defektow TGP dla roznych amplitud A=C", fontsize=11)
ax.legend(fontsize=9, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 15)
ax.set_ylim(0.5, 1.05)

ax = axes[1]
if CE_ratios:
    A_plot = [r['A'] for r in CE_ratios]
    E_plot = [r['E'] for r in CE_ratios]
    ax.loglog(E_plot, A_plot, 'bo-', ms=8, lw=2, label="dane TGP")
    # fit
    coeffs = np.polyfit(np.log(E_plot), np.log(A_plot), 1)
    E_fit = np.linspace(min(E_plot)*0.8, max(E_plot)*1.2, 50)
    A_fit = np.exp(coeffs[1]) * E_fit**coeffs[0]
    ax.loglog(E_fit, A_fit, 'r--', lw=2,
              label=f"fit: C ~ E^{{{coeffs[0]:.2f}}}")
    ax.set_xlabel("E_defect (energia)", fontsize=12)
    ax.set_ylabel("C = A (sila zrodlowa)", fontsize=12)
    ax.set_title("Zwiazek C vs E dla defektow TGP", fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "ex10_topological_defect.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print()
print(f"Wykres: {out}")
