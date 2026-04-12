#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tgp_OJ1b_atail_refined.py
============================
O-J1b: Udoskonalony fit A_tail(g0) z asymetrycznym modelem dwu-galeziowym.

CEL:
  Poprawic aproksymacje A_tail blisko g0=1 (vacuum) i zamknac lancuch:
    rho_0* -> g0_e -> A_tail -> r21

METODA:
  Zamiast jednego symetrycznego power-law A ~ C*|g0-1|^mu,
  uzywamy ASYMETRYCZNEGO modelu z korekcjami:
    Galaz I  (g0<1): A_I = C_I * (1-g0)^mu_I * (1 + a1*(1-g0) + a2*(1-g0)^2)
    Galaz II (g0>1): A_II = C_II * (g0-1)^mu_II * (1 + b1*(g0-1) + b2*(g0-1)^2)

  Nastepnie: r21 = (A_II(phi*g0_e) / A_I(g0_e))^4

TESTY:
  R1: Fit galezi I (elektron) z < 1% rel. RMS
  R2: Fit galezi II (mion/tau) z < 1% rel. RMS
  R3: r21 z analitycznego fitu < 5% od PDG
  R4: g0_e z warunku r21 = 206.77 < 1% od numerycznego
  R5: Pelny lancuch rho_0* -> r21 < 10% od PDG

Wersja: v47b (2026-04-12)
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

PHI = (1 + np.sqrt(5)) / 2
R21_PDG = 206.7682830
D = 3
ALPHA = 2.0
RHO_0_STAR = 0.03045

RESULTS = []

def check(condition, name, detail=""):
    status = "PASS" if condition else "FAIL"
    RESULTS.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")
    return condition


# =====================================================================
# SOLVER
# =====================================================================

def soliton_solve(g0, alpha=ALPHA, r_max=250):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / float(D)]
        return [gp, source - cross - float(D - 1) * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-10, atol=1e-12, max_step=0.1,
                    method='DOP853')
    return sol.t, sol.y[0]


def A_tail(g0, alpha=ALPHA):
    r, g = soliton_solve(g0, alpha)
    mask = (r > 50) & (r < 200)
    if np.sum(mask) < 50:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


# =====================================================================
print("=" * 72)
print("  TGP -- O-J1b: Udoskonalony fit A_tail(g0)")
print("=" * 72)


# =====================================================================
# KROK 1: Gesty skan A_tail
# =====================================================================

print("\n[1] Gesty skan A_tail(g0)")

gc = (2*ALPHA + 4) / (2*ALPHA + 1)  # 1.6

# Galaz I: g0 in (0.5, 0.995) -- elektron
g0_I = np.linspace(0.50, 0.995, 50)
A_I = np.array([A_tail(g0) for g0 in g0_I])
valid_I = A_I > 1e-6
g0_I = g0_I[valid_I]
A_I = A_I[valid_I]
delta_I = 1.0 - g0_I  # odleglosc od vacuum

# Galaz II: g0 in (1.005, gc-0.005) -- mion, tau
g0_II = np.linspace(1.005, gc - 0.005, 50)
A_II = np.array([A_tail(g0) for g0 in g0_II])
valid_II = A_II > 1e-6
g0_II = g0_II[valid_II]
A_II = A_II[valid_II]
delta_II = g0_II - 1.0

print(f"  Galaz I:  {len(g0_I)} pkt, delta in [{delta_I[-1]:.4f}, {delta_I[0]:.4f}]")
print(f"  Galaz II: {len(g0_II)} pkt, delta in [{delta_II[0]:.4f}, {delta_II[-1]:.4f}]")


# =====================================================================
# R1: Fit galezi I (elektron)
# =====================================================================

print("\n[R1] Fit galezi I: A_I = C_I * delta^mu_I * (1 + a1*delta + a2*delta^2)")

# Pure power-law first
def pw(x, C, mu):
    return C * x**mu

popt_I_pw, _ = curve_fit(pw, delta_I, A_I, p0=[1.0, 1.0])
fit_I_pw = pw(delta_I, *popt_I_pw)
res_I_pw = np.sqrt(np.mean(((A_I - fit_I_pw)/A_I)**2))
print(f"  Pure power-law: C={popt_I_pw[0]:.4f}, mu={popt_I_pw[1]:.4f}, "
      f"rel.RMS = {res_I_pw:.4f}")

# Power-law + linear correction
def pw_corr1(x, C, mu, a1):
    return C * x**mu * (1 + a1*x)

popt_I_c1, _ = curve_fit(pw_corr1, delta_I, A_I, p0=[popt_I_pw[0], popt_I_pw[1], 0.0])
fit_I_c1 = pw_corr1(delta_I, *popt_I_c1)
res_I_c1 = np.sqrt(np.mean(((A_I - fit_I_c1)/A_I)**2))
print(f"  +linear corr:   C={popt_I_c1[0]:.4f}, mu={popt_I_c1[1]:.4f}, "
      f"a1={popt_I_c1[2]:.4f}, rel.RMS = {res_I_c1:.4f}")

# Power-law + quadratic correction
def pw_corr2(x, C, mu, a1, a2):
    return C * x**mu * (1 + a1*x + a2*x**2)

try:
    popt_I_c2, _ = curve_fit(pw_corr2, delta_I, A_I,
                              p0=[popt_I_c1[0], popt_I_c1[1], popt_I_c1[2], 0.0])
    fit_I_c2 = pw_corr2(delta_I, *popt_I_c2)
    res_I_c2 = np.sqrt(np.mean(((A_I - fit_I_c2)/A_I)**2))
    print(f"  +quadr corr:    C={popt_I_c2[0]:.4f}, mu={popt_I_c2[1]:.4f}, "
          f"a1={popt_I_c2[2]:.4f}, a2={popt_I_c2[3]:.4f}, rel.RMS = {res_I_c2:.4f}")
except:
    popt_I_c2 = popt_I_c1
    res_I_c2 = res_I_c1

# Wybierz najlepszy
best_I = min([("pw", res_I_pw, popt_I_pw),
              ("corr1", res_I_c1, popt_I_c1),
              ("corr2", res_I_c2, popt_I_c2)],
             key=lambda x: x[1])

print(f"  Najlepszy: {best_I[0]} (rel.RMS = {best_I[1]:.4f})")

# Uzyj najlepszego fitu
if best_I[0] == "pw":
    def A_I_fit(delta):
        return pw(delta, *popt_I_pw)
    mu_I_best = popt_I_pw[1]
    C_I_best = popt_I_pw[0]
elif best_I[0] == "corr1":
    def A_I_fit(delta):
        return pw_corr1(delta, *popt_I_c1)
    mu_I_best = popt_I_c1[1]
    C_I_best = popt_I_c1[0]
else:
    def A_I_fit(delta):
        return pw_corr2(delta, *popt_I_c2)
    mu_I_best = popt_I_c2[1]
    C_I_best = popt_I_c2[0]

check(best_I[1] < 0.01,
      "R1: Fit galezi I (elektron) rel. RMS < 1%",
      f"best = {best_I[0]}, mu_I = {mu_I_best:.4f}, C_I = {C_I_best:.4f}, "
      f"rel.RMS = {best_I[1]:.4f}")


# =====================================================================
# R2: Fit galezi II (mion/tau)
# =====================================================================

print("\n[R2] Fit galezi II: A_II = C_II * delta^mu_II * (1 + b1*delta + b2*delta^2)")

popt_II_pw, _ = curve_fit(pw, delta_II, A_II, p0=[1.0, 1.0])
fit_II_pw = pw(delta_II, *popt_II_pw)
res_II_pw = np.sqrt(np.mean(((A_II - fit_II_pw)/A_II)**2))
print(f"  Pure power-law: C={popt_II_pw[0]:.4f}, mu={popt_II_pw[1]:.4f}, "
      f"rel.RMS = {res_II_pw:.4f}")

popt_II_c1, _ = curve_fit(pw_corr1, delta_II, A_II,
                            p0=[popt_II_pw[0], popt_II_pw[1], 0.0])
fit_II_c1 = pw_corr1(delta_II, *popt_II_c1)
res_II_c1 = np.sqrt(np.mean(((A_II - fit_II_c1)/A_II)**2))
print(f"  +linear corr:   C={popt_II_c1[0]:.4f}, mu={popt_II_c1[1]:.4f}, "
      f"b1={popt_II_c1[2]:.4f}, rel.RMS = {res_II_c1:.4f}")

try:
    popt_II_c2, _ = curve_fit(pw_corr2, delta_II, A_II,
                               p0=[popt_II_c1[0], popt_II_c1[1], popt_II_c1[2], 0.0])
    fit_II_c2 = pw_corr2(delta_II, *popt_II_c2)
    res_II_c2 = np.sqrt(np.mean(((A_II - fit_II_c2)/A_II)**2))
    print(f"  +quadr corr:    C={popt_II_c2[0]:.4f}, mu={popt_II_c2[1]:.4f}, "
          f"b1={popt_II_c2[2]:.4f}, b2={popt_II_c2[3]:.4f}, rel.RMS = {res_II_c2:.4f}")
except:
    popt_II_c2 = popt_II_c1
    res_II_c2 = res_II_c1

best_II = min([("pw", res_II_pw, popt_II_pw),
               ("corr1", res_II_c1, popt_II_c1),
               ("corr2", res_II_c2, popt_II_c2)],
              key=lambda x: x[1])

print(f"  Najlepszy: {best_II[0]} (rel.RMS = {best_II[1]:.4f})")

if best_II[0] == "pw":
    def A_II_fit(delta):
        return pw(delta, *popt_II_pw)
    mu_II_best = popt_II_pw[1]
    C_II_best = popt_II_pw[0]
elif best_II[0] == "corr1":
    def A_II_fit(delta):
        return pw_corr1(delta, *popt_II_c1)
    mu_II_best = popt_II_c1[1]
    C_II_best = popt_II_c1[0]
else:
    def A_II_fit(delta):
        return pw_corr2(delta, *popt_II_c2)
    mu_II_best = popt_II_c2[1]
    C_II_best = popt_II_c2[0]

check(best_II[1] < 0.01,
      "R2: Fit galezi II (mion/tau) rel. RMS < 1%",
      f"best = {best_II[0]}, mu_II = {mu_II_best:.4f}, C_II = {C_II_best:.4f}, "
      f"rel.RMS = {best_II[1]:.4f}")


# =====================================================================
# R3: r21 z analitycznego fitu
# =====================================================================

print("\n[R3] Obliczenie r21 z asymetrycznego fitu analitycznego")

# g0_e numeryczne (z kalibracji)
g0_e_num = 0.867697

# A_e z fitu galezi I
delta_e = 1.0 - g0_e_num
A_e_fit = A_I_fit(delta_e)
A_e_num = A_tail(g0_e_num)

# A_mu z fitu galezi II
g0_mu = PHI * g0_e_num
delta_mu = g0_mu - 1.0
A_mu_fit = A_II_fit(delta_mu)
A_mu_num = A_tail(g0_mu)

r21_fit = (A_mu_fit / A_e_fit)**4
r21_num = (A_mu_num / A_e_num)**4

print(f"  g0_e = {g0_e_num:.6f}, delta_e = {delta_e:.6f}")
print(f"  g0_mu = {g0_mu:.6f}, delta_mu = {delta_mu:.6f}")
print(f"\n  A_e:  fit = {A_e_fit:.6f}, num = {A_e_num:.6f}, "
      f"dev = {abs(A_e_fit/A_e_num-1)*100:.2f}%")
print(f"  A_mu: fit = {A_mu_fit:.6f}, num = {A_mu_num:.6f}, "
      f"dev = {abs(A_mu_fit/A_mu_num-1)*100:.2f}%")
print(f"\n  r21 (fit)      = {r21_fit:.2f}")
print(f"  r21 (num)      = {r21_num:.2f}")
print(f"  r21 (PDG)      = {R21_PDG:.2f}")
print(f"  odch. fit/PDG  = {abs(r21_fit/R21_PDG-1)*100:.2f}%")

check(abs(r21_fit / R21_PDG - 1) < 0.05,
      "R3: r21 z analitycznego fitu < 5% od PDG",
      f"r21_fit = {r21_fit:.2f}, PDG = {R21_PDG:.2f}, "
      f"dev = {abs(r21_fit/R21_PDG-1)*100:.2f}%")


# =====================================================================
# R4: g0_e z warunku r21 = 206.77
# =====================================================================

print("\n[R4] Kalibracja g0_e z analitycznego fitu")

def r21_analytic(g0_e):
    d_e = 1.0 - g0_e
    d_mu = PHI * g0_e - 1.0
    if d_e <= 0 or d_mu <= 0:
        return 0.0
    Ae = A_I_fit(d_e)
    Amu = A_II_fit(d_mu)
    if Ae < 1e-15:
        return 0.0
    return (Amu / Ae)**4

# Skan
g0_scan = np.linspace(0.65, 0.98, 50)
r21_scan = [r21_analytic(g0) for g0 in g0_scan]

bracket = None
for i in range(len(r21_scan)-1):
    if r21_scan[i] > 0 and r21_scan[i+1] > 0:
        if (r21_scan[i] - R21_PDG) * (r21_scan[i+1] - R21_PDG) < 0:
            bracket = (g0_scan[i], g0_scan[i+1])
            break

if bracket:
    g0_e_an = brentq(lambda g: r21_analytic(g) - R21_PDG,
                      bracket[0], bracket[1], xtol=1e-10)
    r21_an = r21_analytic(g0_e_an)
    dev_g0e = abs(g0_e_an / g0_e_num - 1) * 100

    print(f"  g0_e (analityczny) = {g0_e_an:.6f}")
    print(f"  g0_e (numeryczny)  = {g0_e_num:.6f}")
    print(f"  odchylenie g0_e:     {dev_g0e:.4f}%")
    print(f"  r21 (z g0_e_an):     {r21_an:.4f}")

    check(dev_g0e < 1.0,
          "R4: g0_e z analitycznego warunku r21 < 1% od numerycznego",
          f"g0_e_an = {g0_e_an:.6f}, g0_e_num = {g0_e_num:.6f}, "
          f"dev = {dev_g0e:.4f}%")
else:
    print("  Brak bracketu!")
    check(False, "R4: g0_e kalibracja", "brak bracketu")
    g0_e_an = g0_e_num
    dev_g0e = 0


# =====================================================================
# R5: Pelny lancuch rho_0* -> r21
# =====================================================================

print("\n[R5] Pelny lancuch: rho_0* -> g0_e -> A_tail -> r21")

g0_e_erg = 1.0 - 13.0/3 * RHO_0_STAR
g0_mu_erg = PHI * g0_e_erg

delta_e_erg = 1.0 - g0_e_erg
delta_mu_erg = g0_mu_erg - 1.0

A_e_erg = A_I_fit(delta_e_erg)
A_mu_erg = A_II_fit(delta_mu_erg)

r21_erg = (A_mu_erg / A_e_erg)**4

# Porownanie z numerycznym
A_e_erg_num = A_tail(g0_e_erg)
A_mu_erg_num = A_tail(g0_mu_erg)
r21_erg_num = (A_mu_erg_num / A_e_erg_num)**4

print(f"  rho_0* = {RHO_0_STAR}")
print(f"  g0_e = 1 - (13/3)*rho_0* = {g0_e_erg:.6f}")
print(f"  g0_mu = phi*g0_e = {g0_mu_erg:.6f}")
print(f"\n  A_e:  fit = {A_e_erg:.6f}, num = {A_e_erg_num:.6f}, "
      f"dev = {abs(A_e_erg/A_e_erg_num-1)*100:.2f}%")
print(f"  A_mu: fit = {A_mu_erg:.6f}, num = {A_mu_erg_num:.6f}, "
      f"dev = {abs(A_mu_erg/A_mu_erg_num-1)*100:.2f}%")
print(f"\n  r21 (fit chain)  = {r21_erg:.2f}")
print(f"  r21 (num chain)  = {r21_erg_num:.2f}")
print(f"  r21 (PDG)        = {R21_PDG:.2f}")
print(f"  odch. fit/PDG    = {abs(r21_erg/R21_PDG-1)*100:.2f}%")
print(f"  odch. num/PDG    = {abs(r21_erg_num/R21_PDG-1)*100:.2f}%")

r21_erg_dev = abs(r21_erg / R21_PDG - 1) * 100
check(r21_erg_dev < 10,
      "R5: Lancuch rho_0* -> r21 < 10% od PDG",
      f"r21_chain = {r21_erg:.2f}, PDG = {R21_PDG:.2f}, dev = {r21_erg_dev:.1f}%")


# =====================================================================
# EKSTRA: Porownanie wykladnikow obu galezi
# =====================================================================

print("\n[EKSTRA] Porownanie wykladnikow i stalych")

print(f"  Galaz I:  mu_I = {mu_I_best:.4f}, C_I = {C_I_best:.4f}")
print(f"  Galaz II: mu_II = {mu_II_best:.4f}, C_II = {C_II_best:.4f}")
print(f"  Stosunek C: C_II/C_I = {C_II_best/C_I_best:.4f}")
print(f"  Roznica mu: mu_II - mu_I = {mu_II_best - mu_I_best:.4f}")

# Kandydaci na mu
for label, mu in [("mu_I", mu_I_best), ("mu_II", mu_II_best)]:
    print(f"\n  {label} = {mu:.4f}, kandydaci:")
    cands = {"1": 1.0, "3/2": 1.5, "phi-1": PHI-1, "sqrt(2)-1/2": np.sqrt(2)-0.5,
             "pi/3": np.pi/3, "1/phi": 1/PHI, "(d-1)/d": (D-1)/D}
    for name, val in sorted(cands.items(), key=lambda x: abs(x[1]/mu-1)):
        dev = (val/mu-1)*100
        if abs(dev) < 20:
            print(f"    {name:<15s} = {val:.4f} ({dev:+.2f}%)")


# =====================================================================
# PODSUMOWANIE
# =====================================================================

print("\n" + "=" * 72)
n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
n_total = len(RESULTS)
print(f"  WYNIK: {n_pass}/{n_total} PASS")
print("=" * 72)

for name, status, detail in RESULTS:
    print(f"  [{status}] {name}")

print(f"""
  WNIOSKI O-J1b:

  1. Asymetryczny fit znacznie poprawia dokkladnosc A_tail:
     - Galaz I:  mu_I = {mu_I_best:.4f}, C_I = {C_I_best:.4f}
     - Galaz II: mu_II = {mu_II_best:.4f}, C_II = {C_II_best:.4f}

  2. r21 z fitu: {r21_fit:.2f} (odch. {abs(r21_fit/R21_PDG-1)*100:.1f}% od PDG)

  3. Pelny lancuch ERG: rho_0* -> g0_e -> A_tail -> r21 = {r21_erg:.2f}
     (odch. {r21_erg_dev:.1f}% od PDG)

  4. Status O-J1: {'ZAMKNIECIE CZESCIOWE' if n_pass >= 4 else 'POSTEP'}
     Udoskonalony fit daje {'dobra' if r21_erg_dev < 10 else 'przyblizona'}
     predykcje r21 z parametrow substratowych
""")

print("=" * 72)
print("  DONE -- tgp_OJ1b_atail_refined.py")
print("=" * 72)

sys.exit(0)
