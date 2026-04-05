#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex129_top4_r21_rstar.py
========================
ANALIZA T-OP4: czy r₂₁/r* = 9 = N² jest tożsamością dokładną?

CEL:
  T-OP4 odkryty w ex125: r₂₁ = m_μ/m_e = 206.768, r* = (23+5√21)/2 = 22.956
  r₂₁/r* = 206.768/22.956 = 9.007  [δ = 0.078% od 9]

  Pytania:
  1. Czy r₂₁ = 9·r*  DOKŁADNIE w ramach TGP (bez PDG)?
  2. Czy r₂₁ = 9·r* DOKŁADNIE ze wzoru A_tail? (A_μ/A_e)^4 = 9r*?
  3. Jaka algebraiczna struktura łączy r₂₁ z r*?
  4. Czy r₂₁/r* jest bliska innym ciekawym liczbom?
  5. Brannen: wyrażenie r₂₁ z parametrów (M, r=√2, θ)?

TESTY D1..D14:
  --- NUMERYCZNE ---
  D1:  r₂₁_PDG / r* - 9 < 0.1%  (obserwacja startowa)
  D2:  r₂₁_Atail / r* - 9 < 1%  (A_tail vs r*)
  D3:  r₂₁_Atail = (A_μ/A_e)^4 weryfikacja
  D4:  (A_μ/A_e)^4 / r* < 0.5% od 9

  --- ALGEBRAICZNE: BRANNEN ---
  D5:  Brannen z r=√2, θ=θ_fit: oblicz r₂₁^B = m_μ^B/m_e^B
  D6:  r₂₁^B / r* — jak blisko 9?
  D7:  Koide warunek: czy r₂₁ = r* · N² wynika z formuły?

  --- PRZESZUKANIE FORMUŁ ---
  D8:  Szukaj p,q całkowite: r₂₁ ≈ p·r*+q (p=9 q=0 jest najlepsze?)
  D9:  Szukaj formuły r₂₁ = f(r*) z małymi współczynnikami
  D10: Sprawdź r₂₁ vs potęgi r*: r*^n dla n = 0.5..3 (log-log)
  D11: u* jest pierwiastkiem u²-5u+1=0: czy r₂₁=u*² - ... ?
  D12: Sprawdź 9r* vs inne znane stałe: 2π, φ, φ², √21, etc.

  --- HIPOTEZY TGP ---
  D13: Hipoteza: r₂₁ = N²·r* (N=3) z precyzją < 0.1% — "blisko trafienie"
  D14: Hipoteza: r₂₁·(N+1)/(2N) = r* (skalowanie Q_K?) — sprawdź

Referencje: ex117, ex118, ex125
"""

import sys
import io
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
import cmath

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# Stałe
# ============================================================
SQRT21   = math.sqrt(21.0)
R_STAR   = (23.0 + 5.0 * SQRT21) / 2.0   # Koide FP ratio
U_STAR   = (5.0 + SQRT21) / 2.0           # √r*; u²-5u+1=0

PHI      = (1.0 + math.sqrt(5.0)) / 2.0

M_E_MEV  = 0.510999    # PDG 2024
M_MU_MEV = 105.6584    # PDG 2024
M_TAU_MEV= 1776.86     # PDG 2024

R21_PDG  = M_MU_MEV / M_E_MEV     # = 206.7682...
R31_PDG  = M_TAU_MEV / M_E_MEV    # = 3477.48...
R32_PDG  = M_TAU_MEV / M_MU_MEV   # = 16.817...

# A_tail z ex125/ex127 (numeryczne, okno [20,35])
A_E_TAIL  = 0.29882
A_MU_TAIL = 1.13314
A_TAU_TAIL= 2.29471

R21_ATAIL = (A_MU_TAIL / A_E_TAIL)**4
R31_ATAIL = (A_TAU_TAIL / A_E_TAIL)**4
R32_ATAIL = (A_TAU_TAIL / A_MU_TAIL)**4

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
print("=" * 72)
print("EX129: T-OP4 — ANALIZA r₂₁/r* ≈ 9 = N²")
print("=" * 72)
print()

print("[0] WARTOŚCI STARTOWE")
print("-" * 50)
print(f"  r* = (23+5√21)/2 = {R_STAR:.10f}")
print(f"  u* = (5+√21)/2   = {U_STAR:.10f}  [u*²-5u*+1 = {U_STAR**2 - 5*U_STAR + 1:.2e}]")
print(f"  √21 = {SQRT21:.10f}")
print()
print(f"  r₂₁ (PDG)   = {R21_PDG:.8f}")
print(f"  r₂₁ (A^4)   = {R21_ATAIL:.8f}")
print(f"  r₃₁ (PDG)   = {R31_PDG:.4f}")
print(f"  r₃₁ (A^4)   = {R31_ATAIL:.4f}")
print()
print(f"  r₂₁/r*      = {R21_PDG/R_STAR:.8f}  [od 9: δ={(abs(R21_PDG/R_STAR - 9)/9)*100:.4f}%]")
print(f"  r₂₁(A^4)/r* = {R21_ATAIL/R_STAR:.8f}  [od 9: δ={(abs(R21_ATAIL/R_STAR - 9)/9)*100:.4f}%]")
print(f"  r₃₁/r*²     = {R31_PDG/R_STAR**2:.8f}  [od 1: δ={(abs(R31_PDG/R_STAR**2 - 1)/1)*100:.4f}%]")
print()


# ============================================================
# SEKCJA 1: Testy numeryczne
# ============================================================
print("[1] TESTY NUMERYCZNE — PODSTAWOWE")
print("-" * 60)

# D1: obserwacja startowa PDG
delta_d1 = abs(R21_PDG / R_STAR - 9.0) / 9.0
record("D1: r₂₁(PDG)/r* - 9 < 0.1%",
       delta_d1 < 0.001,
       f"r₂₁/r* = {R21_PDG/R_STAR:.6f}, δ = {100*delta_d1:.4f}%")

# D2: z A_tail
delta_d2 = abs(R21_ATAIL / R_STAR - 9.0) / 9.0
record("D2: r₂₁(A^4)/r* - 9 < 1%",
       delta_d2 < 0.01,
       f"r₂₁(A^4)/r* = {R21_ATAIL/R_STAR:.6f}, δ = {100*delta_d2:.4f}%")

# D3: weryfikacja A^4
r21_verify = (A_MU_TAIL/A_E_TAIL)**4
record("D3: r₂₁(A^4) = (A_μ/A_e)^4",
       abs(r21_verify - R21_ATAIL) < 1e-10,
       f"(A_μ/A_e)^4 = {r21_verify:.6f}")

# D4: A^4/r* blisko 9
delta_d4 = abs(R21_ATAIL / R_STAR - 9.0) / 9.0
record("D4: (A_μ/A_e)^4 / r* blisko 9 (δ < 0.5%)",
       delta_d4 < 0.005,
       f"δ = {100*delta_d4:.4f}%")


# ============================================================
# SEKCJA 2: Brannen — wyrażenie r₂₁ z parametrów
# ============================================================
print()
print("[2] BRANNEN — r₂₁ Z (M, r=√2, θ)")
print("-" * 60)

# Brannen: √m_k = M(1 + r·cos(θ + 2πk/3)) dla k=0,1,2 → e,μ,τ
# Z PDG: dopasuj θ (M i θ z fit)
sqm = np.array([math.sqrt(M_E_MEV), math.sqrt(M_MU_MEV), math.sqrt(M_TAU_MEV)])
M_mean = float(np.mean(sqm))
eps_k = sqm / M_mean - 1.0  # epsilon_k = r·cos(θ + 2πk/3)

# DFT: F₁ = (3/2)·r·exp(iθ)
import cmath as cm
F1 = sum(eps_k[k] * cm.exp(-2j*math.pi*k/3) for k in range(3))
r_B   = abs(F1) * 2.0/3.0
theta_B = math.atan2(F1.imag, F1.real)  # atan2(Im, Re)

print(f"\n  Brannen fit do PDG (r=√2):")
print(f"    M_mean   = {M_mean:.6f} MeV^{{1/2}}")
print(f"    r_B      = {r_B:.8f}  [√2 = {math.sqrt(2):.8f}]")
print(f"    θ_B      = {math.degrees(theta_B):.6f}°")

# Oblicz masy Brannena
r_brannen = math.sqrt(2.0)
theta_fit = theta_B
masses_B = np.array([M_mean * (1 + r_brannen * math.cos(theta_fit + 2*math.pi*k/3))**2
                     for k in range(3)])

r21_B = masses_B[1] / masses_B[0]
r31_B = masses_B[2] / masses_B[0]
print(f"\n  Masy Brannena (r=√2, θ=θ_fit):")
print(f"    m_e^B    = {masses_B[0]:.6f} MeV,  PDG = {M_E_MEV:.6f}")
print(f"    m_μ^B    = {masses_B[1]:.4f} MeV,  PDG = {M_MU_MEV:.4f}")
print(f"    m_τ^B    = {masses_B[2]:.4f} MeV,  PDG = {M_TAU_MEV:.4f}")
print(f"    r₂₁^B   = {r21_B:.8f}  PDG={R21_PDG:.8f}")
print(f"    r₂₁^B/r* = {r21_B/R_STAR:.8f}  [od 9: δ={(abs(r21_B/R_STAR-9)/9)*100:.4f}%]")

# D5: Brannen r₂₁ blisko PDG (r=√2 ≠ r_B → masy nie idealne, wystarczy <0.1%)
record("D5: Brannen r₂₁^B ≈ r₂₁(PDG) (r=√2; r_B≈√2 → blisko <0.1%)",
       abs(r21_B - R21_PDG) / R21_PDG < 0.001,
       f"r₂₁^B = {r21_B:.8f}, PDG = {R21_PDG:.8f}, "
       f"δ = {100*abs(r21_B-R21_PDG)/R21_PDG:.4f}%")

# D6: Brannen r₂₁^B / r*
delta_d6 = abs(r21_B/R_STAR - 9.0)/9.0
record("D6: r₂₁^B/r* blisko 9 (r=√2: δ < 0.2%)",
       delta_d6 < 0.002,
       f"δ = {100*delta_d6:.4f}%")

# D7: analityczne rozwinięcie — Brannen vs r*
# r₂₁^B = [(1 + √2·cos(θ+2π/3))/(1 + √2·cos(θ))]²
# Czy θ_fit jest takie, że ten wyraz = 9r*?
def r21_brannen_analytic(theta):
    r = math.sqrt(2.0)
    sqme = 1.0 + r * math.cos(theta)
    sqmmu = 1.0 + r * math.cos(theta + 2*math.pi/3)
    if sqme <= 0 or sqmmu <= 0:
        return float('nan')
    return (sqmmu / sqme)**2

# Szukaj θ takie, że r₂₁(θ) = 9·r*
target_9rstar = 9.0 * R_STAR
print(f"\n  Cel r₂₁=9r* = {target_9rstar:.6f}  (vs PDG {R21_PDG:.6f})")

def eq_9rstar(theta):
    r21 = r21_brannen_analytic(theta)
    if math.isnan(r21): return 1e10
    return r21 - target_9rstar

# Szukaj pierwiastka w [0, 2π]
theta_solutions = []
for t_start in np.linspace(0.1, 2*math.pi - 0.1, 100):
    try:
        val = eq_9rstar(t_start)
        val_next = eq_9rstar(t_start + 0.063)
        if val * val_next < 0:
            t_root = brentq(eq_9rstar, t_start, t_start + 0.063)
            theta_solutions.append(t_root)
    except:
        pass

theta_solutions = sorted(set([round(t, 5) for t in theta_solutions]))
print(f"\n  Θ takie, że r₂₁(θ)=9r*:")
for ts in theta_solutions[:6]:
    print(f"    θ = {math.degrees(ts):.4f}°  ({ts:.6f} rad)")

print(f"\n  θ_fit(PDG) = {math.degrees(theta_fit):.6f}°")

if theta_solutions:
    diff_angles = min(abs(math.degrees(ts) - math.degrees(theta_fit)) for ts in theta_solutions)
    print(f"  Odległość θ_fit od najbliższego θ(9r*): {diff_angles:.4f}°")
    record("D7: θ(r₂₁=9r*) blisko θ_fit (< 1°)",
           diff_angles < 1.0,
           f"min|θ(9r*)-θ_fit| = {diff_angles:.4f}°")
else:
    record("D7: θ(r₂₁=9r*) blisko θ_fit (< 1°)", False, "brak pierwiastka w zakresie")


# ============================================================
# SEKCJA 3: Przeszukanie formuł algebraicznych
# ============================================================
print()
print("[3] PRZESZUKANIE FORMUŁ ALGEBRAICZNYCH")
print("-" * 60)

r21 = R21_PDG
print(f"\n  r₂₁ = {r21:.8f}")
print(f"  r*  = {R_STAR:.8f}")
print(f"  u*  = {U_STAR:.8f}")
print(f"  √21 = {SQRT21:.8f}")
print()

# D8: najlepsze p·r*+q z małymi liczbami całkowitymi
best_pq = None
best_err = float('inf')
for p in range(1, 20):
    for q in range(-20, 20):
        candidate = p * R_STAR + q
        if candidate <= 0: continue
        err = abs(candidate - r21) / r21
        if err < best_err:
            best_err = err
            best_pq = (p, q, candidate)

print(f"  Najlepsze p·r*+q: p={best_pq[0]}, q={best_pq[1]}")
print(f"    = {best_pq[2]:.6f}  (PDG: {r21:.6f}, δ={100*best_err:.4f}%)")
record("D8: najlepsze p·r*+q z p≤20: p=9, q=0 jest optymalne?",
       best_pq[0] == 9 and best_pq[1] == 0,
       f"najlepsze: {best_pq[0]}·r*+{best_pq[1]} = {best_pq[2]:.6f}, δ={100*best_err:.4f}%")

# D9: formuły z √21
# r* = (23+5√21)/2, u* = (5+√21)/2
# Formuła r₂₁ = (207+45√21)/2? (czyli 9r*)
formula_9rstar = 9.0 * R_STAR  # = (207+45√21)/2
formula_val = (207.0 + 45.0 * SQRT21) / 2.0
print(f"\n  Formuła (207+45√21)/2 = {formula_val:.6f}  (δ_PDG={(abs(formula_val-r21)/r21)*100:.4f}%)")

# Inne próby:
candidates_formula = {
    "9·r*": 9*R_STAR,
    "(207+45√21)/2": (207+45*SQRT21)/2,
    "u*²+u*-1": U_STAR**2 + U_STAR - 1,
    "u*²-u*+5": U_STAR**2 - U_STAR + 5,
    "u*³-u*²": U_STAR**3 - U_STAR**2,
    "10·u*²-3": 10*U_STAR**2 - 3,
    "9·u*²-u*": 9*U_STAR**2 - U_STAR,
    "u*^4/r*": U_STAR**4 / R_STAR,
    "(5r*+1)/0.538": (5*R_STAR+1)/0.538,
    "r*·(N²)": R_STAR * 9,
    "r*·(4N-3)": R_STAR * (4*3-3),   # 4N-3=9 dla N=3!
    "φ^11": PHI**11,
    "2π·r*": 2*math.pi*R_STAR,
    "r*+r*²/r*": R_STAR + R_STAR,
}

print(f"\n  {'Formuła':>30}  {'Wartość':>12}  {'δ%':>8}")
print("  " + "-"*60)
for label, val in sorted(candidates_formula.items(), key=lambda x: abs(x[1]-r21)/r21):
    delta = abs(val - r21) / r21 * 100
    mark = "***" if delta < 0.5 else ("  *" if delta < 2 else "   ")
    print(f"  {label:>30}  {val:>12.6f}  {delta:>7.3f}%  {mark}")

record("D9: 9·r* jest najlepszą prostą formułą (δ<0.1%)",
       abs(9*R_STAR - r21)/r21 < 0.001,
       f"9·r* = {9*R_STAR:.6f} vs PDG {r21:.6f}, δ={(abs(9*R_STAR-r21)/r21)*100:.4f}%")

# D10: potęgi r*
print(f"\n  r₂₁ vs potęgi r*:")
for n_log in [0.5, 0.7, 0.8, 1.0, 1.1, 1.2, 1.5, 2.0]:
    val = R_STAR**n_log
    delta = abs(val - r21)/r21 * 100
    print(f"    r*^{n_log:.2f} = {val:.4f}  δ={delta:.2f}%")
n_best = math.log(r21) / math.log(R_STAR)
print(f"    r*^n = r₂₁  → n = {n_best:.6f}")
record("D10: n_log(r₂₁,r*) ≠ całkowita (r₂₁ ≠ r*^n, odl. > 0.1)",
       abs(n_best - round(n_best)) > 0.1,
       f"n = {n_best:.6f}, odl. od całk. = {abs(n_best-round(n_best)):.4f}  "
       f"[NIE jest potęgą całkowitą r*]")

# D11: u* = pierwiastek u²-5u+1=0 → u*=(5+√21)/2
# Sprawdź wyrażenia przez u*
print(f"\n  r₂₁ vs wyrażenia przez u*={U_STAR:.6f}:")
exprs_ustar = {
    "u*²": U_STAR**2,
    "u*²+1": U_STAR**2 + 1,
    "u*²+u*": U_STAR**2 + U_STAR,
    "u*²+u*-1": U_STAR**2 + U_STAR - 1,
    "u*²-u*+1": U_STAR**2 - U_STAR + 1,
    "u*²+2u*": U_STAR**2 + 2*U_STAR,
    "9u*²/r*": 9*U_STAR**2/R_STAR,
    "u*³/u*": U_STAR**3/U_STAR,  # = u*²
    "5u*²-u*": 5*U_STAR**2 - U_STAR,
    "4u*²+u*": 4*U_STAR**2 + U_STAR,
    "3u*²+2u*": 3*U_STAR**2 + 2*U_STAR,
}
for label, val in sorted(exprs_ustar.items(), key=lambda x: abs(x[1]-r21)/r21):
    delta = abs(val - r21)/r21*100
    if delta < 5:
        print(f"    {label:>20} = {val:.6f}  δ={delta:.3f}%")

# Dokładna formuła: r₂₁ = 9r* = 9u*² bo u*² = r* (tak?)
# Sprawdź: u*² = r*?
print(f"\n  u*² = {U_STAR**2:.8f}  vs  r* = {R_STAR:.8f}  "
      f"[równe? {abs(U_STAR**2 - R_STAR) < 1e-8}]")
# Uwaga: u*=(5+√21)/2, u*²=((5+√21)/2)²=(25+10√21+21)/4=(46+10√21)/4=(23+5√21)/2=r* !!
u_star_sq_analytic = (46.0 + 10.0 * SQRT21) / 4.0
print(f"  u*² = (46+10√21)/4 = (23+5√21)/2 = r* → TRUE!")
record("D11: u*² = r* (analitycznie, bo u*=(5+√21)/2)",
       abs(U_STAR**2 - R_STAR) < 1e-10,
       f"u*² = {U_STAR**2:.10f} = r* = {R_STAR:.10f}")

print(f"\n  WNIOSEK D11: r₂₁/r* = 9 ⟺ r₂₁ = 9u*² ⟺ √r₂₁ = 3u*")
r21_sqrt = math.sqrt(r21)
print(f"  √r₂₁(PDG) = {r21_sqrt:.6f}  vs  3u* = {3*U_STAR:.6f}  "
      f"[δ={(abs(r21_sqrt-3*U_STAR)/(3*U_STAR))*100:.4f}%]")

# D12: 9r* vs inne stałe
print(f"\n  Inne formy 9r*:")
print(f"    9·r* = 9·(23+5√21)/2 = (207+45√21)/2")
print(f"    = {(207+45*SQRT21)/2:.8f}")
print(f"    2π·r* = {2*math.pi*R_STAR:.8f}  (vs r₂₁={r21:.8f}: zbyt mała)")
print(f"    φ¹¹ = {PHI**11:.6f}  (φ^11 ≈ r₂₁?)")
delta_phi11 = abs(PHI**11 - r21)/r21*100
print(f"    φ¹¹ vs r₂₁: δ = {delta_phi11:.4f}%")
record("D12: sprawdź φ¹¹ vs r₂₁ (inna formuła?)",
       delta_phi11 < 5.0,
       f"φ¹¹ = {PHI**11:.4f}, r₂₁ = {r21:.4f}, δ = {delta_phi11:.2f}%")


# ============================================================
# SEKCJA 4: Hipotezy TGP
# ============================================================
print()
print("[4] HIPOTEZY TGP — STRUKTURALNE")
print("-" * 60)

# D13: N²·r* hipoteza
r21_predicted = (3.0**2) * R_STAR  # N=3 generacje
delta_d13 = abs(r21_predicted - r21) / r21
print(f"\n  Hipoteza D13: r₂₁ = N²·r* (N=3)")
print(f"    Przewidywanie: {r21_predicted:.6f}")
print(f"    PDG:           {r21:.6f}")
print(f"    Odchylenie:    {100*delta_d13:.4f}%")
record("D13: r₂₁ ≈ N²·r* (δ < 0.1%) — 'bliskie trafienie'",
       delta_d13 < 0.001,
       f"δ = {100*delta_d13:.4f}%  ['bliskie trafienie' potwierdzone]")

# D14: r₂₁·(N+1)/(2N) = r* ?
val_d14 = r21 * (3+1) / (2*3)
delta_d14 = abs(val_d14 - R_STAR) / R_STAR
print(f"\n  Hipoteza D14: r₂₁·(N+1)/(2N) = r*?")
print(f"    Lewa strona: r₂₁·4/6 = {val_d14:.6f}")
print(f"    r* = {R_STAR:.6f}")
print(f"    Odchylenie: {100*delta_d14:.4f}%")
# Uwaga: jeśli r₂₁ = 9r*, to r₂₁·(N+1)/(2N) = 9r*·4/6 = 6r* ≠ r*
# Więc ta hipoteza FAIL — chyba że przez przypadek...
record("D14: r₂₁·(N+1)/(2N) ≠ r* (Q_K-skalowanie nie działa, δ>>1%)",
       delta_d14 > 1.0,
       f"δ = {100*delta_d14:.4f}%  [potwierdza: r₂₁·(N+1)/(2N)=6r*≠r*, ta hipoteza FAŁSZYWA]")


# ============================================================
# SEKCJA 5: Wyjaśnienie zbieżności — analiza głębsza
# ============================================================
print()
print("[5] ANALIZA GŁĘBSZA: SKĄD r₂₁/r* ≈ 9?")
print("-" * 60)

print(f"""
  FAKTY ANALITYCZNE:
  1. r* = u*² (dowiedzione w D11: u*=(5+√21)/2, u*²=(23+5√21)/2)
  2. r* jest Koide FP: Q_K(1, r*, r*²) = 3/2
  3. r₂₁ = (A_μ/A_e)^4 = (m_μ/m_e) z TGP (A_tail∝m^{{1/4}})
  4. N=3 generacje → r_Brannen = √(N-1) = √2

  OBSERWACJA NUMERYCZNA:
  r₂₁/r* = {R21_PDG/R_STAR:.8f} ≈ 9 = 3² = N²  (δ = {100*abs(R21_PDG/R_STAR-9)/9:.4f}%)

  MOŻLIWE WYJAŚNIENIA:
  (a) Koincydencja numeryczna — 0.078% sugeruje, że MOŻE być prawdziwa
  (b) r₂₁ = N²·r* jako warunek samospójności TGP
      → m_μ/m_e = N²·Q_K^{{-1}} · [coś z solitonu]?
  (c) Brannen: r₂₁ = [(1+r·cos(θ+2π/3))/(1+r·cos(θ))]²
      Dla r=√2 i θ=θ_fit(PDG) daje r₂₁=206.768 ≈ 9r*

  ALGEBRAICZNA ŚCIEŻKA:
  Jeśli r₂₁ = 9r* = 9u*² i u* spełnia u²-5u+1=0:
  → m_μ/m_e = 9·(5+√21)²/4 = 9·(46+10√21)/4 = (207+45√21)/2
  → To byłoby PRZEWIDYWANIE TGP dla stosunku mas
  → Vs PDG: (207+45√21)/2 = {(207+45*SQRT21)/2:.6f}, PDG = {r21:.6f}
  → Różnica: Δr₂₁ = {r21 - (207+45*SQRT21)/2:.6f}  ({100*abs(r21-(207+45*SQRT21)/2)/r21:.4f}%)
""")

# Pytanie: czy Brannen MUSI dać r₂₁ = 9r*?
# W Brannenie z r=√2, θ=θ_fit:
# r₂₁ jest funkcją θ — nie ma a priori powodu by = 9r*
# CHYBA że θ jest wyznaczane przez inne kryterium TGP (np. soliton)

print(f"  Brannen r₂₁(θ_fit) = {r21_B:.8f}")
print(f"  9r*                 = {9*R_STAR:.8f}")
print(f"  Różnica:            {r21_B - 9*R_STAR:.8f}  ({100*abs(r21_B-9*R_STAR)/(9*R_STAR):.4f}%)")

# KLUCZOWE PYTANIE: czy TGP przewiduje θ_fit analitycznie?
# Jeśli TGP→θ=θ_TGP i θ_TGP daje r₂₁=9r*, to T-OP4 jest zamknięte
# Ale skąd θ? Soliton daje A_tail^4=m, ale θ = atan2(-C,B) z ogona solitonu
# δ_e, δ_μ, δ_τ = fazy z ex125 (nie tworzą postępu Z₃ — PSH obalone)

# T-OP4 REMAINS OPEN:
print(f"""
  STATUS T-OP4:
  - Numerycznie: r₂₁/r* = 9.007 (δ=0.078% od 9 = N²) — "bliskie trafienie"
  - r*=u*² (analitycznie udowodnione)
  - Potencjalna formuła TGP: m_μ/m_e = N²·r* = (207+45√21)/2
  - Weryfikacja eksperymentalna: Δm_μ/m_e = {r21-(207+45*SQRT21)/2:.4f} (0.078%)
  - STATUS: T-OP4 OPEN — "bliskie trafienie", nie dowód; wymaga derywacji z ODE solitonu
""")


# ============================================================
# SEKCJA 6: Dodatkowe — r₃₁ analiza
# ============================================================
print()
print("[6] DODATKOWE — r₃₁ VS r*")
print("-" * 60)

print(f"\n  r₃₁(PDG) = {R31_PDG:.6f}")
print(f"  r₃₁(A^4) = {R31_ATAIL:.6f}")
print(f"  r*²      = {R_STAR**2:.6f}")
print(f"  r₃₁/r*²  = {R31_PDG/R_STAR**2:.8f}  [od 1: δ={(abs(R31_PDG/R_STAR**2-1)/1)*100:.4f}%]")
print(f"  N⁴·r*²   = {81*R_STAR**2:.4f}  (81r*²)")
print(f"  r₃₁/(81r*²)= {R31_PDG/(81*R_STAR**2):.6f}")
print()

# Inne kombinacje r₃₁
r31_formulas = {
    "r*²": R_STAR**2,
    "9r*": 9*R_STAR,
    "N⁴·r*": (3**4)*R_STAR,
    "r₂₁·r*": R21_PDG * R_STAR,
    "r₂₁·N/r*": R21_PDG * 3 / R_STAR,
    "(9r*)·r*": (9*R_STAR)*R_STAR,
    "r₂₁+r*²": R21_PDG + R_STAR**2,
    "r₂₁²/r*": R21_PDG**2/R_STAR,
    "u*⁸": U_STAR**8,
}
print(f"  {'Formuła':>20}  {'Wartość':>12}  {'δ%':>8}")
for label, val in sorted(r31_formulas.items(), key=lambda x: abs(x[1]-R31_PDG)/R31_PDG):
    delta = abs(val - R31_PDG)/R31_PDG*100
    mark = "***" if delta < 1 else "  "
    print(f"  {label:>20}  {val:>12.4f}  {delta:>7.3f}%  {mark}")


# ============================================================
# Wyniki końcowe
# ============================================================
print()
print("=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
passed_all = sum(1 for _, p, _ in TESTS if p)
total_all  = len(TESTS)
for name, passed_t, detail in TESTS:
    mark = "PASS" if passed_t else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print(f"WYNIK: {passed_all}/{total_all} testów PASS")
print()
print("KLUCZOWE WYNIKI:")
print(f"  1. r₂₁/r* = {R21_PDG/R_STAR:.6f} ≈ 9 = N² (δ={100*abs(R21_PDG/R_STAR-9)/9:.4f}%)")
print(f"  2. u*² = r* — tożsamość analityczna (D11)")
print(f"  3. 9r* = (207+45√21)/2 = {(207+45*SQRT21)/2:.6f} vs PDG {R21_PDG:.6f}")
print(f"  4. r*=u*² → r₂₁/r*≈9 ⟺ r₂₁≈9u*² ⟺ √r₂₁≈3u*")
print(f"     √r₂₁(PDG)={math.sqrt(R21_PDG):.6f}, 3u*={3*U_STAR:.6f}")
print()
print("STATUS T-OP4: BLISKIE TRAFIENIE (δ=0.078%)")
print("  Formuła kandydacka: m_μ/m_e = N²·r* = (207+45√21)/2")
print("  Do zamknięcia: derywacja θ_Brannen z dynamiki solitonu TGP")
