#!/usr/bin/env python3
"""
Wilson RG 2-loop: weryfikacja lambda_Koide vs lambda_WF
========================================================
Rozszerzenie wilson_rg_vmod.py o poprawki 2-petlowe.

Cel: sprawdzic czy lambda_WF(2-loop) ~ lambda_Koide = 5.47e-6,
co zamknaloby O-K1 (Koide z Wilson RG substratu).

Funkcje beta 2-loop dla O(N=1) w d=4-eps:
  beta(u4) = -eps*u4 + (3/(16pi^2))*u4^2 - (17/(3*(16pi^2)^2))*u4^3
  beta(u6) = -(3-d)*u6 + (15/(16pi^2))*u4*u6 + (10/(16pi^2))*u4^2
             - 2-loop corrections to u6

Refs: Pelissetto & Vicari, Phys. Rep. 368 (2002);
      Kleinert & Schulte-Frohlinde, Critical Properties of phi^4 theories

TGP v1 -- 2026-03-31
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

eps = 1.0  # d = 3
S = 1.0 / (8.0 * np.pi**2)  # 1/(8pi^2) = 0.01267

# ============================================================
# 1. Funkcje beta: 1-loop i 2-loop
# ============================================================
# 1-LOOP (standard):
# beta_4^(1) = -eps*u4 + 3*S*u4^2
# beta_6^(1) = -2*eps*u6 + 15*S*u4*u6 + 10*S*u4^2

# 2-LOOP (O(1) Ising universality class):
# beta_4^(2) = -eps*u4 + 3*S*u4^2 - 17/3*S^2*u4^3
# beta_6^(2) = -2*eps*u6 + 15*S*u4*u6 + 10*S*u4^2
#              - (89/3)*S^2*u4^2*u6 - (272/9)*S^2*u4^3
#
# Refs: Kleinert-Schulte-Frohlinde eq.(15.4), (25.17) for O(N=1)
# 2-loop coefficient for u4: -(17/3)*u4^3 term from sunset diagram
# 2-loop coefficient for u6: corrections from overlapping diagrams

# Precyzyjne wspolczynniki 2-loop (O(1), d=4-eps):
# beta(g) = -eps*g + 3*g^2 - 17/3*g^3  (g = S*u4, normalized)
# Punkt staly: g* = eps/3 + 17*eps^2/81 + O(eps^3)
# W naszych zmiennych u4: u4* = g*/S

print("=" * 65)
print("  WILSON RG 2-LOOP DLA SUBSTRATU TGP")
print("=" * 65)

# --- 1-loop ---
u4_1loop = eps / (3 * S)
u6_num_1 = -10 * S * u4_1loop**2
u6_den_1 = -2 * eps + 15 * S * u4_1loop
u6_1loop = u6_num_1 / u6_den_1

print(f"\n--- 1-LOOP ---")
print(f"  u4* = {u4_1loop:.4f}")
print(f"  u6* = {u6_1loop:.4f}")

# --- 2-loop: punkt staly u4 ---
# Przy eps=1 proste rozwiazanie 2-loop nie istnieje (dyskryminant < 0).
# To znany problem rozwinienia epsilon. Standardowe rozwiazanie:
# Pade resummation lub dane z conformal bootstrap / MC.
#
# Rozwinienie epsilon dla g* = u4* * S:
#   g* = eps/3 + 17*eps^2/81 + O(eps^3)    [1-loop + 2-loop]
# Pade[1,1]:
#   g*_Pade = (eps/3) / (1 - 17*eps/(27))
#           = (1/3) / (1 - 17/27) = (1/3) / (10/27) = 9/10 = 0.900
# To jest resumowany punkt staly!
#
# Porownanie z conformal bootstrap (Kos, Poland, Simmons-Duffin 2016):
#   g* = 1.1064 (w konwencji g = S*u4*16pi^2 = u4/(2*pi^2))
#   -> u4* w naszej konwencji = g* / S = ... zalezy od normalizacji

# Metoda Pade[1,1] dla u4*:
g_star_1L = eps / 3        # = 0.333
g_star_2L_coeff = 17.0 / 81  # wspolczynnik eps^2
# Pade[1,1]: g*(eps) = (a1*eps) / (1 - a2*eps)
# g* = a1*eps + a1*a2*eps^2 + ...
# Matching: a1 = 1/3, a1*a2 = 17/81 -> a2 = (17/81)/(1/3) = 17/27
a1_pade = 1.0 / 3.0
a2_pade = 17.0 / 27.0
g_star_pade = a1_pade * eps / (1 - a2_pade * eps)

u4_2loop = g_star_pade / S  # przeksztalcenie na u4

# Dodatkowe: dane z Monte Carlo 3D Isinga (Pelissetto & Vicari 2002)
# Resumowany u4* w naszej konwencji:
# g* (resumed, Borel-Pade) ~ 1.411 (Guida & Zinn-Justin 1998)
g_star_resumed = 1.411  # Borel-Pade piec petli
u4_resumed = g_star_resumed / S

print(f"\n--- 2-LOOP (PADE + RESUMMACJA) ---")
print(f"  g*(1-loop) = eps/3 = {g_star_1L:.4f}")
print(f"  g*(Pade[1,1]) = {g_star_pade:.4f}")
print(f"  g*(Borel-Pade 5L, GZJ98) = {g_star_resumed:.4f}")
print(f"  u4*(1-loop) = {u4_1loop:.4f}")
print(f"  u4*(Pade) = {u4_2loop:.4f}")
print(f"  u4*(resumed) = {u4_resumed:.4f}")
print(f"  Stosunek Pade/1L = {u4_2loop/u4_1loop:.4f}")
print(f"  Stosunek resumed/1L = {u4_resumed/u4_1loop:.4f}")

# --- 2-loop: punkt staly u6 ---
# beta_6 = -2*eps*u6 + 15*S*u4*u6 + 10*S*u4^2
#           - (89/3)*S^2*u4^2*u6 - (272/9)*S^2*u4^3
# Na punkcie stalym u4 = u4*(2L):
# u6*[-2*eps + 15*S*u4* - (89/3)*S^2*u4*^2] = -10*S*u4*^2 + (272/9)*S^2*u4*^3

# u6 na punkcie stalym u4 = u4*(Pade):
u4s = u4_2loop
coeff_u6 = -2*eps + 15*S*u4s - (89.0/3.0)*S**2*u4s**2
rhs_u6 = -10*S*u4s**2 + (272.0/9.0)*S**2*u4s**3
u6_pade = rhs_u6 / coeff_u6 if abs(coeff_u6) > 1e-10 else np.nan

# u6 na punkcie stalym u4 = u4*(resumed):
u4r = u4_resumed
coeff_u6_r = -2*eps + 15*S*u4r - (89.0/3.0)*S**2*u4r**2
rhs_u6_r = -10*S*u4r**2 + (272.0/9.0)*S**2*u4r**3
u6_resumed = rhs_u6_r / coeff_u6_r if abs(coeff_u6_r) > 1e-10 else np.nan

# Literaturowa wartosc: Pelissetto & Vicari (2002), Table 3
# r6 = u6*/(u4*)^2 ~ 1.649 +/- 0.024 (3D Ising, MC)
r6_PV = 1.649  # uniwersalny stosunek (w pewnej konwencji)

u6_2loop = u6_pade  # uzyjemy Pade

print(f"  u6*(Pade) = {u6_pade:.4f}")
print(f"  u6*(resumed) = {u6_resumed:.4f}")
print(f"  Stosunek u6(Pade)/u6(1L) = {u6_pade/u6_1loop:.4f}" if not np.isnan(u6_pade) else "  u6*(Pade) = NaN")
print(f"  r6 = u6*/u4*^2 (Pade) = {u6_pade/u4s**2:.6f}" if not np.isnan(u6_pade) else "")
print(f"  r6 = u6*/u4*^2 (resum) = {u6_resumed/u4r**2:.6f}" if not np.isnan(u6_resumed) else "")
print(f"  r6 (Pelissetto-Vicari 2002 MC) = {r6_PV:.3f}")

# ============================================================
# 2. Mapowanie na lambda_eff TGP
# ============================================================
Phi0 = 24.66
a_Gamma = 0.040049

# lambda_eff ~ |u6*| * normalizacja^2 / Phi0^2
# Normalizacja: a_Gamma^2 (skala substratowa do Plancka)
# Precyzyjniej: lambda = C_norm * a_Gamma^2 / Phi0^2
# gdzie C_norm zawiera stosunek u6*/u4*^{3/2} i czynniki geometryczne

# 1-loop: C_1 = |u6*_1L / u4*_1L^{3/2}| * factor
# 2-loop: C_2 = |u6*_2L / u4*_2L^{3/2}| * factor
# Stosunek C_2/C_1 daje korekre do lambda

ratio_u6u4_1L = abs(u6_1loop) / u4_1loop**1.5
ratio_u6u4_2L = abs(u6_2loop) / u4_2loop**1.5

print(f"\n--- MAPOWANIE NA lambda ---")
print(f"  |u6*/u4*^(3/2)| (1-loop) = {ratio_u6u4_1L:.6f}")
print(f"  |u6*/u4*^(3/2)| (2-loop) = {ratio_u6u4_2L:.6f}")
print(f"  Stosunek 2L/1L = {ratio_u6u4_2L/ratio_u6u4_1L:.4f}")

lambda_1L = a_Gamma**2 / Phi0**2  # = 2.637e-6
lambda_2L = lambda_1L * (ratio_u6u4_2L / ratio_u6u4_1L)

print(f"\n  lambda_WF (1-loop) = {lambda_1L:.4e}")
print(f"  lambda_WF (2-loop) = {lambda_2L:.4e}")

# ============================================================
# 3. Porownanie z lambda_Koide
# ============================================================
lambda_Koide = 5.4709e-6

print(f"\n  lambda_Koide = {lambda_Koide:.4e}")
print(f"\n  Stosunek lambda_Koide / lambda_WF:")
print(f"    1-loop: {lambda_Koide/lambda_1L:.4f}")
print(f"    2-loop: {lambda_Koide/lambda_2L:.4f}")

# ============================================================
# 4. Alternatywne mapowanie: bezposrednie
# ============================================================
# Moze lambda_eff = u6* * (a_Gamma/Phi0)^2 bezposrednio?
# (bo u6* juz jest w jednostkach cutoff)
# lambda ~ |u6*| * (a_Gamma)^2 = ??

lambda_alt_1L = abs(u6_1loop) * a_Gamma**2
lambda_alt_2L = abs(u6_2loop) * a_Gamma**2

print(f"\n  Alternatywne mapowanie lambda = |u6*| * a_G^2:")
print(f"    1-loop: {lambda_alt_1L:.4e}")
print(f"    2-loop: {lambda_alt_2L:.4e}")

# Normalizacja z Phi0
lambda_norm_1L = abs(u6_1loop) * (a_Gamma/Phi0)**2
lambda_norm_2L = abs(u6_2loop) * (a_Gamma/Phi0)**2

print(f"\n  Normalizowane lambda = |u6*| * (a_G/Phi0)^2:")
print(f"    1-loop: {lambda_norm_1L:.4e}")
print(f"    2-loop: {lambda_norm_2L:.4e}")
print(f"    Stosunek do lambda_Koide:")
print(f"      1-loop: {lambda_Koide/lambda_norm_1L:.4f}")
print(f"      2-loop: {lambda_Koide/lambda_norm_2L:.4f}")

# ============================================================
# 5. Scanning: jaka normalizacja C daje lambda_Koide?
# ============================================================
print(f"\n--- INWERSJA: jaki C daje lambda_Koide? ---")
# lambda_Koide = C * a_Gamma^2 / Phi0^2
C_needed = lambda_Koide * Phi0**2 / a_Gamma**2
print(f"  C_needed = lambda_K * Phi0^2 / a_G^2 = {C_needed:.4f}")
print(f"  |u6*(1L)| / u4*(1L)^(3/2) = {ratio_u6u4_1L:.4f}")
print(f"  |u6*(2L)| / u4*(2L)^(3/2) = {ratio_u6u4_2L:.4f}")
print(f"  |u6*(1L)| = {abs(u6_1loop):.4f}")
print(f"  |u6*(2L)| = {abs(u6_2loop):.4f}")

# Jaki czynnik z normalizacji calkowej?
# lambda = |u6*| * g_6 / (5! * 64) * ... ?
# Z dodatekI: lambda_eff = |u6*| * v0^6 / (5! * 64)
# v0 = a_Gamma^{1/3}? v0 = Phi0^{1/2}?

# Proba: lambda = |u6*| / (5! * 64) * (a_Gamma)^k
for k in [0, 1, 2, 3, 4]:
    lam_try = abs(u6_2loop) / (120 * 64) * a_Gamma**k
    ratio_try = lambda_Koide / lam_try if lam_try > 0 else float('inf')
    print(f"  |u6*(2L)|/(5!*64) * a_G^{k} = {lam_try:.4e}  (ratio to target: {ratio_try:.2f})")

# ============================================================
# 6. Numeryczny RG flow 2-loop
# ============================================================
def beta_4_2L(u4, u6):
    return -eps*u4 + 3*S*u4**2 - (17.0/3.0)*S**2*u4**3

def beta_6_2L(u4, u6):
    return (-2*eps*u6 + 15*S*u4*u6 + 10*S*u4**2
            - (89.0/3.0)*S**2*u4**2*u6 - (272.0/9.0)*S**2*u4**3)

def rg_2L(t, y):
    return [beta_4_2L(y[0], y[1]), beta_6_2L(y[0], y[1])]

t_span = (0, 20)
t_eval = np.linspace(0, 20, 500)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for u4b, style in [(5.0, '-'), (20.0, '--'), (50.0, ':')]:
    sol1 = solve_ivp(lambda t,y: [
        -eps*y[0]+3*S*y[0]**2,
        -2*eps*y[1]+15*S*y[0]*y[1]+10*S*y[0]**2
    ], t_span, [u4b, 0], t_eval=t_eval, rtol=1e-10, atol=1e-12)

    sol2 = solve_ivp(rg_2L, t_span, [u4b, 0], t_eval=t_eval,
                     rtol=1e-10, atol=1e-12)

    axes[0].plot(sol1.t, sol1.y[0], f'b{style}', alpha=0.5)
    axes[0].plot(sol2.t, sol2.y[0], f'r{style}', alpha=0.8)
    axes[1].plot(sol1.t, sol1.y[1], f'b{style}', alpha=0.5,
                 label=f'1L u4b={u4b}' if style=='-' else None)
    axes[1].plot(sol2.t, sol2.y[1], f'r{style}', alpha=0.8,
                 label=f'2L u4b={u4b}' if style=='-' else None)

axes[0].axhline(u4_1loop, color='blue', ls='-.', alpha=0.5, label=f'u4* 1L={u4_1loop:.1f}')
axes[0].axhline(u4_2loop, color='red', ls='-.', alpha=0.5, label=f'u4* 2L={u4_2loop:.1f}')
axes[0].set_xlabel('t')
axes[0].set_ylabel('u4')
axes[0].set_title('u4 flow: 1-loop (blue) vs 2-loop (red)')
axes[0].legend(fontsize=8)

axes[1].axhline(u6_1loop, color='blue', ls='-.', alpha=0.5, label=f'u6* 1L={u6_1loop:.1f}')
axes[1].axhline(u6_2loop, color='red', ls='-.', alpha=0.5, label=f'u6* 2L={u6_2loop:.1f}')
axes[1].set_xlabel('t')
axes[1].set_ylabel('u6')
axes[1].set_title('u6 flow: 1-loop (blue) vs 2-loop (red)')
axes[1].legend(fontsize=8)

plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'wilson_rg_2loop.png'), dpi=150)
print(f"\n  Wykres: scripts/wilson_rg_2loop.png")

# ============================================================
# 7. Podsumowanie
# ============================================================
print(f"\n{'='*65}")
print(f"  PODSUMOWANIE 2-LOOP")
print(f"{'='*65}")
print(f"  u4*(1L) = {u4_1loop:.4f}  ->  u4*(2L) = {u4_2loop:.4f}  (zmiana {(u4_2loop/u4_1loop-1)*100:+.1f}%)")
print(f"  u6*(1L) = {u6_1loop:.4f}  ->  u6*(2L) = {u6_2loop:.4f}  (zmiana {(u6_2loop/u6_1loop-1)*100:+.1f}%)")
print(f"  |u6/u4^1.5| ratio 2L/1L = {ratio_u6u4_2L/ratio_u6u4_1L:.4f}")
print(f"")
print(f"  lambda_WF:  1L = {lambda_1L:.4e}  2L = {lambda_2L:.4e}")
print(f"  lambda_Koide      = {lambda_Koide:.4e}")
print(f"  Gap 1L: factor {lambda_Koide/lambda_1L:.2f}")
print(f"  Gap 2L: factor {lambda_Koide/lambda_2L:.2f}")
gap_2L = lambda_Koide/lambda_2L
if abs(gap_2L - 1) < 0.3:
    print(f"  >>> 2-loop ZAMYKA luke O-K1! (gap < 30%)")
elif abs(gap_2L - 1) < 0.5:
    print(f"  >>> 2-loop ZMNIEJSZA luke (gap {abs(gap_2L-1)*100:.0f}%)")
else:
    print(f"  >>> 2-loop NIE zamyka luki (gap {abs(gap_2L-1)*100:.0f}%)")
    print(f"  >>> Potrzebne: 3-loop, Pade resummation, lub MC")

print(f"\nGOTOWE.")
