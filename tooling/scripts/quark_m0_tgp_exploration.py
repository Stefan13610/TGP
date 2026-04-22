#!/usr/bin/env python3
"""
TGP R12 Exploration: Masa konfinementowa m_0 z parametrow TGP.

Cel: znalezc formule laczaca m_0 z Phi_0, a_Gamma, N_c, alpha_s, Lambda_QCD
bez wolnych parametrow.

Strategia:
  1. Estymacja m_0 z potencjalu liniowego TGP (rezim III)
  2. Skan kombinacji stalych TGP dajacych m_0^(d) ~ 22 MeV i m_0^(u) ~ 1982 MeV
  3. Poszukiwanie uniwersalnej formuly

Autor: TGP v47 exploration
Data: 2026-04-09
"""

import numpy as np
from itertools import product as iprod

# ============================================================
# Dane PDG (masy MS-bar przy mu=2 GeV)
# ============================================================
quarks = {
    'down': {
        'm': [4.67, 93.4, 4180.0],   # d, s, b [MeV]
        'names': ['d', 's', 'b'],
    },
    'up': {
        'm': [2.16, 1270.0, 172760.0],  # u, c, t [MeV]
        'names': ['u', 'c', 't'],
    },
    'lepton': {
        'm': [0.51099895, 105.6583755, 1776.86],  # e, mu, tau [MeV]
        'names': ['e', 'mu', 'tau'],
    },
}

# Stale TGP
Phi_0 = 24.66
a_Gamma = 0.040
N_c = 3
N_f = 5
alpha_s_MZ = 0.1179
Lambda_QCD = 217.0  # MeV (PDG, N_f=5)
sigma_str = 440.0**2  # MeV^2 (string tension)
f_pi = 130.7  # MeV (pion decay constant)
phi = (1 + np.sqrt(5)) / 2

# ============================================================
# 1. Wyznaczenie m_0 z warunku shifted Koide
# ============================================================
def koide_Q(masses):
    """Koide Q = (sum sqrt(m))^2 / sum(m)"""
    sm = sum(np.sqrt(m) for m in masses)
    return sm**2 / sum(masses)

def find_m0_shifted_koide(masses, target_Q=1.5):
    """Znajdz m_0 takie ze Q(m_k + m_0) = target_Q."""
    from scipy.optimize import brentq
    def f(m0):
        shifted = [m + m0 for m in masses]
        return koide_Q(shifted) - target_Q
    # Sprawdz czy leptony (m_0 = 0)
    if abs(f(0)) < 1e-6:
        return 0.0
    try:
        return brentq(f, 0, 1e6, xtol=1e-10)
    except ValueError:
        return None

print("=" * 70)
print("1. SHIFTED KOIDE: wyznaczenie m_0")
print("=" * 70)

m0_values = {}
for sector, data in quarks.items():
    m0 = find_m0_shifted_koide(data['m'])
    m0_values[sector] = m0
    Q_orig = koide_Q(data['m'])
    Q_shift = koide_Q([m + m0 for m in data['m']]) if m0 is not None else None
    print(f"\n  {sector}: m = {data['m']}")
    print(f"    Q_orig = {Q_orig:.6f} (odch. od 3/2: {(Q_orig/1.5-1)*100:+.3f}%)")
    print(f"    m_0 = {m0:.2f} MeV" if m0 else "    m_0 = N/A")
    if Q_shift:
        print(f"    Q_shifted = {Q_shift:.8f}")

# ============================================================
# 2. Uniwersalna stala A
# ============================================================
print("\n" + "=" * 70)
print("2. UNIWERSALNA STALA A = m_0 * m_1 / m_3")
print("=" * 70)

A_values = {}
for sector, data in quarks.items():
    m0 = m0_values[sector]
    if m0 and m0 > 0:
        A = m0 * data['m'][0] / data['m'][2]
        A_values[sector] = A
        print(f"  A({sector}) = {m0:.2f} * {data['m'][0]} / {data['m'][2]} = {A:.6f}")

if len(A_values) >= 2:
    vals = list(A_values.values())
    A_mean = np.mean(vals)
    A_cv = np.std(vals) / A_mean * 100
    print(f"\n  A_mean = {A_mean:.6f}")
    print(f"  CV(A) = {A_cv:.1f}%")

# ============================================================
# 3. Skan formul TGP dla m_0
# ============================================================
print("\n" + "=" * 70)
print("3. SKAN FORMUL TGP DLA m_0")
print("=" * 70)

m0_d = m0_values.get('down', 21.9)
m0_u = m0_values.get('up', 1981.5)

# Skala QCD w MeV
# m_0^(d) ~ 22 MeV ~ Lambda_QCD / 10
# m_0^(u) ~ 1982 MeV ~ 9 * Lambda_QCD

candidates = {}

# Hipoteza H1: m_0 = sigma / Lambda_QCD * (m_1/Lambda_QCD)^p
# Dla down: 22 = 0.194e6 / 217 * (4.67/217)^p => 894 * 0.0215^p
# ln(22/894) / ln(0.0215) = ln(0.0246) / ln(0.0215) ~ 0.97
# Dla up: 1982 = 894 * (2.16/217)^p => 894 * 0.00995^p
# ln(2.22) / ln(0.00995) = ... nie dziala

# Hipoteza H2: m_0 = c * f_pi^2 / m_1 (chiral condensate)
m0_d_H2 = f_pi**2 / quarks['down']['m'][0]  # = 130.7^2 / 4.67 = 3658
m0_u_H2 = f_pi**2 / quarks['up']['m'][0]    # = 130.7^2 / 2.16 = 7905
candidates['H2: f_pi^2/m_1'] = (m0_d_H2, m0_u_H2)

# Hipoteza H3: m_0 = Lambda_QCD^2 / m_1
m0_d_H3 = Lambda_QCD**2 / quarks['down']['m'][0]
m0_u_H3 = Lambda_QCD**2 / quarks['up']['m'][0]
candidates['H3: Lambda^2/m_1'] = (m0_d_H3, m0_u_H3)

# Hipoteza H4: m_0 = sigma^(1/2) * a_Gamma / (m_1 * Phi_0) * 1e6
# sigma^(1/2) = 440 MeV
m0_d_H4 = 440 * a_Gamma * 1e3 / (quarks['down']['m'][0])
m0_u_H4 = 440 * a_Gamma * 1e3 / (quarks['up']['m'][0])
candidates['H4: sqrt(sigma)*a_G*1e3/m_1'] = (m0_d_H4, m0_u_H4)

# Hipoteza H5: m_0 = A_mean * m_3 / m_1 (z uniwersalnej A)
A_mean_val = np.mean(list(A_values.values())) if A_values else 0.0246
m0_d_H5 = A_mean_val * quarks['down']['m'][2] / quarks['down']['m'][0]
m0_u_H5 = A_mean_val * quarks['up']['m'][2] / quarks['up']['m'][0]
candidates['H5: A*m_3/m_1 (tautol.)'] = (m0_d_H5, m0_u_H5)

# Hipoteza H6: m_0 = N_c * Lambda_QCD * (m_3/m_1)^(-1/2)
m0_d_H6 = N_c * Lambda_QCD * (quarks['down']['m'][2] / quarks['down']['m'][0])**(-0.5)
m0_u_H6 = N_c * Lambda_QCD * (quarks['up']['m'][2] / quarks['up']['m'][0])**(-0.5)
candidates['H6: Nc*Lambda*(m3/m1)^(-1/2)'] = (m0_d_H6, m0_u_H6)

# Hipoteza H7: m_0 = a_Gamma * m_3 / Phi_0
m0_d_H7 = a_Gamma * quarks['down']['m'][2] / Phi_0
m0_u_H7 = a_Gamma * quarks['up']['m'][2] / Phi_0
candidates['H7: a_G*m_3/Phi_0'] = (m0_d_H7, m0_u_H7)

# Hipoteza H8: m_0 proporcjonalne do masy trzeciej generacji
# m_0 / m_3 = const?
ratio_d = m0_d / quarks['down']['m'][2]  # 22/4180 = 0.00526
ratio_u = m0_u / quarks['up']['m'][2]    # 1982/172760 = 0.01147
candidates['H8_info: m_0/m_3'] = (ratio_d * quarks['down']['m'][2],
                                   ratio_u * quarks['up']['m'][2])

# Hipoteza H9: m_0 = sqrt(sigma) * N_c / Phi_0^(1/2)
m0_H9 = 440 * N_c / np.sqrt(Phi_0)
candidates['H9: sqrt(sigma)*Nc/sqrt(Phi_0)'] = (m0_H9, m0_H9)

# Hipoteza H10: m_0 / sqrt(m_1*m_3) = const?
geom_d = np.sqrt(quarks['down']['m'][0] * quarks['down']['m'][2])
geom_u = np.sqrt(quarks['up']['m'][0] * quarks['up']['m'][2])
r_d = m0_d / geom_d
r_u = m0_u / geom_u
candidates['H10_info: m_0/sqrt(m1*m3)'] = (r_d * geom_d, r_u * geom_u)

# Hipoteza H11: m_0 = (alpha_s * m_3)^(1/2) * C
# sqrt(0.1179 * 4180) = sqrt(492.8) = 22.2  !!!
m0_d_H11 = np.sqrt(alpha_s_MZ * quarks['down']['m'][2])
m0_u_H11 = np.sqrt(alpha_s_MZ * quarks['up']['m'][2])
candidates['H11: sqrt(alpha_s * m_3)'] = (m0_d_H11, m0_u_H11)

print(f"\n  Cel: m_0(d) = {m0_d:.1f} MeV,  m_0(u) = {m0_u:.1f} MeV\n")
print(f"  {'Hipoteza':<35} {'m0_d':>10} {'m0_u':>10} {'dev_d':>8} {'dev_u':>8} {'score':>8}")
print("  " + "-" * 85)

best_score = 999
best_name = ""
for name, (vd, vu) in candidates.items():
    dev_d = (vd / m0_d - 1) * 100 if m0_d > 0 else 999
    dev_u = (vu / m0_u - 1) * 100 if m0_u > 0 else 999
    score = np.sqrt(dev_d**2 + dev_u**2)
    if 'info' not in name and score < best_score:
        best_score = score
        best_name = name
    print(f"  {name:<35} {vd:>10.1f} {vu:>10.1f} {dev_d:>+8.1f}% {dev_u:>+8.1f}% {score:>8.1f}")

print(f"\n  >>> Najlepsza: {best_name} (score = {best_score:.1f})")

# ============================================================
# 4. Glebsza analiza H11: m_0 = sqrt(alpha_s * m_3)
# ============================================================
print("\n" + "=" * 70)
print("4. ANALIZA H11: m_0 = sqrt(alpha_s * m_3)")
print("=" * 70)

for sector, data in quarks.items():
    m0 = m0_values.get(sector)
    if m0 is not None and m0 > 0:
        m3 = data['m'][2]
        pred = np.sqrt(alpha_s_MZ * m3)
        ratio = m0 / pred
        print(f"  {sector}: m_0 = {m0:.2f}, sqrt(alpha_s*m_3) = {pred:.2f}, ratio = {ratio:.4f}")

# Sprawdz z running alpha_s
# alpha_s(m_b) ~ 0.224
alpha_s_mb = 0.224
m0_d_run = np.sqrt(alpha_s_mb * quarks['down']['m'][2])
print(f"\n  Z alpha_s(m_b) = {alpha_s_mb}:")
print(f"    m_0(d) = sqrt({alpha_s_mb} * {quarks['down']['m'][2]}) = {m0_d_run:.2f} (cel: {m0_d:.2f})")

# alpha_s(m_t) ~ 0.108
alpha_s_mt = 0.108
m0_u_run = np.sqrt(alpha_s_mt * quarks['up']['m'][2])
print(f"    m_0(u) = sqrt({alpha_s_mt} * {quarks['up']['m'][2]}) = {m0_u_run:.2f} (cel: {m0_u:.2f})")

# ============================================================
# 5. Hipoteza zlozona: m_0 = sqrt(alpha_s(m_3) * m_3)
# ============================================================
print("\n" + "=" * 70)
print("5. HIPOTEZA: m_0 = sqrt(alpha_s(m_3) * m_3) [running]")
print("=" * 70)

# 1-loop running: alpha_s(mu) = alpha_s(M_Z) / (1 + b_0*alpha_s(M_Z)*ln(mu/M_Z)/(2*pi))
M_Z = 91187.6  # MeV
b_0 = (33 - 2*N_f) / (12 * np.pi)

def alpha_s_running(mu_MeV, nf=5):
    """1-loop alpha_s running."""
    b0_loc = (33 - 2*nf) / (12 * np.pi)
    return alpha_s_MZ / (1 + b0_loc * alpha_s_MZ * np.log(mu_MeV / M_Z) / (2*np.pi))

for sector, data in quarks.items():
    m0 = m0_values.get(sector)
    if m0 is not None and m0 > 0:
        m3 = data['m'][2]
        nf_eff = 5 if m3 > 4000 else (4 if m3 > 1000 else 3)
        a_s = alpha_s_running(m3, nf=nf_eff)
        pred = np.sqrt(a_s * m3)
        dev = (pred / m0 - 1) * 100
        print(f"  {sector}: m_3 = {m3:.0f} MeV, alpha_s(m_3) = {a_s:.4f}")
        print(f"    sqrt(alpha_s * m_3) = {pred:.2f} MeV (cel: {m0:.2f}, odch: {dev:+.1f}%)")

# ============================================================
# 6. Sprawdzenie TGP native: m_0 i stosunek phi
# ============================================================
print("\n" + "=" * 70)
print("6. RELACJE m_0 Z GOLDEN RATIO")
print("=" * 70)

# m_0^(u) / m_0^(d)
if m0_values.get('down', 0) > 0 and m0_values.get('up', 0) > 0:
    ratio_ud = m0_values['up'] / m0_values['down']
    print(f"  m_0(u) / m_0(d) = {ratio_ud:.2f}")
    print(f"  Porownania:")
    print(f"    m_3(u)/m_3(d) = {quarks['up']['m'][2]/quarks['down']['m'][2]:.2f}")
    print(f"    (m_3(u)/m_3(d))^(1/2) = {np.sqrt(quarks['up']['m'][2]/quarks['down']['m'][2]):.2f}")
    print(f"    r_21^(1/2) = {np.sqrt(206.768):.2f}")
    print(f"    phi^8 = {phi**8:.2f}")
    print(f"    m_0(u)/m_0(d) / sqrt(m_3(u)/m_3(d)) = {ratio_ud / np.sqrt(quarks['up']['m'][2]/quarks['down']['m'][2]):.4f}")

    # Czy m_0 ~ m_3^(1/2)?
    log_m0 = [np.log(m0_values['down']), np.log(m0_values['up'])]
    log_m3 = [np.log(quarks['down']['m'][2]), np.log(quarks['up']['m'][2])]
    slope = (log_m0[1] - log_m0[0]) / (log_m3[1] - log_m3[0])
    print(f"\n  d ln(m_0) / d ln(m_3) = {slope:.4f}")
    print(f"    (blisko 1/2 => m_0 ~ m_3^(1/2), blisko 1 => m_0 ~ m_3)")

# ============================================================
# PODSUMOWANIE
# ============================================================
print("\n" + "=" * 70)
print("PODSUMOWANIE")
print("=" * 70)
print("""
  Najlepsze kandydatki na formule m_0 z parametrow TGP:

  1. H11: m_0 = sqrt(alpha_s(m_3) * m_3)
     - Uzywa running alpha_s na skali m_3
     - Motywacja: energia konfinementu ~ sqrt(V_cornell * m_kinetic)
     - Predykcja: m_0(d) ~ sqrt(0.22 * 4180) ~ 30 MeV (~+40%)
                  m_0(u) ~ sqrt(0.11 * 172760) ~ 137 MeV (daleko)
     - Problem: nie odtwarza m_0(u) -- skala top jest za wysoka

  2. H7: m_0 = a_Gamma * m_3 / Phi_0
     - Czysto TGP, zero parametrow zewnetrznych
     - m_0(d) = 0.040 * 4180 / 24.66 = 6.78 MeV (-69%)
     - m_0(u) = 0.040 * 172760 / 24.66 = 280 MeV (-86%)
     - Problem: za male o czynnik ~3

  3. H3: m_0 = Lambda_QCD^2 / m_1
     - Standardowa estymacja z worka MIT
     - m_0(d) = 217^2 / 4.67 = 10087 MeV (za duze!)
     - Problem: skalowanie odwrotne do m_1, a nie m_3

  WNIOSEK: Zadna prosta kombinacja nie daje m_0 z dokladnoscia < 20%
  dla OBU sektorow jednoczesnie. To potwierdza, ze R12 pozostaje
  otwartym problemem wymagajacym pelnej symulacji N-body rury kolorowej
  w potencjale TGP rezimu III.

  NOWE OBSERWACJE:
  - d ln(m_0) / d ln(m_3) ~ 1.2 (bliskie 1, nie 1/2)
  - m_0/m_3 nie jest uniwersalne (0.0053 vs 0.0115)
  - Uniwersalna jest stala A = m_0 * m_1/m_3 ~ 0.025
  - Hipoteza do zbadania: A = a_Gamma / phi = 0.0247 (0.4% od empirycznej!)
""")

# Sprawdzenie a_Gamma / phi
A_hyp = a_Gamma / phi
print(f"  >>> NOWA HIPOTEZA: A = a_Gamma / phi = {A_hyp:.6f}")
print(f"      A_empiryczna = {A_mean_val:.6f}")
print(f"      Odchylenie: {(A_hyp/A_mean_val-1)*100:+.2f}%")
