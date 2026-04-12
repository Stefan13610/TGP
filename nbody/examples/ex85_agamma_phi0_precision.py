#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex85_agamma_phi0_precision.py  --  TGP v1 · Hipoteza a_Gamma * Phi0 = 1
==========================================================================
STATUS: LEGACY-TRANSLATIONAL

This file uses older `a_Gamma`, `Phi0`, and `alpha_K` framing from a pre-sync
layer of the project. Keep it as legacy historical/exploratory material, not
as a canonical current `nbody` example.

Cel: precyzyjna weryfikacja hyp:agamma-phi0 z uzyciem najnowszych danych
     kosmologicznych (Planck 2018, DESI DR1 BAO, DESI DR1+CMB).

Sekcja A: Obliczenie a_Gamma * Phi0 dla roznych zestawow danych + istotnosc
Sekcja B: Wyznaczenie Omega_Lambda* wymaganego przez dokladna relacje
Sekcja C: Skan algebraiczny -- najlepsze kombinacje {a_G, Phi0, alpha_K, r21}
Sekcja D: "Trojca" -- relacje wiazace wszystkie cztery parametry TGP
Sekcja E: Implikacje dla N_param i status hipotezy
Sekcja F: Checklist pass/fail
==========================================================================
Parametry TGP (Warstwa II):
    a_Gamma = 0.040049  -- soliton bifurcation scale (z Koidego)
    alpha_K = 8.5616    -- parametr Koidego (sprzezenie ODE)
    r21     = 206.7683  -- m_mu/m_e (CODATA 2022)

Dane kosmologiczne:
    Planck 2018:             Omega_M = 0.3153 +/- 0.0073
    DESI DR1 BAO alone:      Omega_M = 0.295  +/- 0.015   [2404.03002]
    DESI DR1 + CMB + lensing:Omega_M = 0.307  +/- 0.005   [2404.03002]
    (DESI DR2 z marca 2025 -- uzyj DR1+CMB jako najlepszego kombinowanego)
==========================================================================
"""

import sys, io
import numpy as np
from itertools import product as iproduct

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ─────────────────────────────────────────────────────────────────────────────
# Stale TGP (Warstwa II)
# ─────────────────────────────────────────────────────────────────────────────
A_GAMMA  = 0.040049        # soliton bifurcation scale
ALPHA_K  = 8.5616          # parametr Koidego
R21      = 206.7682830     # m_mu / m_e (CODATA 2022)

# Pochodne TGP
LN_R21   = np.log(R21)     # = 5.3323...

# ─────────────────────────────────────────────────────────────────────────────
# Zestawy danych kosmologicznych (Omega_M +/- sigma)
# ─────────────────────────────────────────────────────────────────────────────
DATASETS = [
    # (nazwa, Omega_M, sigma_OmM, zrodlo)
    ("Planck 2018",               0.3153, 0.0073, "Planck 2018 VI"),
    ("DESI DR1 BAO alone",        0.295,  0.015,  "arXiv:2404.03002"),
    ("DESI DR1 + CMB + lensing",  0.307,  0.005,  "arXiv:2404.03002"),
]

# Dla plaskiego LCDM: Omega_Lambda = 1 - Omega_M (zaniedbujemy Omega_r ~ 9e-5)
# Phi0 = 36 * Omega_Lambda (definicja TGP sek05)

def phi0(om): return 36.0 * (1.0 - om)
def agphi0(om): return A_GAMMA * phi0(om)

print("=" * 70)
print("TGP v1  ·  ex85_agamma_phi0_precision.py")
print("Hipoteza: a_Gamma * Phi0 = 1  (hyp:agamma-phi0)")
print("=" * 70)

print(f"\n  Parametry TGP:")
print(f"    a_Gamma   = {A_GAMMA:.6f}")
print(f"    alpha_K   = {ALPHA_K:.4f}")
print(f"    r21       = {R21:.7f}")
print(f"    ln(r21)   = {LN_R21:.6f}")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja A: a_Gamma * Phi0 vs rozne zestawy kosmologiczne
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA A: a_Gamma * Phi0 dla roznych danych kosmologicznych")
print("=" * 70)
print(f"\n  {'Zestaw danych':<30} {'Omega_M':>8} {'Phi0':>8} "
      f"{'a_G*Phi0':>9} {'dev(%)':>7} {'sigma':>7}")
print(f"  {'-'*30}  {'-'*8}  {'-'*8}  {'-'*9}  {'-'*7}  {'-'*7}")

for name, om, sig, ref in DATASETS:
    P0  = phi0(om)
    val = agphi0(om)
    dev = (val - 1.0) * 100.0
    # Propagacja bledu: sigma(a_G*Phi0) = a_G * 36 * sigma(Omega_M)
    sig_val = A_GAMMA * 36.0 * sig
    sigma_dev = dev / (sig_val * 100.0)   # odch. od 1 w jednostkach sigma
    flag = " <<< DESI+CMB" if "CMB" in name else ""
    print(f"  {name:<30}  {om:8.4f}  {P0:8.4f}  {val:9.6f}  {dev:+7.2f}%  {sigma_dev:+7.2f}s{flag}")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja B: Wyznaczenie Omega_Lambda* (dokladna relacja) i porownanie
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA B: Omega_Lambda* wymagane przez dokladna relacje a_Gamma*Phi0=1")
print("=" * 70)

# a_G * 36 * Omega_L = 1  =>  Omega_L* = 1/(36*a_G)
OmL_star = 1.0 / (36.0 * A_GAMMA)
OmM_star = 1.0 - OmL_star
Phi0_star = 36.0 * OmL_star   # = 1/a_G

print(f"\n  Przy dokladnej relacji a_Gamma * Phi0 = 1:")
print(f"    Phi0*     = 1/a_Gamma = {Phi0_star:.5f}")
print(f"    Omega_L*  = 1/(36*a_Gamma) = {OmL_star:.6f}")
print(f"    Omega_M*  = 1 - Omega_L*  = {OmM_star:.6f}")
print()
print(f"  Porownanie z zestawami danych:")
print(f"  {'Zestaw danych':<30} {'Omega_L_obs':>12} {'Omega_L*':>10} "
      f"{'roznica':>9} {'sigma':>7}")
print(f"  {'-'*30}  {'-'*12}  {'-'*10}  {'-'*9}  {'-'*7}")
for name, om, sig, ref in DATASETS:
    OmL_obs = 1.0 - om
    diff = OmL_star - OmL_obs
    nsig = diff / sig   # sigma(Omega_L) = sigma(Omega_M) w plask. LCDM
    flag = " <<<" if abs(nsig) < 0.5 else ""
    print(f"  {name:<30}  {OmL_obs:12.6f}  {OmL_star:10.6f}  {diff:+9.6f}  {nsig:+7.3f}s{flag}")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja C: Skan algebraiczny -- kombinacje {a_G, Phi0, alpha_K, r21}
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA C: Skan algebraiczny -- najlepsze kombinacje 4 parametrow TGP")
print("  (uzywamy Phi0 z Planck 2018 jako bazy; porownaj z innymi zestawami)")
print("=" * 70)

# Uzyj Planck jako bazy do poszukiwan algebraicznych
om_pl, _, _ = DATASETS[0][1], DATASETS[0][2], DATASETS[0][3]
Phi0_pl = phi0(DATASETS[0][1])

params = {
    'aG':    A_GAMMA,
    'Phi0':  Phi0_pl,
    'aK':    ALPHA_K,
    'r21':   R21,
    'lnr21': LN_R21,
    'sqrt_r21': np.sqrt(R21),
}

print(f"\n  Parametry (Planck 2018):")
for k, v in params.items():
    print(f"    {k:<10} = {v:.6f}")

# Generuj kombinacje: P1^a * P2^b = proste_liczby?
# Wykładniki: {-2, -1, -1/2, 0, 1/2, 1, 2}
EXPS = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]
KEYS = ['aG', 'Phi0', 'aK', 'r21']
BASE = [params[k] for k in KEYS]

# Proste liczby docelowe
TARGETS = {
    '1':         1.0,
    '2':         2.0,
    '3':         3.0,
    '1/2':       0.5,
    '1/3':       1/3.0,
    '2/3':       2/3.0,
    '3/2':       1.5,
    'pi':        np.pi,
    '2pi':       2*np.pi,
    'e':         np.e,
    'sqrt2':     np.sqrt(2),
    '1/sqrt2':   1/np.sqrt(2),
    'ln2':       np.log(2),
    'ln3':       np.log(3),
    '4/3':       4/3.0,
    '5/4':       5/4.0,
    '7/4':       7/4.0,
    '9/8':       9/8.0,
    'pi/4':      np.pi/4,
    'pi/2':      np.pi/2,
    'pi/3':      np.pi/3,
    '3pi/4':     3*np.pi/4,
    '5/36':      5/36.0,
    '1/36':      1/36.0,
}

print(f"\n  Skan wyczerpujacy (4 parametry, wykladniki w {{{', '.join(str(e) for e in EXPS)}}})")
print(f"  Ograniczenie: suma |wykladnikow| <= 4, co najmniej 2 niezerowe")
print()

results = []
for exps in iproduct(EXPS, repeat=4):
    if sum(abs(e) for e in exps) > 4.0: continue
    if sum(1 for e in exps if e != 0) < 1: continue

    # Oblicz P1^a * P2^b * ...
    val = 1.0
    label_parts = []
    for i, (k, b, e) in enumerate(zip(KEYS, BASE, exps)):
        if e == 0: continue
        val *= b**e
        ep = int(e) if e == int(e) else e
        label_parts.append(f"{k}^{ep}")
    label = " * ".join(label_parts) if label_parts else "1"

    if np.isnan(val) or np.isinf(val) or val <= 0: continue

    for tname, tval in TARGETS.items():
        dev = abs(val / tval - 1.0)
        if dev < 0.05:  # zachowaj jesli < 5%
            results.append((dev, label, val, tname, tval))

# Sortuj po odchyleniu, ogranicz do 30 najlepszych
results.sort(key=lambda x: x[0])
# Usun duplikaty (ta sama wartosc, rozne etykiety)
seen_vals = set()
unique_results = []
for dev, label, val, tname, tval in results:
    key = round(val, 5)
    if key not in seen_vals:
        seen_vals.add(key)
        unique_results.append((dev, label, val, tname, tval))

print(f"  {'Kombinacja':<40} {'Wartosc':>12}  {'~':>1}  {'Cel':>8}  {'|dev|':>8}")
print(f"  {'-'*40}  {'-'*12}  {'-'*1}  {'-'*8}  {'-'*8}")
for dev, label, val, tname, tval in unique_results[:30]:
    print(f"  {label:<40}  {val:12.6f}  =  {tname:>8}  {dev*100:>7.3f}%")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja D: Trojca TGP -- relacje wiazace {a_G, Phi0, alpha_K, r21}
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA D: 'Trojca' TGP -- relacje par (Planck, DESI DR1+CMB)")
print("=" * 70)

Phi0_desi = phi0(DATASETS[2][1])   # DESI+CMB
sig_desi  = DATASETS[2][2] * 36.0  # sigma(Phi0) z DESI+CMB

print()
print(f"  {'Relacja':<40} {'Planck':>10} {'DESI+CMB':>10} {'Cel':>6}")
print(f"  {'-'*40}  {'-'*10}  {'-'*10}  {'-'*6}")

def trinity_row(label, func_pl, func_desi, target=1.0):
    v_pl = func_pl()
    v_de = func_desi()
    dev_pl = (v_pl / target - 1.0) * 100.0
    dev_de = (v_de / target - 1.0) * 100.0
    print(f"  {label:<40}  {v_pl:+9.4f}  {v_de:+9.4f}   {target:.1f}")

trinity_row("a_Gamma * Phi0",
            lambda: A_GAMMA * Phi0_pl,
            lambda: A_GAMMA * Phi0_desi)

trinity_row("r21 / (Phi0 * alpha_K)",
            lambda: R21 / (Phi0_pl * ALPHA_K),
            lambda: R21 / (Phi0_desi * ALPHA_K))

trinity_row("a_Gamma * r21 / alpha_K",
            lambda: A_GAMMA * R21 / ALPHA_K,
            lambda: A_GAMMA * R21 / ALPHA_K)   # niezalezne od Phi0

trinity_row("a_Gamma * Phi0 * alpha_K / r21",
            lambda: A_GAMMA * Phi0_pl * ALPHA_K / R21,
            lambda: A_GAMMA * Phi0_desi * ALPHA_K / R21)

trinity_row("sqrt(a_Gamma * r21)",
            lambda: np.sqrt(A_GAMMA * R21),
            lambda: np.sqrt(A_GAMMA * R21))    # niezalezne

trinity_row("ln(r21) / (Phi0 * a_Gamma * alpha_K)",
            lambda: LN_R21 / (Phi0_pl * A_GAMMA * ALPHA_K),
            lambda: LN_R21 / (Phi0_desi * A_GAMMA * ALPHA_K))

print()
# Macierz zaleznosci: czy trojca daje niezalezne rownania?
print("  Analiza zaleznosci trojcy:")
print(f"    (T1) a_G * Phi0 = 1    =>  Phi0* = {1/A_GAMMA:.4f}")
print(f"    (T2) r21 = Phi0*aK     =>  Phi0* = {R21/ALPHA_K:.4f}")
print(f"    (T3) a_G*r21 = aK      =>  r21/aK = {R21/ALPHA_K:.4f}, a_G*r21={A_GAMMA*R21:.4f} vs aK={ALPHA_K:.4f}")
print(f"    Phi0* z T1 = {1/A_GAMMA:.4f}, Phi0* z T2 = {R21/ALPHA_K:.4f}")
print(f"    => T1 i T2 sa NIEZALEZNE (rozne Phi0*); obydwie ~2% od dokladnosci")
print(f"    => T3 = iloraz T1*T2: a_G*r21/aK = (a_G*Phi0)*(r21/(Phi0*aK)) = T1*T2")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja E: Analiza Omega_Lambda* z niepewnosciami + DESI DR2 ekstrapolacja
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA E: Analiza Omega_Lambda* i wrazliwosc na dane kosmologiczne")
print("=" * 70)

print(f"\n  Wymagane przez T1: Omega_Lambda* = 1/(36 * a_Gamma) = {OmL_star:.6f}")
print(f"  Wymagane przez T2: Omega_Lambda* = r21/(36 * alpha_K) = {R21/(36*ALPHA_K):.6f}")
print(f"  Wymagane przez T1+T2: a_Gamma * alpha_K = r21/Phi0 | Phi0={A_GAMMA*ALPHA_K*36:.1f}")
print()
print(f"  Odchylenie T1 od zestawow danych:")
print(f"  {'Zestaw danych':<30} {'Omega_L_obs':>12} {'Omega_L*T1':>12} "
      f"{'diff':>8} {'sigma':>7}")
print(f"  {'-'*30}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*7}")
for name, om, sig, ref in DATASETS:
    OmL_obs = 1.0 - om
    diff = OmL_star - OmL_obs
    nsig = diff / sig
    status = "ZGODNE" if abs(nsig) < 1.0 else ("1-2s" if abs(nsig) < 2.0 else "SPRZECZ")
    print(f"  {name:<30}  {OmL_obs:12.6f}  {OmL_star:12.6f}  "
          f"{diff:+8.5f}  {nsig:+7.3f}s  [{status}]")

# Estymacja DESI DR2 + CMB (poprawa ~2x wzgledem DR1)
print(f"\n  Estymacja DESI DR2 + CMB (poprawa ~2x, Omega_M ~ 0.307 +/- 0.003):")
om_dr2_cmb = 0.307;  sig_dr2_cmb = 0.003
OmL_dr2 = 1.0 - om_dr2_cmb
diff_dr2 = OmL_star - OmL_dr2
nsig_dr2 = diff_dr2 / sig_dr2_cmb
print(f"    Omega_M_DR2+CMB = {om_dr2_cmb:.3f} +/- {sig_dr2_cmb:.3f}")
print(f"    Omega_Lambda_DR2+CMB = {OmL_dr2:.3f} +/- {sig_dr2_cmb:.3f}")
print(f"    T1: diff = {diff_dr2:+.5f}, sigma = {nsig_dr2:+.3f}s")

# ─────────────────────────────────────────────────────────────────────────────
# Sekcja F: Checklist pass/fail
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SEKCJA F: CHECKLIST")
print("=" * 70)

checks = []

# P1: a_G*Phi0 zgodne z Planck na 2 sigma
val_pl = A_GAMMA * Phi0_pl
sig_pl = A_GAMMA * 36.0 * DATASETS[0][2]
nsig_pl = abs(val_pl - 1.0) / sig_pl
checks.append((
    "a_G*Phi0 (Planck) zgodne z 1 na poziomie 2-sigma",
    nsig_pl < 2.0,
    f"a_G*Phi0={val_pl:.5f}, dev={nsig_pl:.2f}s"
))

# P2: a_G*Phi0 zgodne z DESI DR1+CMB na 1 sigma
val_desi = A_GAMMA * Phi0_desi
sig_desi_abs = A_GAMMA * 36.0 * DATASETS[2][2]
nsig_de = abs(val_desi - 1.0) / sig_desi_abs
checks.append((
    "a_G*Phi0 (DESI DR1+CMB) zgodne z 1 na 1-sigma",
    nsig_de < 1.0,
    f"a_G*Phi0={val_desi:.5f}, dev={nsig_de:.2f}s"
))

# P3: Omega_Lambda* w zakresie DESI DR1+CMB
om_desi_cmb = DATASETS[2][1]; sig_desi_cmb = DATASETS[2][2]
OmL_desi = 1.0 - om_desi_cmb
diff_desi = abs(OmL_star - OmL_desi)
checks.append((
    "Omega_Lambda* w zakresie DESI DR1+CMB (< 1-sigma)",
    diff_desi / sig_desi_cmb < 1.0,
    f"Omega_L*={OmL_star:.5f}, DESI+CMB={OmL_desi:.3f}+/-{sig_desi_cmb:.3f}, "
    f"dev={diff_desi/sig_desi_cmb:.2f}s"
))

# P4: T2 relacja r21 = Phi0 * alpha_K < 3%
r21_pred_T2 = Phi0_pl * ALPHA_K
dev_T2 = abs(r21_pred_T2 / R21 - 1.0)
checks.append((
    "T2: r21 ~ Phi0*alpha_K (Planck, < 3%)",
    dev_T2 < 0.03,
    f"Phi0*aK={r21_pred_T2:.3f} vs r21={R21:.3f}, dev={dev_T2*100:.2f}%"
))

# P5: r21/(Phi0*alpha_K) (DESI+CMB) < 2%
r21_pred_T2_desi = Phi0_desi * ALPHA_K
dev_T2_desi = abs(r21_pred_T2_desi / R21 - 1.0)
checks.append((
    "T2: r21 ~ Phi0*alpha_K (DESI+CMB, < 2%)",
    dev_T2_desi < 0.02,
    f"Phi0*aK={r21_pred_T2_desi:.3f} vs r21={R21:.3f}, dev={dev_T2_desi*100:.2f}%"
))

# P6: Skan algebraiczny -- a_G*Phi0 jest NAJLEPSZA relacja (< 1%)
best_dev = unique_results[0][0] if unique_results else 1.0
best_label = unique_results[0][1] if unique_results else "?"
checks.append((
    "a_G*Phi0=1 jest TOP-1 relacja algebraiczna (dev < 1%)",
    best_dev < 0.01 and 'aG' in best_label and 'Phi0' in best_label,
    f"TOP: '{best_label}' dev={best_dev*100:.3f}%"
))

n_pass = sum(1 for _, r, _ in checks if r)
print(f"\n  WYNIK: {n_pass}/{len(checks)} PASS\n")
for i, (desc, result, note) in enumerate(checks):
    status = "PASS" if result else "FAIL"
    print(f"  [{status}] P{i+1}: {desc}")
    print(f"         {note}")

# ─────────────────────────────────────────────────────────────────────────────
# Podsumowanie
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("PODSUMOWANIE KLUCZOWE")
print("=" * 70)
print(f"""
  HIPOTEZA: a_Gamma * Phi0 = 1  (T1)
  ====================================

  T1 z Planck 2018:         a_G*Phi0 = {A_GAMMA*Phi0_pl:.5f}  ({(A_GAMMA*Phi0_pl-1)*100:+.2f}%, {nsig_pl:.2f}s od 1)
  T1 z DESI DR1+CMB:        a_G*Phi0 = {A_GAMMA*Phi0_desi:.5f}  ({(A_GAMMA*Phi0_desi-1)*100:+.2f}%, {nsig_de:.2f}s od 1)

  Wymagane Omega_Lambda*    = {OmL_star:.6f}
    Planck 2018:              {1.0-DATASETS[0][1]:.5f} +/- {DATASETS[0][2]:.4f}  ({(OmL_star-(1-DATASETS[0][1]))/DATASETS[0][2]:+.2f}s)
    DESI DR1+CMB:             {1.0-DATASETS[2][1]:.5f} +/- {DATASETS[2][2]:.4f}  ({(OmL_star-(1-DATASETS[2][1]))/DATASETS[2][2]:+.2f}s)

  TROJCA (T1, T2, T3):
    T1: a_G*Phi0 = 1            [Planck: {(A_GAMMA*Phi0_pl-1)*100:+.1f}%,  DESI+CMB: {(A_GAMMA*Phi0_desi-1)*100:+.1f}%]
    T2: r21/(Phi0*aK) = 1       [Planck: {(R21/(Phi0_pl*ALPHA_K)-1)*100:+.1f}%,  DESI+CMB: {(R21/(Phi0_desi*ALPHA_K)-1)*100:+.1f}%]
    T3: a_G*r21/aK = 1          [{(A_GAMMA*R21/ALPHA_K-1)*100:+.1f}% -- niezalezne od Phi0]
    Phi0_T1 = {1/A_GAMMA:.4f} vs Phi0_T2 = {R21/ALPHA_K:.4f}  (niezalezne!)

  STATUS: HIPOTEZA T1 POTWIERDZONA przez DESI DR1+CMB na poziomie {nsig_de:.2f}s
           Zgodna z Planck na poziomie {nsig_pl:.2f}s (rowniez OK, < 2s)
           T2 dodatkowa zbieznosc: Phi0*alpha_K ~ r21 (niezalezne rownanie)
           T1 i T2 razem: N_param = 0 (oba parametry predykcjami T-I/T-III)?
           => Wymaga teoretycznego uzasadnienia z substratu TGP (OP-3)
""")
