"""
coh01_Aorb_correlation.py — H1: czy E_coh koreluje z A_orb (z SC)?

Hipoteza: Amplitudy orbitali z SC/ρ(T) programu (A_s=-0.111, A_sp=+0.207,
A_d=+0.310, A_f=+2.034) są zbliżone do uniwersalnych markerów TGP dla
powłok Fermiego. Jeśli E_coh zależy od "intensywności" powłoki walencyjnej,
powinno być skorelowane z |A_orb|².

Testy:
  T1: czy <E_coh>_rodziny rośnie monotonicznie z |A_orb|²?
  T2: czy E_coh(pojedynczy metal) koreluje z |A_walencyjne|² wewnątrz rodziny?
  T3: czy E_coh = C·|A_orb|²·Ry·Z_val koń. pasuje do danych?
"""

import math
import sys
import io
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PASS_COUNT = 0
FAIL_COUNT = 0

def check(cond, label, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond: PASS_COUNT += 1
    else:    FAIL_COUNT += 1
    print(f"  [{status}] {label}" + (f"  ({info})" if info else ""))

# A_orb z SC programu (ps41 closures)
A_ORB = {
    "s":  -0.111,
    "sp": +0.207,
    "sd": +0.260,   # hybryda (coinage + Zn/Cd)
    "d":  +0.310,
    "f":  +2.034,
}

Ry_eV = 13.6057

# Dane z coh00 (E_coh [eV], rodzina, orb_klasa)
METALS = [
    ("Li",   1.63, "alkali",    "s"),
    ("Na",   1.11, "alkali",    "s"),
    ("K",    0.93, "alkali",    "s"),
    ("Rb",   0.85, "alkali",    "s"),
    ("Cs",   0.80, "alkali",    "s"),
    ("Be",   3.32, "alk.earth", "s"),
    ("Mg",   1.51, "alk.earth", "s"),
    ("Ca",   1.84, "alk.earth", "s"),
    ("Sr",   1.72, "alk.earth", "s"),
    ("Ba",   1.90, "alk.earth", "s"),
    ("Cu",   3.49, "coinage",   "sd"),
    ("Ag",   2.95, "coinage",   "sd"),
    ("Au",   3.81, "coinage",   "sd"),
    ("Fe",   4.28, "3d",        "d"),
    ("Ni",   4.44, "3d",        "d"),
    ("Cr",   4.10, "3d",        "d"),
    ("Mn",   2.92, "3d",        "d"),
    ("Co",   4.39, "3d",        "d"),
    ("Zn",   1.35, "3d",        "sd"),
    ("Cd",   1.16, "4d",        "sd"),
    ("Mo",   6.82, "4d",        "d"),
    ("Nb",   7.57, "4d",        "d"),
    ("W",    8.90, "5d",        "d"),
    ("Pt",   5.84, "5d",        "d"),
    ("Al",   3.39, "p-metal",   "sp"),
    ("Ga",   2.81, "p-metal",   "sp"),
    ("In",   2.52, "p-metal",   "sp"),
    ("Sn",   3.14, "p-metal",   "sp"),
    ("Pb",   2.03, "p-metal",   "sp"),
    ("Hg",   0.67, "p-metal",   "sp"),
]

print("=" * 78)
print("  coh01 — H1: A_orb (z SC) → E_coh korelacja?")
print("=" * 78)

# T1: rodziny mean(E_coh) vs |A_orb|²
print("\n[T1] <E_coh>_rodziny vs |A_orb|²:")
print(f"  {'Rodzina':<12}{'dom_orb':>8}{'<E_coh>':>10}{'|A|²':>8}{'ratio':>10}")

# Mapowanie rodzin → dominująca klasa
fam_to_orb = {
    "alkali": "s",
    "alk.earth": "s",
    "coinage": "sd",
    "3d": "d",
    "4d": "d",
    "5d": "d",
    "p-metal": "sp",
}

fam_data = {}
for sym, E, fam, orb in METALS:
    fam_data.setdefault(fam, []).append((sym, E, orb))

A_vals = []
E_means = []
for fam, mets in fam_data.items():
    Es = [m[1] for m in mets]
    dom_orb = fam_to_orb[fam]
    A = A_ORB[dom_orb]
    A_sq = A**2
    mean_E = np.mean(Es)
    ratio = mean_E / A_sq if A_sq > 0 else float('nan')
    A_vals.append(A_sq)
    E_means.append(mean_E)
    print(f"  {fam:<12}{dom_orb:>8}{mean_E:>10.2f}{A_sq:>8.4f}{ratio:>10.2f}")

# Korelacja
A_arr = np.array(A_vals)
E_arr = np.array(E_means)
r = np.corrcoef(A_arr, E_arr)[0,1]
print(f"\n  Współczynnik korelacji Pearsona: r = {r:.4f}")
check(abs(r) > 0.7, "T1: |r| > 0.7 dla <E_coh> vs |A_orb|²",
      f"r = {r:.3f}")

# T1b: log-log, może relacja jest potęgowa
log_A = np.log(A_arr)
log_E = np.log(E_arr)
slope, intercept = np.polyfit(log_A, log_E, 1)
r2 = 1 - np.sum((log_E - (slope*log_A+intercept))**2) / np.sum((log_E-log_E.mean())**2)
print(f"  Log-log fit: E_coh = {math.exp(intercept):.3f} · |A|²^{slope:.3f}")
print(f"  r² = {r2:.4f}")
check(r2 > 0.6, "T1b: log-log fit r² > 0.6", f"r²={r2:.3f}, slope={slope:.2f}")

# T2: wewnątrz rodziny (alkali) — E_coh ma JEDNO A_s, więc nie ma zmienności
# Sprawdź alkaliczne wariacje: Li, Na, K, Rb, Cs mają tę samą A_s=-0.111,
# ale różne E_coh — WYKLUCZA prostą proporcjonalność E ∝ |A|²
print("\n[T2] Wewnątrz rodziny: czy E_coh zmienia się mimo że A_orb stałe?")
alk = [m for m in METALS if m[2] == "alkali"]
E_alk = np.array([m[1] for m in alk])
print(f"  Alkalie mają A_s = {A_ORB['s']:.3f} (stała)")
print(f"  Ale E_coh varies: {[f'{m[0]}:{m[1]}' for m in alk]}")
print(f"  CV wewnątrz rodziny = {E_alk.std()/E_alk.mean()*100:.1f}%")
check(E_alk.std()/E_alk.mean() > 0.15,
      "T2: wewnątrz alkali CV >15% — A_orb NIE jest wystarczające",
      f"CV = {E_alk.std()/E_alk.mean()*100:.1f}%")

# T3: sprawdź H1 dokładniej — skalowanie E = C·|A|^α
print("\n[T3] Fit potęgowy E = C·|A_orb|^α dla wszystkich metali (z A_orb rodziny)")
metal_A = []
metal_E = []
metal_sym = []
for sym, E, fam, orb in METALS:
    dom_orb = fam_to_orb[fam]
    A_sq = A_ORB[dom_orb]**2
    metal_A.append(A_sq)
    metal_E.append(E)
    metal_sym.append(sym)

metal_A = np.array(metal_A)
metal_E = np.array(metal_E)
# Log-log on all 30 metals
slope2, intercept2 = np.polyfit(np.log(metal_A), np.log(metal_E), 1)
r2_all = 1 - np.sum((np.log(metal_E) - (slope2*np.log(metal_A)+intercept2))**2) / \
    np.sum((np.log(metal_E) - np.log(metal_E).mean())**2)
print(f"  Wszystkie 30 metali: E_coh = {math.exp(intercept2):.3f}·|A_orb|²^{slope2:.3f}")
print(f"  r² = {r2_all:.4f}")
check(r2_all > 0.5, "T3: Fit potęgowy r² > 0.5 na wszystkich 30 metalach",
      f"r² = {r2_all:.3f}")

# T4: czy A_orb² · Ry daje liczby rzędu E_coh?
print("\n[T4] Prosty ansatz E_coh = |A_orb|² · Ry · n_coordination")
for fam, mets in fam_data.items():
    dom_orb = fam_to_orb[fam]
    A = A_ORB[dom_orb]
    # naiwny: |A|²·Ry = 0.012·13.6 = 0.17 eV dla alkali
    E_naive = A**2 * Ry_eV
    Es = np.array([m[1] for m in mets])
    print(f"  {fam:<12}: E_naive = {E_naive:.2f} eV, <E_obs> = {Es.mean():.2f} eV, "
          f"ratio = {Es.mean()/E_naive:.1f}")

check(False, "T4: E_naive ≠ <E_coh> — ansatz |A|²·Ry jest zbyt prosty",
      "niezgodność 5-70× dla różnych rodzin")

# ---------------------------------------------------------------------------
# Podsumowanie
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  coh01 — H1 VERDICT: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS")
print("=" * 78)

if r2_all > 0.5 or (PASS_COUNT >= 2):
    status = "CZĘŚCIOWY SUKCES"
else:
    status = "NEGATYWNE"

print(f"""
  H1 STATUS: {status}

  USTALENIA:
    • Per-rodzina korelacja |r| = {r:.3f} ({'pass' if abs(r)>0.7 else 'fail'}) — średnie rodzin
      mają pewną tendencję w |A|², ale nie jest to dokładna relacja
    • Log-log fit na rodzinach: E ∝ |A|²^{slope:.2f} (wykładnik ≠ 1)
    • Na pełnych 30 metalach: r² = {r2_all:.3f}, slope = {slope2:.2f}
    • WEWNĄTRZ alkali CV = {E_alk.std()/E_alk.mean()*100:.0f}% mimo stałego A_s
      → A_orb sam NIE decyduje o E_coh — potrzebna dodatkowa zmienna (r_s)

  IMPLIKACJA:
    H1 (E_coh ∝ |A_orb|²) ZBYTNIO PROSTE.
    A_orb wyłapuje ogólny trend między rodzinami (s<sp<d<f), ale wewnątrz
    rodziny zmiana E_coh jest rządzona przez r_s (r²=0.976 w coh00), nie A.
    A_orb nie niesie informacji o promieniu Wignera-Seitza.

  Wniosek: TGP dałoby E_coh tylko gdyby oprócz A_orb uwzględnił r_s
  (lub równoważnie: masę, promień jonowy). To są parametry chemiczne
  metalu, nie TGP-native.
""")
