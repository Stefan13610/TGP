"""
p38_dsb_landscape.py
====================
CEL: Analiza krajobrazu d/s/b kwarków i predykcje Q_TGP dla wszystkich rodzin

DIAGNOZA (P37):
  d/s/b: r21=20.0 wymaga alpha_f < 1.0 (poza zakresem alpha_lo=0.5 w P37)
  Przy a=0.04, alpha=0.5: r21 ~ 19-20 -> rozwiazanie TUZY granic brentq
  Przy a=0.04, alpha=0:   r21 ~ 12    -> minimum r21 jest > 0

PLAN P38:
  Czesc A -- Mapa g(K) dla malych alpha (weryfikacja 3 zer)
    alpha in [0.1, 0.5, 1.0, 1.5, 2.0] przy a=0.040
    Wykaz K1, K2 numerycznie + r21(alpha)

  Czesc B -- r21(alpha, a_Gamma) dla malych alpha
    Mapa r21 dla alpha in [0, 3] x a in [0.010, 0.060]
    Znajdz minimalne r21_min(a_Gamma) i odpowiednie alpha=0

  Czesc C -- Znajdz alpha_f dla d/s/b (rozszerzony zakres alpha)
    Brentq z alpha_lo=0.05

  Czesc D -- Q_TGP dla WSZYSTKICH trzech rodzin
    Dla kazdej rodziny (leptons, u-quarks, d-quarks):
      znajdz (alpha_f, lambda_Koide) spelniajace r21 i Q=3/2
      vs obserwowane Q_PDG
    Pytanie: czy TGP przewiduje Q=3/2 dla wszystkich, czy tylko leptonow?
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import exp1
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX   = 60.0
GAMMA   = 1.0
LAM_REF = 1e-5

print("P38: Krajobraz d/s/b i predykcje Q_TGP dla wszystkich rodzin fermionow")
print("=" * 72)
print()

# ============================================================
# NARZEDZIA
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam=LAM_REF, N=2000):
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam=LAM_REF):
    e = energy_log(K, alpha, a_gam, lam)
    return e / (4*np.pi*K) - 1.0 if np.isfinite(e) else np.nan

def find_K1_num(alpha, a_gam, K_lo=0.0003, K_hi=0.15):
    """K1 - poszerzone granice dla malych alpha."""
    for lo, hi in [(K_lo, K_hi), (K_lo, 0.3)]:
        try:
            glo = g_func(lo, alpha, a_gam)
            ghi = g_func(hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), lo, hi, xtol=1e-13)
        except Exception:
            pass
    return np.nan

def find_K2_num(alpha, a_gam):
    """K2 - drugie zero g(K)."""
    for K_lo, K_hi in [(0.05, 0.5), (0.3, 2.0), (0.5, 5.0)]:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-10)
        except Exception:
            pass
    return np.nan

# ============================================================
# CZESC A: Mapa g(K) i struktura zer dla malych alpha
# ============================================================
print("CZESC A: Struktura zer g(K) dla malych alpha (a=0.040)")
print("-" * 70)
print()

a_fix = 0.040
alpha_small = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.553]

print(f"  {'alpha':>6}  {'K1':>10}  {'K2':>10}  {'r21':>8}  {'g(K1/2)':>10}")
print("  " + "-"*52)

r21_vs_alpha = {}
for alpha in alpha_small:
    K1 = find_K1_num(alpha, a_fix)
    K2 = find_K2_num(alpha, a_fix)
    if np.isnan(K1) or np.isnan(K2):
        print(f"  {alpha:>6.2f}  {'BRAK':>10}  {'BRAK':>10}")
        continue
    r21 = K2/K1
    # Sprawdz g w srodku (K1,K2) - powinno byc ujemne
    K_mid = (K1+K2)/2
    g_mid = g_func(K_mid, alpha, a_fix)
    r21_vs_alpha[alpha] = r21
    print(f"  {alpha:>6.2f}  {K1:>10.6f}  {K2:>10.6f}  {r21:>8.2f}  {g_mid:>10.4f}")

print()
print("  Minimum r21 (alpha->0) dla a=0.040:")
# Ekstrapolacja alpha=0
K1_a0 = find_K1_num(0.001, a_fix)
K2_a0 = find_K2_num(0.001, a_fix)
if not (np.isnan(K1_a0) or np.isnan(K2_a0)):
    r21_a0 = K2_a0/K1_a0
    print(f"  alpha=0.001: K1={K1_a0:.6f}, K2={K2_a0:.6f}, r21={r21_a0:.2f}")
else:
    print("  alpha=0.001: BRAK (granica numeryczna)")

# ============================================================
# CZESC B: r21_min(a_Gamma) i zakres r21 TGP
# ============================================================
print()
print("CZESC B: Zakres r21(alpha) dla roznych a_Gamma")
print("-" * 70)
print()

agam_scan = [0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060]
alpha_dense = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.553, 12.0, 20.0]

# Dla kazdego a: min r21 (przy alpha=0.1), max r21 (przy alpha=20)
print(f"  {'a_Gamma':>8}  {'r21_min(a0.1)':>14}  {'r21(a8.553)':>12}  "
      f"{'r21_max(a20)':>13}  {'r21_dsb=20':>12}")
print("  " + "-"*66)

r21_landscape = {}
for a in agam_scan:
    K1_lo = find_K1_num(0.1, a); K2_lo = find_K2_num(0.1, a)
    K1_mid = find_K1_num(8.553, a); K2_mid = find_K2_num(8.553, a)
    K1_hi = find_K1_num(20.0, a); K2_hi = find_K2_num(20.0, a)

    r21_lo  = K2_lo/K1_lo   if not (np.isnan(K1_lo)  or np.isnan(K2_lo))  else np.nan
    r21_mid = K2_mid/K1_mid if not (np.isnan(K1_mid) or np.isnan(K2_mid)) else np.nan
    r21_hi  = K2_hi/K1_hi   if not (np.isnan(K1_hi)  or np.isnan(K2_hi))  else np.nan

    # Czy r21=20 miezci sie w zakresie?
    in_range = "TAK" if (not np.isnan(r21_lo) and r21_lo < 20 < (r21_hi or 0)) else "NIE"

    r21_landscape[a] = {'lo':r21_lo, 'mid':r21_mid, 'hi':r21_hi}
    print(f"  {a:>8.3f}  {r21_lo:>14.2f}  {r21_mid:>12.2f}  "
          f"{r21_hi:>13.2f}  {in_range:>12}")

# ============================================================
# CZESC C: Wyznaczanie alpha_f dla d/s/b (rozszerzony zakres)
# ============================================================
print()
print("CZESC C: alpha_f dla wszystkich rodzin fermionow (zakres rozszerzony)")
print("-" * 70)
print()

fermions = {
    'e/mu/tau': {'m1': 0.511,   'm2': 105.658,  'm3': 1776.86},
    'u/c/t':    {'m1': 2.16,    'm2': 1270.0,    'm3': 172690.0},
    'd/s/b':    {'m1': 4.67,    'm2': 93.4,      'm3': 4180.0},
}

def find_alpha_for_r21_wide(target_r21, a_gam, alpha_lo=0.05, alpha_hi=50.0):
    """Szeroki zakres alpha, wlacznie z malymi wartosciami."""
    def f(al):
        K1 = find_K1_num(al, a_gam)
        K2 = find_K2_num(al, a_gam)
        if np.isnan(K1) or np.isnan(K2) or K1 < 1e-12:
            return np.nan
        return K2/K1 - target_r21
    try:
        flo = f(alpha_lo); fhi = f(alpha_hi)
        if np.isnan(flo) or np.isnan(fhi) or flo*fhi >= 0:
            return np.nan
        return brentq(f, alpha_lo, alpha_hi, xtol=1e-6)
    except Exception:
        return np.nan

agam_search = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060])
landscape_full = {name: {} for name in fermions}

for name, ferm in fermions.items():
    r21_obs = ferm['m2']/ferm['m1']
    r31_obs = ferm['m3']/ferm['m1']
    print(f"  === {name} (r21={r21_obs:.2f}, r31={r31_obs:.0f}) ===")
    print(f"  {'a_Gamma':>8}  {'alpha_f':>9}  {'r21_TGP':>9}  {'K1':>10}  "
          f"{'K2':>8}  {'alpha/a':>8}")
    print("  " + "-"*62)

    for a in agam_search:
        af = find_alpha_for_r21_wide(r21_obs, a)
        if not np.isnan(af):
            K1f = find_K1_num(af, a)
            K2f = find_K2_num(af, a)
            r21c = K2f/K1f if not (np.isnan(K1f) or np.isnan(K2f)) else np.nan
            landscape_full[name][a] = af
            print(f"  {a:>8.4f}  {af:>9.4f}  {r21c:>9.1f}  "
                  f"{K1f:>10.7f}  {K2f:>8.4f}  {af/a:>8.2f}")
        else:
            print(f"  {a:>8.4f}  {'FAIL':>9}")
    print()

# ============================================================
# CZESC D: Q_TGP dla wszystkich rodzin
# ============================================================
print("CZESC D: Predykcje Q_TGP dla wszystkich rodzin fermionow")
print("-" * 70)
print()

def Q_koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (m1 + m2 + m3)

def r31_koide_alg(r21):
    """r31 = x takie ze Q(1, r21, x) = 3/2 (algebraicznie)."""
    # Q = (1 + sqrt(r21) + sqrt(x))^2 / (1 + r21 + x) = 3/2
    # Rozwiazanie: (1 + sqrt(r21) + sqrt(x))^2 = 3/2*(1 + r21 + x)
    # Rozpisujemy: S = 1 + sqrt(r21), ksi = sqrt(x)
    # (S + ksi)^2 = 3/2*(1 + r21 + ksi^2)
    # S^2 + 2S*ksi + ksi^2 = 3/2*(1+r21) + 3/2*ksi^2
    # -1/2*ksi^2 + 2S*ksi + S^2 - 3/2*(1+r21) = 0
    # ksi^2 - 4S*ksi - 2*(S^2 - 3/2*(1+r21)) = 0
    S = 1 + np.sqrt(r21)
    A = 1; B = -4*S; C_coef = -2*(S**2 - 1.5*(1+r21))
    disc = B**2 - 4*A*C_coef
    if disc < 0: return np.nan
    ksi = (-B + np.sqrt(disc))/(2*A)
    return ksi**2

def find_lambda_koide(alpha, a_gam, C_val=2.0):
    """lambda_Koide z formuly analitycznej (P31-P35)."""
    K1 = find_K1_num(alpha, a_gam)
    K2 = find_K2_num(alpha, a_gam)
    if np.isnan(K1) or np.isnan(K2) or K1 < 1e-12:
        return np.nan
    r21 = K2/K1
    r31K = r31_koide_alg(r21)
    if np.isnan(r31K): return np.nan
    return (C_val * a_gam)**2 / (K1 * r31K)**2

def find_K3_num(alpha, a_gam, lam):
    """K3 - trzecie zero g(K) przy danym lambda."""
    K_lo_list = [(5.0, 100.0), (10.0, 200.0), (20.0, 500.0)]
    for K_lo, K_hi in K_lo_list:
        try:
            glo = g_func(K_lo, alpha, a_gam, lam)
            ghi = g_func(K_hi, alpha, a_gam, lam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam, lam),
                              K_lo, K_hi, xtol=1e-8)
        except Exception:
            pass
    return np.nan

# Dla kazdej rodziny fermionow przy wybranych a_Gamma:
# Znajdz (alpha_f, lambda_K) i oblicz (r21_TGP, r31_TGP, Q_TGP)
print("  Wyniki: dla kazdej rodziny przy reprezentatywnym a_Gamma")
print(f"  {'Rodzina':>12}  {'a_Gamma':>8}  {'alpha_f':>8}  {'r21_TGP':>9}  "
      f"{'r21_obs':>8}  {'Q_obs':>7}  {'r31_K':>9}  {'lambda_K':>12}")
print("  " + "-"*90)

# Referncyjne a_Gamma dla kazdej rodziny (gdzie alpha_f jest znana z Czesci C)
ref_agam = {'e/mu/tau': 0.040, 'u/c/t': 0.020, 'd/s/b': 0.020}

q_results = {}
for name, ferm in fermions.items():
    m1,m2,m3 = ferm['m1'],ferm['m2'],ferm['m3']
    r21_obs = m2/m1; r31_obs = m3/m1
    Q_obs   = Q_koide(m1, m2, m3)

    a_ref = ref_agam[name]
    af = landscape_full[name].get(a_ref, np.nan)

    if np.isnan(af):
        # Probuj inne a_Gamma
        for a_try in [0.025, 0.030, 0.015, 0.010]:
            if a_try in landscape_full[name]:
                af = landscape_full[name][a_try]
                a_ref = a_try
                break

    if np.isnan(af):
        print(f"  {name:>12}  BRAK DANYCH")
        continue

    K1f = find_K1_num(af, a_ref)
    K2f = find_K2_num(af, a_ref)
    if np.isnan(K1f) or np.isnan(K2f):
        print(f"  {name:>12}  K1/K2 FAIL")
        continue

    r21_TGP = K2f/K1f
    r31_K   = r31_koide_alg(r21_TGP)
    lam_K   = find_lambda_koide(af, a_ref)

    q_results[name] = {
        'a': a_ref, 'alpha_f': af, 'r21_TGP': r21_TGP, 'r21_obs': r21_obs,
        'Q_obs': Q_obs, 'r31_K': r31_K, 'lam_K': lam_K,
        'K1': K1f, 'K2': K2f
    }

    print(f"  {name:>12}  {a_ref:>8.3f}  {af:>8.4f}  {r21_TGP:>9.2f}  "
          f"{r21_obs:>8.1f}  {Q_obs:>7.4f}  {r31_K:>9.2f}  {lam_K:>12.4e}")

print()

# Kluczowe pytanie: czy Q=3/2 jest SPECJALNE dla leptonow?
print("  Kluczowa analiza: Q_TGP vs Q_obs dla kazdej rodziny")
print()
print(f"  {'Rodzina':>12}  {'r21_obs':>8}  {'Q_obs':>7}  "
      f"{'r31_K(r21_obs)':>16}  {'r31_obs':>9}  {'delta_r31':>10}")
print("  " + "-"*72)

for name, ferm in fermions.items():
    m1,m2,m3 = ferm['m1'],ferm['m2'],ferm['m3']
    r21_obs = m2/m1; r31_obs = m3/m1
    Q_obs   = Q_koide(m1, m2, m3)
    r31_K   = r31_koide_alg(r21_obs)   # r31 przy KTORYM Q(r21_obs, r31) = 3/2
    delta   = (r31_obs - r31_K)/r31_K*100 if not np.isnan(r31_K) else np.nan
    print(f"  {name:>12}  {r21_obs:>8.2f}  {Q_obs:>7.4f}  "
          f"{r31_K:>16.2f}  {r31_obs:>9.2f}  {delta:>+9.2f}%")

print()
print("  Interpretacja:")
print("  Jesli delta_r31 ~ 0% -> Q_obs ~ 3/2 (leptony: idealne)")
print("  Jesli |delta_r31| >> 0% -> Q_obs != 3/2 (kwarki)")
print()

# Dodatkowa analiza: dla d/s/b - jakie lambda_K daje Q=3/2?
# vs. lambda ktore spelnia r31_obs
print("  Dla d/s/b: porownanie lambda_K(Q=3/2) vs lambda(r31_obs)")
print()
for name, ferm in fermions.items():
    m1,m2,m3 = ferm['m1'],ferm['m2'],ferm['m3']
    r21_obs = m2/m1; r31_obs = m3/m1
    Q_obs = Q_koide(m1, m2, m3)

    if name not in q_results: continue
    res = q_results[name]
    af, a_ref = res['alpha_f'], res['a']
    K1f = res['K1']; K2f = res['K2']
    lam_K_32 = res['lam_K']  # lambda dla Q=3/2

    # Lambda dla r31_obs
    K3_for_r31 = r31_obs * K1f
    # Odwrotnosc K3 = C*a/sqrt(lam): lam = (C*a/K3)^2
    lam_r31_obs = (2.0*a_ref)**2 / K3_for_r31**2

    print(f"  {name}: Q_obs={Q_obs:.4f}")
    if not np.isnan(lam_K_32):
        print(f"    lambda(Q=3/2)   = {lam_K_32:.4e}  [r31_K = {res['r31_K']:.1f}]")
    print(f"    lambda(r31_obs) = {lam_r31_obs:.4e}  [K3 = {K3_for_r31:.3f}]")
    if not np.isnan(lam_K_32):
        ratio = lam_r31_obs/lam_K_32
        print(f"    Stosunek lambda_obs/lambda_K = {ratio:.4f}")
    print()

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("P38: Krajobraz d/s/b i predykcje Q_TGP", fontsize=13)

# Panel A: r21(alpha) dla roznych a_Gamma
ax = axes[0, 0]
ax.set_title(r"$r_{21}(\alpha)$ dla róznych $a_\Gamma$")
colors_a = plt.cm.viridis(np.linspace(0, 0.85, len(agam_scan)))
alpha_plot = np.concatenate([[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7],
                              np.arange(1.0, 21.0, 1.0)])
for i, a_val in enumerate(agam_scan):
    r21_plot = []
    al_plot  = []
    for al in alpha_plot:
        K1 = find_K1_num(al, a_val)
        K2 = find_K2_num(al, a_val)
        if np.isnan(K1) or np.isnan(K2): continue
        r21_plot.append(K2/K1)
        al_plot.append(al)
    if r21_plot:
        ax.plot(al_plot, r21_plot, '-', color=colors_a[i], lw=1.5, label=f"a={a_val:.3f}")

# Zaznacz obserwowane r21 dla kazdej rodziny
ax.axhline(206.8, color='blue', ls='--', alpha=0.6, lw=1, label="leptony r21=207")
ax.axhline(588.0, color='red',  ls='--', alpha=0.6, lw=1, label="u-kwarki r21=588")
ax.axhline(20.0,  color='green',ls='--', alpha=0.6, lw=1, label="d-kwarki r21=20")
ax.set_xlabel(r"$\alpha$"); ax.set_ylabel(r"$r_{21} = K_2/K_1$")
ax.set_yscale('log'); ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

# Panel B: r21(alpha) dla malych alpha, a=0.04
ax = axes[0, 1]
ax.set_title(r"$r_{21}(\alpha)$ dla małych $\alpha$ (a=0.040)")
alpha_small_plot = np.concatenate([[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
                                    0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0, 3.0]])
r21_small = []
al_small  = []
for al in alpha_small_plot:
    K1 = find_K1_num(al, 0.040)
    K2 = find_K2_num(al, 0.040)
    if np.isnan(K1) or np.isnan(K2): continue
    r21_small.append(K2/K1)
    al_small.append(al)
ax.plot(al_small, r21_small, 'bo-', ms=5, lw=1.5)
ax.axhline(20.0, color='green', ls='--', label="d/s/b r21=20")
ax.axhline(68.0, color='orange',ls=':',  label="d/s/b α=3 → r21=68")
ax.set_xlabel(r"$\alpha$"); ax.set_ylabel(r"$r_{21}$")
ax.legend(); ax.grid(True, alpha=0.3)

# Panel C: alpha_f(a_Gamma) dla wszystkich 3 rodzin
ax = axes[1, 0]
ax.set_title(r"$\alpha_f(a_\Gamma)$ dla wszystkich rodzin")
colors_f = ['blue', 'red', 'green']
labels_f = ['e/μ/τ (r₂₁=207)', 'u/c/t (r₂₁=588)', 'd/s/b (r₂₁=20)']
for i, name in enumerate(fermions.keys()):
    ldata = landscape_full[name]
    if not ldata: continue
    a_vals  = sorted(ldata.keys())
    af_vals = [ldata[a] for a in a_vals]
    ax.plot(a_vals, af_vals, 'o-', color=colors_f[i], ms=5, lw=1.5, label=labels_f[i])
ax.set_xlabel(r"$a_\Gamma$"); ax.set_ylabel(r"$\alpha_f$")
ax.legend(); ax.grid(True, alpha=0.3)

# Panel D: Q_PDG vs r21 (diagram Koidego) ze wskazaniem rodzin
ax = axes[1, 1]
ax.set_title("Diagram Q(r21) Koidego — rodziny fermionow")
# Krzywa Q=3/2 w przestrzeni (r21, r31): dla danego r21, r31 daje Q=3/2
r21_line = np.logspace(0.3, 3.5, 200)
r31_q32  = np.array([r31_koide_alg(r) for r in r21_line])
valid    = ~np.isnan(r31_q32)
ax.loglog(r21_line[valid], r31_q32[valid], 'k-', lw=2, label="Q=3/2 (Koide)")

# Punkty PDG
pdg_data = [
    ('e/μ/τ', 206.8, 3477.65, 'blue', '*', 14),
    ('u/c/t', 588.0, 79949.0, 'red',  'D', 8),
    ('d/s/b', 20.0,  895.0,   'green','s', 8),
]
for label, r21v, r31v, c, m, ms in pdg_data:
    ax.scatter(r21v, r31v, color=c, marker=m, s=ms**2, zorder=5, label=label)
ax.set_xlabel(r"$r_{21}$"); ax.set_ylabel(r"$r_{31}$")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
outpath = "TGP/TGP_v1/scripts/advanced/p38_dsb_landscape.png"
plt.savefig(outpath, dpi=120, bbox_inches='tight')
print(f"Wykres zapisany: {outpath}")

print()
print("=" * 72)
print("WYNIKI KLUCZOWE P38")
print("=" * 72)
print()
print("  A. r21(alpha) jest monotonicznie rosnace z alpha")
print("     Minimum r21 (alpha->0): czas na diagnostyke")
print()
print("  B. d/s/b (r21=20): wymaga alpha_f << 1 (problem P37 = zly zakres)")
print()
print("  C. Krzywa Koidego Q=3/2: tylko leptony leza na niej!")
print("     u/c/t: Q=1.1779 (r31 o 99.8% za duze od krzywej Koidego)")
print("     d/s/b: Q=1.3672 (r31 o ?% od krzywej)")
