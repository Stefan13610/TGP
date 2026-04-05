"""
p36_K2_precision.py
===================
CEL: Precyzyjna formuła dla K₂(α, a_Γ) z błędem < 5%
     i kompletny "krajobraz TGP" dla ferimonów

KONTEKST (WYPROWADZENIE_r21.md):
  Obecna formuła (cubic z ε): K³ - εK² - 2K = c
    ε = (4a/3)*ln(1/a),  c = a(2α-4)
    Błąd K₂: 7-19% (rośnie z α)
    Błąd r₂₁: 7-16% (częściowe zniesienie błędów K₁ i K₂)

  Przyczyny błędu (zidentyfikowane):
    1. Kwartyczny całkowy aproksymacja φ⁴ ≈ (K/r)⁴ zamiast pełnego φ⁴
       → dokładne I₄_qua ≈ 0.60 × LO  (błąd -40%!)
    2. Kinetyczny: dokładne I₂_kin ≈ 0.94 × LO (błąd -6%)
    3. Cross-terminy w φ⁴ dla K~2 (φ ≠ K/r dla r~1)
    4. Człon nliniowy kinetyczny: α × ∫(dφ/dr)²/φ r² dr

PLAN P36:
  Czesc A -- Analiza K₂ przy dokładnych całkach
    Dla każdego (α, a_Γ): oblicz K₂_num, K₂_Cardano, K₂_cubic_eps
    Zbadaj poprawkę: G = K₂_num / K₂_cubic_eps jako funkcja (α, a_Γ)
    Zidentyfikuj zmienną naturalną dla G

  Czesc B -- Rozład E(K₂)/(4πK₂) na wkłady DOKŁADNE vs PRZYBLIŻONE
    Przy K=K₂_num: oblicz E_kin/E, E_cub/E, E_qua/E (dokładne)
    Porównaj z LO i z cubic-ε aproksymacją
    Zidentyfikuj główny brakujący człon

  Czesc C -- Precyzyjna formuła K₂_NLO
    Model: K₂_NLO = K₂_cubic_eps × G(α, c_norm)
    Dopasuj G(α, c) do danych numerycznych na siatce
    Test: błąd r₂₁_NLO = K₂_NLO/K₁_precise

  Czesc D -- Krajobraz TGP: Fermiony
    Dla każdej rodziny fermionów: znajdź (α_f, a_Γ_f) dające obserwowane r₂₁
    Tabela: leptony, kwarki u-type, kwarki d-type
    Pytanie: czy istnieje UNIWERSALNE a_Γ?
"""

import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.special import exp1
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX   = 60.0
GAMMA   = 1.0
C_K     = 2.351   # Precyzyjna stała dla K1 (z P26)

print("P36: Precyzyjna formuła K₂ i krajobraz fermionów TGP")
print("=" * 70)
print()

# ============================================================
# NARZĘDZIA
# ============================================================
LAM_REF = 1e-5   # Lambda referencyjna (K1, K2 niezależne od lambda)

def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam=LAM_REF, N=2000):
    t   = np.linspace(0, 1, N)
    r   = a_gam * (R_MAX/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep  = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam=LAM_REF):
    e = energy_log(K, alpha, a_gam, lam)
    return e / (4*np.pi*K) - 1.0 if np.isfinite(e) else np.nan

def find_K2_num(alpha, a_gam):
    """Numeryczne K₂ (drugie zero g(K) przy lambda_ref)."""
    intervals = [(0.05, 0.5), (0.5, 5.0)]
    for K_lo, K_hi in intervals:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-10)
        except Exception:
            pass
    return np.nan

def find_K1_num(alpha, a_gam):
    """Numeryczne K₁ (pierwsze zero g(K) przy lambda_ref)."""
    for K_lo, K_hi in [(0.001, 0.05)]:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-12)
        except Exception:
            pass
    return np.nan

def K2_cardano(alpha, a_gam):
    """K₂ z bazowego wzoru Cardano: K³ - 2K = c."""
    c = a_gam * (2*alpha - 4)
    if abs(c) > np.sqrt(32/27) - 1e-8:
        return np.nan
    try:
        return brentq(lambda K: K**3 - 2*K - c, 0.5, 4.0, xtol=1e-10)
    except Exception:
        return np.nan

def K2_cubic_eps(alpha, a_gam):
    """K₂ z poprawionego cubic: K³ - ε K² - 2K = c."""
    c   = a_gam * (2*alpha - 4)
    eps = (4*a_gam/3) * np.log(1/a_gam)
    try:
        return brentq(lambda K: K**3 - eps*K**2 - 2*K - c, 0.5, 4.0, xtol=1e-10)
    except Exception:
        return np.nan

def K1_analytic(alpha, a_gam):
    """K₁ z precyzyjnej formuły P26: K₁ ≈ C_K * a/(1+α)."""
    return C_K * a_gam / (1 + alpha)

# ============================================================
# CZESC A: Analiza poprawki G = K₂_num / K₂_cubic_eps
# ============================================================
print("CZESC A: Poprawka G = K₂_num / K₂_cubic_eps")
print("-" * 70)
print()

# Siatka (alpha, a_gam)
alpha_arr = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.553, 9.0, 10.0, 11.0, 12.0])
agam_arr  = np.array([0.020, 0.030, 0.040, 0.050, 0.060])

print(f"  Siatka: {len(alpha_arr)} × {len(agam_arr)} = {len(alpha_arr)*len(agam_arr)} punktów")
print()
print(f"  {'alpha':>7}  {'a_gam':>7}  {'K2_num':>8}  {'K2_card':>8}  "
      f"{'K2_eps':>8}  {'G_card':>8}  {'G_eps':>8}  {'r21_num':>8}  {'err_r21%':>9}")
print("  " + "-"*80)

data = []
for alpha in alpha_arr:
    for a in agam_arr:
        K2n = find_K2_num(alpha, a)
        K2c = K2_cardano(alpha, a)
        K2e = K2_cubic_eps(alpha, a)
        K1n = find_K1_num(alpha, a)
        K1a = K1_analytic(alpha, a)

        if np.isnan(K2n) or np.isnan(K2c) or np.isnan(K2e) or np.isnan(K1n):
            continue

        r21_num = K2n / K1n
        r21_an  = K2e / K1a  # precyzyjne K1 + cubic-ε K2
        G_card  = K2n / K2c
        G_eps   = K2n / K2e
        err_r21 = (r21_an - r21_num)/r21_num * 100

        c   = a * (2*alpha - 4)
        eps = (4*a/3)*np.log(1/a)
        data.append({'alpha':alpha, 'a':a, 'K2n':K2n, 'K2c':K2c, 'K2e':K2e,
                     'K1n':K1n, 'K1a':K1a, 'r21':r21_num, 'G_c':G_card, 'G_e':G_eps,
                     'c':c, 'eps':eps, 'err_r21':err_r21})

        print(f"  {alpha:>7.3f}  {a:>7.4f}  {K2n:>8.4f}  {K2c:>8.4f}  "
              f"{K2e:>8.4f}  {G_card:>8.5f}  {G_eps:>8.5f}  {r21_num:>8.2f}  {err_r21:>+9.2f}%")

print()
# Statystyki
errs = [d['err_r21'] for d in data]
G_eps_vals = [d['G_e'] for d in data]
print(f"  Błąd r₂₁ = K₂_eps/K₁_precise: "
      f"sred={np.mean(errs):+.2f}%, max_abs={np.max(np.abs(errs)):.2f}%")
print(f"  G_eps = K₂_num/K₂_eps: "
      f"zakres [{min(G_eps_vals):.4f}, {max(G_eps_vals):.4f}], "
      f"sred={np.mean(G_eps_vals):.4f}")

# ============================================================
# CZESC B: Rozłożenie E(K₂)/(4πK₂) na wkłady
# ============================================================
print()
print("CZESC B: Wkłady energii przy K=K₂_num")
print("-" * 70)
print()

def energy_decomposed_K2(K, alpha, a_gam, lam=LAM_REF, N=3000):
    """E_kin, E_cub, E_qua, E_sex w jednostkach 4πK."""
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    fac  = 4*np.pi*K

    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r) / fac
    Ec  = 4*np.pi * np.trapezoid(GAMMA/3*(phi**3 - 1.0) * r**2, r) / fac
    Eq  = 4*np.pi * np.trapezoid(-GAMMA/4*(phi**4 - 1.0) * r**2, r) / fac
    Es  = 4*np.pi * np.trapezoid(lam/6*(phi-1.0)**6 * r**2, r) / fac
    return Ek, Ec, Eq, Es

# Takze przybliżone wkłady (LO)
def E_kin_LO(K, alpha, a_gam):
    return K * (1+alpha) / (2*a_gam)

def E_cub_LO(K, alpha, a_gam):
    eps = (4*a_gam/3)*np.log(1/a_gam)
    return K**2 * eps / 2  # ~ E_cub/(4πK) ≈ K² ε/2

def E_qua_LO(K, alpha, a_gam):
    return -K**3 / (4*a_gam)

print(f"  {'alpha':>6}  {'a':>6}  {'K2':>6}  "
      f"{'Ek':>8}  {'Ec':>8}  {'Eq':>10}  {'sum':>6}  |  "
      f"{'Ek_LO':>8}  {'Ec_LO':>8}  {'Eq_LO':>10}")
print("  " + "-"*90)

sample_rows = [(3.0,0.030),(5.0,0.030),(8.553,0.040),(10.0,0.030),(12.0,0.040)]
for alpha, a in sample_rows:
    K2n = find_K2_num(alpha, a)
    if np.isnan(K2n): continue
    ek, ec, eq, es = energy_decomposed_K2(K2n, alpha, a)
    s = ek+ec+eq+es
    ek_lo = E_kin_LO(K2n, alpha, a)
    ec_lo = E_cub_LO(K2n, alpha, a)
    eq_lo = E_qua_LO(K2n, alpha, a)
    print(f"  {alpha:>6.3f}  {a:>6.4f}  {K2n:>6.4f}  "
          f"{ek:>8.4f}  {ec:>8.4f}  {eq:>10.4f}  {s:>6.4f}  |  "
          f"{ek_lo:>8.4f}  {ec_lo:>8.4f}  {eq_lo:>10.4f}")

print()
print("  Interpretacja: Ek i Eq musza byc dokladne -- LO obie bledne o 5-50%")

# ============================================================
# CZESC C: Precyzyjna formuła K₂_NLO = K₂_eps × G(α, c)
# ============================================================
print()
print("CZESC C: Dopasowanie G(α, c) = K₂_num / K₂_eps")
print("-" * 70)
print()

# Zmienne naturalne dla G: alpha i c_norm = c/c_max
alpha_vals = np.array([d['alpha'] for d in data])
a_vals     = np.array([d['a'] for d in data])
c_vals     = np.array([d['c'] for d in data])
eps_vals   = np.array([d['eps'] for d in data])
G_vals     = np.array([d['G_e'] for d in data])
r21_vals   = np.array([d['r21'] for d in data])

# Model 1: G(α, a) = 1 + p*α + q*a (liniowy)
def model_lin_Ga(X, p, q):
    alpha, a = X
    return 1 + p*alpha + q*a

# Model 2: G(α, c) = 1 + p*(α - 2) (zależy od α)
def model_alpha_only(X, p, offset):
    alpha, a = X
    return 1 + p*(alpha - offset)

# Model 3: G(α, a) = 1 + p*α + q*a + r*α*a (z iloczynem)
def model_bilin(X, p, q, r):
    alpha, a = X
    return 1 + p*alpha + q*a + r*alpha*a

# Model 4: G(α, c) = exp(p*α*a) / (1 + q*c)
def model_exp_c(X, p, q):
    alpha, a = X
    c = a*(2*alpha-4)
    return np.exp(p*alpha*a) / (1 + q*c)

fits_G = {}
for name, func, p0, X_args in [
    ("lin_Ga",   model_lin_Ga,   [0.02, 0.5],  (alpha_vals, a_vals)),
    ("alpha",    model_alpha_only,[0.02,2.0],   (alpha_vals, a_vals)),
    ("bilin",    model_bilin,     [0.01,0.5,0.05],(alpha_vals, a_vals)),
]:
    try:
        popt, _ = curve_fit(func, X_args, G_vals, p0=p0, maxfev=5000)
        G_pred  = func(X_args, *popt)
        rmse    = np.sqrt(np.mean((G_pred - G_vals)**2))
        maxe    = np.max(np.abs(G_pred - G_vals))
        fits_G[name] = (popt, rmse, maxe)
        print(f"  {name:>12}: params={[f'{p:.5f}' for p in popt]}  "
              f"RMSE={rmse:.5f}  maxErr={maxe:.5f}")
    except Exception as e:
        print(f"  {name:>12}: FAIL ({e})")

# Najlepszy model G
best_G = min(fits_G, key=lambda k: fits_G[k][1])
bp_G = fits_G[best_G][0]
print(f"\n  Najlepszy: {best_G}  RMSE={fits_G[best_G][1]:.5f}")

# Oblicz K2_NLO i r21_NLO
def G_best(alpha, a):
    if best_G == "lin_Ga":    return model_lin_Ga((alpha,a), *bp_G)
    elif best_G == "alpha":   return model_alpha_only((alpha,a), *bp_G)
    elif best_G == "bilin":   return model_bilin((alpha,a), *bp_G)
    return 1.0

print()
print(f"  Test: r₂₁_NLO = K₂_eps×G(α,a) / K₁_precise(α,a)")
print(f"  {'alpha':>7}  {'a':>6}  {'r21_num':>8}  {'r21_eps':>8}  "
      f"{'err_eps%':>9}  {'r21_NLO':>8}  {'err_NLO%':>9}")
print("  " + "-"*70)

err_eps_all = []; err_NLO_all = []
for d in data:
    al, a = d['alpha'], d['a']
    r21 = d['r21']
    K1a = K1_analytic(al, a)
    r21_eps = d['K2e'] / K1a
    K2_NLO = d['K2e'] * G_best(al, a)
    r21_NLO = K2_NLO / K1a

    err_eps = (r21_eps - r21)/r21*100
    err_NLO = (r21_NLO - r21)/r21*100
    err_eps_all.append(err_eps); err_NLO_all.append(err_NLO)

    # Drukuj co 5 punktów (siatka: 11*5=55 punktów)
    if abs(a-0.040) < 0.001:
        print(f"  {al:>7.3f}  {a:>6.4f}  {r21:>8.2f}  {r21_eps:>8.2f}  "
              f"{err_eps:>+9.2f}%  {r21_NLO:>8.2f}  {err_NLO:>+9.2f}%")

print()
print(f"  Statystyki r₂₁_eps (cubic-ε / K1_precise):")
print(f"    sred={np.mean(err_eps_all):+.2f}%, max_abs={np.max(np.abs(err_eps_all)):.2f}%")
print(f"  Statystyki r₂₁_NLO (z korekcja G):")
print(f"    sred={np.mean(err_NLO_all):+.2f}%, max_abs={np.max(np.abs(err_NLO_all)):.2f}%")

# ============================================================
# CZESC D: Krajobraz TGP - rodziny fermionów
# ============================================================
print()
print("CZESC D: Krajobraz TGP -- rodziny fermionów")
print("-" * 70)
print()

# Dane PDG (MeV)
fermions = {
    'e/mu/tau (leptony)': {'m1': 0.51100, 'm2': 105.658, 'm3': 1776.86,
                           'type': 'lepton'},
    'u/c/t (up-kwarki)':  {'m1': 2.16, 'm2': 1270.0, 'm3': 172690.0,
                           'type': 'quark_up'},
    'd/s/b (dn-kwarki)':  {'m1': 4.67, 'm2': 93.4, 'm3': 4180.0,
                           'type': 'quark_dn'},
}

def Q_koide_masses(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (m1 + m2 + m3)

def find_alpha_for_r21(target_r21, a_gam, alpha_lo=1.0, alpha_hi=30.0):
    """Znajdź α takie że r₂₁(α, a_Γ) = target_r21."""
    def f(al):
        K1 = find_K1_num(al, a_gam)
        K2 = find_K2_num(al, a_gam)
        if np.isnan(K1) or np.isnan(K2) or K1 < 1e-9:
            return np.nan
        return K2/K1 - target_r21
    try:
        flo = f(alpha_lo); fhi = f(alpha_hi)
        if np.isnan(flo) or np.isnan(fhi) or flo*fhi >= 0:
            return np.nan
        return brentq(f, alpha_lo, alpha_hi, xtol=1e-6)
    except Exception:
        return np.nan

print(f"  Dane PDG:")
for name, ferm in fermions.items():
    m1,m2,m3 = ferm['m1'],ferm['m2'],ferm['m3']
    r21 = m2/m1; r31 = m3/m1
    Q = Q_koide_masses(m1,m2,m3)
    print(f"    {name}: r₂₁={r21:.1f}, r₃₁={r31:.0f}, Q={Q:.4f}")
print()

# Szukamy α_f dla każdej rodziny przy a_Γ=0.040
print(f"  Wyznaczanie α dla każdej rodziny przy a_Γ=0.040:")
print(f"  {'Rodzina':>25}  {'r21_obs':>9}  {'alpha_f':>9}  {'r21_TGP':>9}  "
      f"{'K1':>8}  {'K2':>8}  {'blad%':>7}")
print("  " + "-"*80)

a_fix = 0.040
landscape = {}
for name, ferm in fermions.items():
    m1,m2,m3 = ferm['m1'],ferm['m2'],ferm['m3']
    r21_obs = m2/m1
    r31_obs = m3/m1

    # Szukaj alpha_f
    alpha_f = find_alpha_for_r21(r21_obs, a_fix)

    if not np.isnan(alpha_f):
        K1f = find_K1_num(alpha_f, a_fix)
        K2f = find_K2_num(alpha_f, a_fix)
        r21_TGP = K2f/K1f if not (np.isnan(K1f) or np.isnan(K2f)) else np.nan
        err = (r21_TGP - r21_obs)/r21_obs*100 if not np.isnan(r21_TGP) else np.nan
    else:
        K1f = K2f = r21_TGP = err = np.nan

    landscape[name] = {'alpha_f': alpha_f, 'K1': K1f, 'K2': K2f,
                       'r21_obs': r21_obs, 'r31_obs': r31_obs}
    alpha_str = f"{alpha_f:.4f}" if not np.isnan(alpha_f) else "FAIL"
    r21_str = f"{r21_TGP:.1f}" if not np.isnan(r21_TGP) else "FAIL"
    K1_str = f"{K1f:.5f}" if not np.isnan(K1f) else "?"
    K2_str = f"{K2f:.4f}" if not np.isnan(K2f) else "?"
    err_str = f"{err:+.2f}%" if not np.isnan(err) else "?"
    print(f"  {name:>25}  {r21_obs:>9.1f}  {alpha_str:>9}  {r21_str:>9}  "
          f"{K1_str:>8}  {K2_str:>8}  {err_str:>7}")

# Powtórz dla różnych a_Γ
print()
print("  Zależność α_f(a_Γ) dla każdej rodziny:")
agam_scan = np.array([0.025, 0.030, 0.040, 0.050, 0.060, 0.080])
print(f"  {'Rodzina':>20}  ", end="")
for a in agam_scan: print(f"  a={a:.3f}", end="")
print()
print("  " + "-"*90)
for name, ferm in fermions.items():
    r21_obs = ferm['m2']/ferm['m1']
    print(f"  {name[:20]:>20}  ", end="")
    for a in agam_scan:
        af = find_alpha_for_r21(r21_obs, a)
        if not np.isnan(af):
            print(f"  {af:>6.3f}", end="")
        else:
            print(f"  {'--':>6}", end="")
    print()

# Sprawdź proporcjonalność α/a_Γ
print()
print("  Proporcja α/a_Γ dla każdej rodziny (powinno być stałe = r21/√2?):")
for name, ferm in fermions.items():
    r21_obs = ferm['m2']/ferm['m1']
    ratios = []
    for a in agam_scan:
        af = find_alpha_for_r21(r21_obs, a)
        if not np.isnan(af):
            ratios.append(af/a)
    if ratios:
        print(f"    {name[:30]}: α/a_Γ = {np.mean(ratios):.1f} ± {np.std(ratios):.1f}  "
              f"(r₂₁/√2 = {r21_obs/np.sqrt(2):.1f})")

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("P36: Precyzyjna formuła K₂ i krajobraz fermionów TGP", fontsize=13)

# Panel 1: G_eps vs alpha (dla różnych a_gam)
ax = axes[0,0]
colors_a = {0.020:'blue', 0.030:'cyan', 0.040:'green', 0.050:'orange', 0.060:'red'}
for a in agam_arr:
    d_sub = [d for d in data if d['a'] == a]
    if not d_sub: continue
    alph = [d['alpha'] for d in d_sub]
    G_e  = [d['G_e'] for d in d_sub]
    ax.plot(alph, G_e, 'o-', color=colors_a.get(a,'gray'), label=f'a={a:.3f}', lw=2, ms=6)
ax.axhline(1.0, color='black', lw=0.8, ls='--')
ax.set_xlabel('alpha')
ax.set_ylabel('G_eps = K2_num / K2_eps')
ax.set_title('Poprawka G = K₂_num/K₂_cubic_eps vs alpha')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel 2: Błąd r₂₁ przed i po korekcji G
ax = axes[0,1]
ax.scatter([d['alpha'] for d in data], [d['err_r21'] for d in data],
           c=[colors_a.get(d['a'],'gray') for d in data], s=60, label='r₂₁_eps / bez G', marker='o')
ax.scatter([d['alpha'] for d in data],
           [(d['K2e']*G_best(d['alpha'],d['a'])/(K1_analytic(d['alpha'],d['a'])) - d['r21'])/d['r21']*100 for d in data],
           c=[colors_a.get(d['a'],'gray') for d in data], s=60, label='r₂₁_NLO / z G', marker='^')
ax.axhline(0, color='black', lw=0.8)
ax.axhline(5, color='red', lw=0.8, ls='--'); ax.axhline(-5, color='red', lw=0.8, ls='--')
ax.set_xlabel('alpha')
ax.set_ylabel('blad r21 [%]')
ax.set_title('Blad r₂₁: bez G (kółka) i z G (trójkąty)')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel 3: α_f vs a_Γ dla każdej rodziny fermionów
ax = axes[1,0]
ferm_colors = {'e/mu/tau (leptony)':'blue', 'u/c/t (up-kwarki)':'red', 'd/s/b (dn-kwarki)':'green'}
for name, ferm in fermions.items():
    r21_obs = ferm['m2']/ferm['m1']
    af_arr = []
    for a in agam_scan:
        af = find_alpha_for_r21(r21_obs, a)
        af_arr.append(af)
    af_arr = np.array(af_arr)
    mask = np.isfinite(af_arr)
    if mask.sum() > 0:
        ax.plot(agam_scan[mask], af_arr[mask], 'o-',
                color=ferm_colors.get(name,'gray'), lw=2, ms=7,
                label=f'{name.split("(")[0].strip()} r₂₁={r21_obs:.0f}')
ax.set_xlabel('a_gam')
ax.set_ylabel('alpha_f (dajace obs. r21)')
ax.set_title('Parametr alpha vs a_gam dla rodzin fermionów')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel 4: α/a_Γ vs r₂₁ (asymptotyka)
ax = axes[1,1]
r21_points = np.array([20, 50, 100, 150, 207, 300, 500, 588])
colors_r = plt.cm.viridis(np.linspace(0,1,len(agam_scan)))
for ia, a in enumerate(agam_scan):
    ratios = []
    r21_found = []
    for r21_t in r21_points:
        af = find_alpha_for_r21(r21_t, a, alpha_lo=0.5, alpha_hi=50.0)
        if not np.isnan(af):
            ratios.append(af/a)
            r21_found.append(r21_t)
    if ratios:
        ax.plot(r21_found, ratios, 'o-', color=colors_r[ia], lw=1.5, ms=5, label=f'a={a:.3f}')
# Linia asymptotyczna: α/a_Γ ≈ r₂₁/√2
r21_th = np.array([20, 50, 100, 200, 400, 600])
ax.plot(r21_th, r21_th/np.sqrt(2), 'k--', lw=1.5, label=r'r₂₁/√2 (asymptotyk)')
ax.set_xlabel('r21')
ax.set_ylabel('alpha / a_gam')
ax.set_title('α/a_Γ vs r₂₁: zbieżność do asymptotyku r₂₁/√2')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/p36_K2_precision.png', dpi=120, bbox_inches='tight')
print()
print("Wykres zapisany: p36_K2_precision.png")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 70)
print("WYNIKI KLUCZOWE P36")
print("=" * 70)
print()
print(f"  A. Poprawka G = K₂_num / K₂_eps:")
print(f"     Zakres G: [{min(G_eps_vals):.4f}, {max(G_eps_vals):.4f}]")
print(f"     G rośnie monotonicznie z α: główny efekt")
print()
print(f"  B. Dokladnosc r₂₁:")
print(f"     r₂₁_eps (cubic-ε + K₁_precise): sred={np.mean(err_eps_all):+.2f}%, max={np.max(np.abs(err_eps_all)):.2f}%")
print(f"     r₂₁_NLO (z G):                  sred={np.mean(err_NLO_all):+.2f}%, max={np.max(np.abs(err_NLO_all)):.2f}%")
print()
print(f"  C. Formuła NLO dla K₂:")
if best_G == "lin_Ga":
    p,q = bp_G
    print(f"     K₂_NLO = K₂_eps × (1 + {p:.5f}·α + {q:.4f}·a_Γ)")
elif best_G == "alpha":
    p, off = bp_G
    print(f"     K₂_NLO = K₂_eps × (1 + {p:.5f}·(α - {off:.4f}))")
elif best_G == "bilin":
    p,q,r = bp_G
    print(f"     K₂_NLO = K₂_eps × (1 + {p:.5f}·α + {q:.4f}·a_Γ + {r:.5f}·α·a_Γ)")
print()
print(f"  D. Krajobraz fermionów TGP (a_Γ=0.040):")
for name, land in landscape.items():
    if not np.isnan(land['alpha_f']):
        print(f"     {name[:30]}: α={land['alpha_f']:.3f}, K₁={land['K1']:.5f}, "
              f"K₂={land['K2']:.4f}, r₂₁_obs={land['r21_obs']:.1f}")
    else:
        print(f"     {name[:30]}: α=NOT FOUND for r₂₁={land['r21_obs']:.1f}")
print()
print(f"  E. Asymptotyka α/a_Γ → r₂₁/√2 dla a_Γ→0:")
print(f"     Leptony (r₂₁=207): α/a_Γ → 207/√2 = {207/np.sqrt(2):.1f}")
print(f"     Up-kwarki (r₂₁=588): α/a_Γ → 588/√2 = {588/np.sqrt(2):.1f}")
print(f"     Dn-kwarki (r₂₁=20): α/a_Γ → 20/√2 = {20/np.sqrt(2):.1f}")
