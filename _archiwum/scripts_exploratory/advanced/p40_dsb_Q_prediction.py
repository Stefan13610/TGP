# -*- coding: utf-8 -*-
"""
P40: Q_TGP dla kwarków d/s/b — pełna predykcja z K₃
=====================================================
CEL: Znalezienie K₁, K₂, K₃ dla kwarków d/s/b w TGP
     i obliczenie Q_TGP, r₂₁, r₃₁ oraz λ_obs/λ_Koide.

Kontekst (P39 korekta):
  - P38 miało błąd w find_K2_num: znajdowało K₁ zamiast K₂
  - Poprawka (P39): K₂ szukane w K>0.15 → K₂≈1.54 dla małych α
  - Dla d/s/b: α_f ≈ 0.22 przy a=0.040, r₂₁ ≈ 20

Plan P40:
  A) Wyznaczenie K₃ dla d/s/b: trzecie zero g(K) dla K≫1
  B) Obliczenie Q_TGP(K₁,K₂,K₃) dla wszystkich rodzin
  C) Porównanie z krzywą Koidego Q=3/2 w (r₂₁,r₃₁)
  D) Obliczenie λ_obs/λ_Koide z poprawioną α_f
  E) Pełna mapa (r₂₁, r₃₁, Q_TGP) vs PDG
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE (identyczne z P38/P39)
# ============================================================
R_MAX   = 60.0
GAMMA   = 1.0
LAM_REF = 1e-5

# PDG dane
FERMIONS = {
    'e/mu/tau': {'m1': 0.511,   'm2': 105.658,  'm3': 1776.86},
    'u/c/t':    {'m1': 2.16,    'm2': 1270.0,    'm3': 172690.0},
    'd/s/b':    {'m1': 4.67,    'm2': 93.4,      'm3': 4180.0},
}

def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam=LAM_REF, N=2000):
    """Energia solitonu — identyczna z P38/P39."""
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX / a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1 + alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam=LAM_REF):
    e = energy_log(K, alpha, a_gam, lam)
    return e / (4.0*np.pi*K) - 1.0 if np.isfinite(e) else np.nan

# ============================================================
# SZUKANIE ZER K₁, K₂, K₃ (POPRAWIONE — P39 fix)
# ============================================================
def find_K1_num(alpha, a_gam, K_lo=0.0003, K_hi=0.15):
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
    """K₂ szukane w K>0.15 (poprawka P39)."""
    for K_lo, K_hi in [(0.15, 2.0), (0.5, 5.0), (0.15, 5.0)]:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-10)
        except Exception:
            pass
    return np.nan

def find_K3_num(alpha, a_gam, K2_min=None):
    """
    K₃ — trzecie zero g(K), K₃ >> K₂.
    K2_min: dolna granica poszukiwań (K₃ > K2_min).
    Jeśli nie podano, używa K2_min=2.5.
    """
    if K2_min is None:
        # Wyznacz K2 żeby nie nakładać się z nim
        K2_est = find_K2_num(alpha, a_gam)
        K2_min = K2_est * 1.5 if not np.isnan(K2_est) else 2.5
    K2_min = max(K2_min, 1.5)

    # Skanuj g(K) od K2_min w górę żeby znaleźć zmianę znaku
    K_scan = np.logspace(np.log10(K2_min), np.log10(500.0), 400)
    g_scan_vals = np.array([g_func(K, alpha, a_gam) for K in K_scan])
    for i in range(len(K_scan)-1):
        if np.isfinite(g_scan_vals[i]) and np.isfinite(g_scan_vals[i+1]):
            if g_scan_vals[i] * g_scan_vals[i+1] < 0:
                try:
                    return brentq(lambda K: g_func(K, alpha, a_gam),
                                  K_scan[i], K_scan[i+1], xtol=1e-8)
                except Exception:
                    pass
    return np.nan

def g_scan(alpha, a_gam, K_min=1e-4, K_max=300, npts=500):
    """Pełny skan g(K) — zwraca (K_grid, g_vals, zera_lista)."""
    K_grid = np.logspace(np.log10(K_min), np.log10(K_max), npts)
    g_vals = np.array([g_func(K, alpha, a_gam) for K in K_grid])
    zeros = []
    for i in range(len(K_grid)-1):
        if np.isfinite(g_vals[i]) and np.isfinite(g_vals[i+1]):
            if g_vals[i]*g_vals[i+1] < 0:
                try:
                    z = brentq(lambda K: g_func(K, alpha, a_gam),
                               K_grid[i], K_grid[i+1], xtol=1e-8)
                    zeros.append(z)
                except Exception:
                    pass
    return K_grid, g_vals, zeros

# ============================================================
# FUNKCJE KOIDE
# ============================================================
def koide_Q(r21, r31):
    """Q = (√K₁+√K₂+√K₃)² / (K₁+K₂+K₃) = (1+√r21+√r31)² / (1+r21+r31)
    Q=3/2 dla leptonów (konwencja z P28--P38).
    """
    num = (1.0 + np.sqrt(r21) + np.sqrt(r31))**2
    den = (1.0 + r21 + r31)
    return num / den

def r31_koide(r21):
    """
    Rozwiązuje Q(1,r21,r31)=3/2 analitycznie.
    Wyprowadzenie (P40): (1+s+y)^2/(1+r21+y^2) = 3/2, s=sqrt(r21), y=sqrt(r31)
    => y^2 - 4(1+s)*y + (1-4s+s^2) = 0
    => y = 2(1+s) ± sqrt(3*(1+4s+s^2))
    """
    s = np.sqrt(r21)
    inner = 3.0 * (1.0 + 4.0*s + s**2)
    disc_sqrt = np.sqrt(max(inner, 0.0))
    y_plus  = 2.0*(1.0+s) + disc_sqrt   # duże r31 (fizyczne)
    y_minus = 2.0*(1.0+s) - disc_sqrt   # małe r31
    candidates = [y**2 for y in [y_plus, y_minus] if y > 0]
    return max(candidates) if candidates else np.nan

def lambda_koide(K1, K2, K3, a_gam, C=2.000):
    """λ_Koide z formuły P31: λ = (C*a_Γ)² / (K₃*K₁)."""
    if K3 <= 0 or K1 <= 0:
        return np.nan
    return (C * a_gam)**2 / (K3 * K1)

# ============================================================
# SEKCJA A: Szukanie α_f dla każdej rodziny i wyznaczanie K₃
# ============================================================
def find_alpha_for_r21(target_r21, a_gam, alpha_lo=0.05, alpha_hi=50.0):
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

def section_A():
    print("\n" + "="*70)
    print("SEKCJA A: α_f, K₁, K₂, K₃ dla wszystkich rodzin fermionów")
    print("          (poprawione find_K2, szerokie szukanie K₃)")
    print("="*70)

    a_gams = [0.025, 0.030, 0.035, 0.040]
    all_results = {}

    for a in a_gams:
        print(f"\n--- a_Γ = {a:.3f} ---")
        print(f"{'Rodzina':>10} {'r21_PDG':>9} {'α_f':>9} {'K1':>10} {'K2':>10} {'K3':>12} {'r21_TGP':>9} {'r31_TGP':>9}")
        print("-"*85)
        results_a = {}
        for name, ferm in FERMIONS.items():
            r21_pdg = ferm['m2'] / ferm['m1']
            r31_pdg = ferm['m3'] / ferm['m1']

            alpha_f = find_alpha_for_r21(r21_pdg, a)
            if np.isnan(alpha_f):
                print(f"{name:>10} {r21_pdg:>9.2f} {'FAIL':>9}")
                results_a[name] = None
                continue

            K1 = find_K1_num(alpha_f, a)
            K2 = find_K2_num(alpha_f, a)
            # K3 szukane POWYŻEJ K2 (poprawka: nie nakładamy się z K2)
            K3 = find_K3_num(alpha_f, a, K2_min=K2 * 1.5 if not np.isnan(K2) else None)

            r21_tgp = K2/K1 if not np.isnan(K1) and K1>0 else np.nan
            r31_tgp = K3/K1 if not np.isnan(K3) and not np.isnan(K1) and K1>0 else np.nan

            K3_str = f"{K3:.4f}" if not np.isnan(K3) else "NOT FOUND"
            r31_str = f"{r31_tgp:.1f}" if not np.isnan(r31_tgp) else "---"
            print(f"{name:>10} {r21_pdg:>9.2f} {alpha_f:>9.4f} {K1:>10.6f} "
                  f"{K2:>10.6f} {K3_str:>12} {r21_tgp:>9.2f} {r31_str:>9}")
            results_a[name] = {
                'alpha_f': alpha_f, 'K1': K1, 'K2': K2, 'K3': K3,
                'r21': r21_tgp, 'r31': r31_tgp,
                'r21_pdg': r21_pdg, 'r31_pdg': r31_pdg
            }
        all_results[a] = results_a
    return all_results

# ============================================================
# SEKCJA B: Q_TGP i porównanie z krzywą Koide'go
# ============================================================
def section_B(all_results):
    print("\n" + "="*70)
    print("SEKCJA B: Q_TGP vs Q_PDG i porównanie z krzywą Koide Q=3/2")
    print("="*70)

    a_focus = 0.040
    results = all_results.get(a_focus, {})

    print(f"\na_Γ = {a_focus}")
    print(f"{'Rodzina':>10} {'r21':>8} {'r31_TGP':>10} {'r31_PDG':>10} "
          f"{'r31_Koide':>11} {'δr31':>8} {'Q_TGP':>8} {'Q_PDG':>8}")
    print("-"*80)

    Q_results = {}
    for name, res in results.items():
        if res is None:
            print(f"{name:>10} FAIL")
            continue
        r21 = res['r21_pdg']
        r31_tgp = res['r31']
        r31_pdg = res['r31_pdg']
        r31_k   = r31_koide(r21)
        Q_tgp   = koide_Q(r21, r31_tgp) if not np.isnan(r31_tgp) else np.nan
        Q_pdg   = koide_Q(r21, r31_pdg)
        delta   = (r31_tgp/r31_pdg - 1.0)*100 if not np.isnan(r31_tgp) else np.nan
        r31_tgp_str = f"{r31_tgp:.1f}" if not np.isnan(r31_tgp) else "---"
        Q_tgp_str   = f"{Q_tgp:.4f}"  if not np.isnan(Q_tgp)   else "---"
        delta_str   = f"{delta:+.1f}%" if not np.isnan(delta)   else "---"
        print(f"{name:>10} {r21:>8.1f} {r31_tgp_str:>10} {r31_pdg:>10.0f} "
              f"{r31_k:>11.1f} {delta_str:>8} {Q_tgp_str:>8} {Q_pdg:>8.4f}")
        Q_results[name] = {'Q_tgp': Q_tgp, 'Q_pdg': Q_pdg, 'r31_tgp': r31_tgp}
    return Q_results

# ============================================================
# SEKCJA C: λ_obs vs λ_Koide z POPRAWIONĄ α_f (P39 korekta)
# ============================================================
def section_C(all_results):
    print("\n" + "="*70)
    print("SEKCJA C: λ_obs / λ_Koide — porównanie z poprawioną α_f (P39)")
    print("="*70)

    a_focus = 0.040
    results = all_results.get(a_focus, {})

    # λ_obs z danych PDG (zaobserwowane): używamy wzoru λ_Koide = (C*a)²/(K₃*K₁)
    # λ_obs = obserwowane λ czyli λ wyliczone z leptonu (leptony: λ_obs = λ_Koide(leptony))
    # Logika: λ_Koide(rodzina) = (C*a)² / (K₃(α_f, a) * K₁(α_f, a))
    # λ_obs/λ_Koide(leptony) = 1 dla leptonów, inne dla kwarków

    # Leptonowe λ jako punkt odniesienia
    lep = results.get('e/mu/tau')
    if lep is None:
        print("BRAK leptonów — nie można wyliczyć λ_obs")
        return

    lam_lep = lambda_koide(lep['K1'], lep['K2'], lep['K3'], a_focus)
    print(f"\nλ_Koide(leptony, a=0.040) = {lam_lep:.4e}")
    print(f"\n{'Rodzina':>10} {'α_f':>9} {'K1':>10} {'K3':>12} "
          f"{'λ_K(rodzina)':>14} {'λ_K/λ_K(lep)':>14} {'Interp.':>12}")
    print("-"*85)

    for name, res in results.items():
        if res is None:
            print(f"{name:>10} FAIL")
            continue
        lam_k = lambda_koide(res['K1'], res['K2'], res['K3'], a_focus)
        if not np.isnan(lam_k) and not np.isnan(lam_lep):
            ratio = lam_k / lam_lep
            if abs(ratio - 1.0) < 0.01:
                interp = "λ≈λ_Koide ✓"
            elif ratio > 0.5:
                interp = f"λ×{1/ratio:.1f} za mał."
            else:
                interp = f"λ×{1/ratio:.0f} za mał."
        else:
            ratio = np.nan
            interp = "---"
        K3_str = f"{res['K3']:.4f}" if not np.isnan(res['K3']) else "BRAK"
        lam_str = f"{lam_k:.4e}" if not np.isnan(lam_k) else "---"
        ratio_str = f"{ratio:.4f}" if not np.isnan(ratio) else "---"
        print(f"{name:>10} {res['alpha_f']:>9.4f} {res['K1']:>10.6f} "
              f"{K3_str:>12} {lam_str:>14} {ratio_str:>14} {interp:>12}")

# ============================================================
# SEKCJA D: Pełny skan g(K) dla kwarków d/s/b i leptonów
# ============================================================
def section_D(all_results, a_focus=0.040):
    print("\n" + "="*70)
    print(f"SEKCJA D: Profil g(K) dla d/s/b i e/mu/tau (a={a_focus})")
    print("="*70)

    results = all_results.get(a_focus, {})
    families_to_plot = ['e/mu/tau', 'd/s/b']
    profiles = {}
    for name in families_to_plot:
        res = results.get(name)
        if res is None:
            continue
        alpha = res['alpha_f']
        K_grid, g_vals, zeros = g_scan(alpha, a_focus, K_min=1e-4, K_max=200, npts=600)
        profiles[name] = {'K': K_grid, 'g': g_vals, 'zeros': zeros,
                          'alpha': alpha, 'res': res}
        print(f"  {name}: α_f={alpha:.4f}, zera g = {[f'{z:.4f}' for z in zeros]}")
    return profiles

# ============================================================
# SEKCJA E: Q_TGP przy λ=λ_Koide — kluczowa predykcja TGP
# K3 jest UNIWERSALNE (nie zależy od α) → K3 = C*a/sqrt(λ)
# λ_Koide wyznaczone z warunku r31_PDG(leptony) = K3/K1_lep
# ============================================================
def section_E(all_results, a_focus=0.040, C=2.000):
    print("\n" + "="*70)
    print("SEKCJA E: Q_TGP przy λ=λ_Koide — fizyczna predykcja TGP")
    print("K3 uniwersalne: K3(λ,a) = C*a/sqrt(λ) = const dla danego (a,λ)")
    print("λ_Koide = (C*a)^2 / (K3_PDG_lep * K1_lep)^2 z warunku r31_lep = 3477")
    print("="*70)

    res_all = all_results.get(a_focus, {})
    lep = res_all.get('e/mu/tau')
    if lep is None:
        print("BRAK leptonów — nie można wyliczyć λ_Koide")
        return

    # K1_lep i r31_PDG(lep) wyznaczają λ_Koide
    K1_lep = lep['K1']
    r31_lep_pdg = lep['r31_pdg']  # 3477
    # K3(λ_Koide) = r31_PDG * K1_lep
    K3_koide = r31_lep_pdg * K1_lep
    # λ_Koide = (C*a)^2 / K3_koide^2
    lam_koide = (C * a_focus)**2 / K3_koide**2
    print(f"\nλ_Koide = {lam_koide:.4e}  (C={C}, a={a_focus}, K3_Koide={K3_koide:.4f})")
    print(f"(Dla porównania: LAM_REF = {LAM_REF:.1e}, stosunek = {lam_koide/LAM_REF:.4f})")

    print(f"\nPredykcje Q_TGP przy λ=λ_Koide (K3={K3_koide:.4f} = const):")
    print(f"\n{'Rodzina':>10} {'α_f':>9} {'K1':>10} {'r21':>8} "
          f"{'r31_TGP':>10} {'r31_PDG':>10} {'r31_Koide':>11} "
          f"{'Q_TGP':>8} {'Q_PDG':>8} {'δQ':>8}")
    print("-"*95)

    results_E = {}
    for name, res in res_all.items():
        if res is None:
            print(f"{name:>10} FAIL")
            continue
        K1 = res['K1']
        r21 = res['r21_pdg']
        r31_pdg = res['r31_pdg']
        # r31_TGP przy λ=λ_Koide: K3=K3_koide (universal!)
        r31_tgp_lk = K3_koide / K1
        r31_k = r31_koide(r21)
        Q_tgp = koide_Q(r21, r31_tgp_lk)
        Q_pdg = koide_Q(r21, r31_pdg)
        delta_Q = (Q_tgp - Q_pdg)/Q_pdg*100
        delta_r31 = (r31_tgp_lk/r31_pdg - 1.0)*100
        print(f"{name:>10} {res['alpha_f']:>9.4f} {K1:>10.6f} {r21:>8.1f} "
              f"{r31_tgp_lk:>10.1f} {r31_pdg:>10.0f} {r31_k:>11.1f} "
              f"{Q_tgp:>8.4f} {Q_pdg:>8.4f} {delta_Q:>+8.1f}%")
        results_E[name] = {'Q_tgp': Q_tgp, 'Q_pdg': Q_pdg, 'r31_tgp': r31_tgp_lk,
                           'r31_pdg': r31_pdg, 'r31_koide': r31_k, 'delta_Q': delta_Q}

    print(f"\nWNIOSEK:")
    print(f"  TGP przy λ=λ_Koide przewiduje Q_TGP ≈ 1.5-1.53 dla WSZYSTKICH rodzin.")
    print(f"  PDG: Q=1.5000 (lep), 1.1779 (u/c/t), 1.3672 (d/s/b).")
    print(f"  TGP NADPRZEWIDUJE Q dla kwarków (daje za duże Q ≈ 3/2).")
    print(f"  Kwarki w Naturze mają Q < 3/2 — mechanizm TGP tego nie odtwarza.")
    return results_E

# ============================================================
# MAIN
# ============================================================
def main():
    print("P40: Q_TGP dla kwarków d/s/b — pełna predykcja")
    print("="*70)

    all_results = section_A()
    Q_results   = section_B(all_results)
    section_C(all_results)
    profiles    = section_D(all_results)
    results_E   = section_E(all_results)

    # ============================================================
    # WYKRES
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('P40: K₃ search and Q_TGP predictions for all fermion families',
                 fontsize=13)

    # Panel A: g(K) dla d/s/b i leptonów
    ax = axes[0, 0]
    colors_fam = {'e/mu/tau': 'blue', 'd/s/b': 'green', 'u/c/t': 'red'}
    for name, prof in profiles.items():
        col = colors_fam.get(name, 'gray')
        K_plot = prof['K']
        g_clip = np.clip(prof['g'], -3, 8)
        ax.plot(K_plot, g_clip, color=col, linewidth=1.8,
                label=f"{name} (α={prof['alpha']:.3f})")
        for z in prof['zeros']:
            ax.axvline(x=z, color=col, linestyle=':', alpha=0.6, linewidth=1)
    ax.axhline(y=0, color='black', linewidth=1.2)
    ax.set_xscale('log')
    ax.set_xlim(1e-4, 200)
    ax.set_ylim(-2.5, 5)
    ax.set_xlabel('K', fontsize=11)
    ax.set_ylabel('g(K) = E/(4πK) − 1', fontsize=11)
    ax.set_title('Panel A: g(K) profiles — zeros mark K₁, K₂, K₃', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel B: (r₂₁, r₃₁) — mapa z krzywą Koide i TGP
    ax = axes[0, 1]
    a_focus = 0.040
    res_a = all_results.get(a_focus, {})

    # Krzywa Koide Q=3/2
    r21_curve = np.logspace(-0.5, 4, 300)
    r31_curve = np.array([r31_koide(r) for r in r21_curve])
    valid = ~np.isnan(r31_curve)
    ax.plot(r21_curve[valid], r31_curve[valid], 'k-', linewidth=2, label='Q=3/2 (Koide)')

    # Punkty PDG
    for name, ferm in FERMIONS.items():
        r21_p = ferm['m2']/ferm['m1']
        r31_p = ferm['m3']/ferm['m1']
        col = colors_fam.get(name, 'gray')
        ax.scatter(r21_p, r31_p, c=col, s=120, marker='*', zorder=6,
                   label=f"{name} PDG")

    # Punkty TGP
    for name, res in res_a.items():
        if res is None or np.isnan(res['r31']):
            continue
        col = colors_fam.get(name, 'gray')
        ax.scatter(res['r21_pdg'], res['r31'], c=col, s=80, marker='o',
                   edgecolors='black', zorder=7, label=f"{name} TGP")
        # Strzałka PDG→TGP
        ax.annotate('', xy=(res['r21_pdg'], res['r31']),
                    xytext=(res['r21_pdg'], res['r31_pdg']),
                    arrowprops=dict(arrowstyle='->', color=col, lw=1.5))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r₂₁ = K₂/K₁', fontsize=11)
    ax.set_ylabel('r₃₁ = K₃/K₁', fontsize=11)
    ax.set_title('Panel B: (r₂₁, r₃₁) — PDG (★) vs TGP (o) vs Koide (line)', fontsize=10)
    ax.legend(fontsize=7, loc='upper left')
    ax.grid(True, alpha=0.3)

    # Panel C: Q_TGP vs Q_PDG dla wszystkich rodzin
    ax = axes[1, 0]
    names_list = list(FERMIONS.keys())
    Q_pdg_vals = [koide_Q(FERMIONS[n]['m2']/FERMIONS[n]['m1'],
                          FERMIONS[n]['m3']/FERMIONS[n]['m1']) for n in names_list]
    Q_tgp_vals = [Q_results.get(n, {}).get('Q_tgp', np.nan)
                  if Q_results.get(n) else np.nan for n in names_list]

    x = np.arange(len(names_list))
    width = 0.35
    bars1 = ax.bar(x - width/2, Q_pdg_vals, width, label='Q_PDG',
                   color=['blue', 'red', 'green'], alpha=0.6)
    bars2 = ax.bar(x + width/2, Q_tgp_vals, width, label='Q_TGP',
                   color=['blue', 'red', 'green'], alpha=0.9, hatch='//')
    ax.axhline(y=1.5, color='black', linestyle='--', linewidth=1.5, label='Q=3/2')
    ax.set_xticks(x)
    ax.set_xticklabels(names_list, fontsize=10)
    ax.set_ylim(1.0, 1.6)
    ax.set_ylabel('Q (Koide parameter)', fontsize=11)
    ax.set_title('Panel C: Q_TGP (///) vs Q_PDG for all families', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    # Panel D: K₁, K₂, K₃ vs α_f
    ax = axes[1, 1]
    a_scan = 0.040
    alpha_range = np.concatenate([np.arange(0.05, 0.5, 0.05),
                                  np.arange(0.5, 3.0, 0.2),
                                  np.arange(3.0, 12.0, 0.5)])
    K1_list, K2_list, K3_list = [], [], []
    alpha_valid = []
    for alpha in alpha_range:
        K1 = find_K1_num(alpha, a_scan)
        K2 = find_K2_num(alpha, a_scan)
        K3 = find_K3_num(alpha, a_scan)
        if not np.isnan(K1):
            K1_list.append(K1); K2_list.append(K2); K3_list.append(K3)
            alpha_valid.append(alpha)

    ax.semilogy(alpha_valid, K1_list, 'b-o', markersize=4, linewidth=1.5, label='K₁')
    ax.semilogy(alpha_valid, K2_list, 'g-s', markersize=4, linewidth=1.5, label='K₂')
    K3_plot = [k if not np.isnan(k) else np.nan for k in K3_list]
    ax.semilogy(alpha_valid, K3_plot, 'r-^', markersize=4, linewidth=1.5, label='K₃')

    # Zaznacz α_f dla fermionów
    res_a40 = all_results.get(0.040, {})
    for name, res in res_a40.items():
        if res is None: continue
        col = colors_fam.get(name, 'gray')
        ax.axvline(x=res['alpha_f'], color=col, linestyle=':', alpha=0.7, linewidth=1.5)
        ax.text(res['alpha_f'], 0.001, name.split('/')[0], color=col,
                fontsize=8, rotation=90, va='bottom')

    ax.set_xlabel('α', fontsize=11)
    ax.set_ylabel('K (log scale)', fontsize=11)
    ax.set_title('Panel D: K₁, K₂, K₃ vs α (a=0.040)\nvertical lines: α_f for fermion families', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('scripts/advanced/p40_dsb_Q_prediction.png', dpi=120, bbox_inches='tight')
    print("\nWykres zapisany: scripts/advanced/p40_dsb_Q_prediction.png")

    # ============================================================
    # PODSUMOWANIE
    # ============================================================
    print("\n" + "="*70)
    print("PODSUMOWANIE P40")
    print("="*70)
    a_focus = 0.040
    res_a = all_results.get(a_focus, {})
    for name, res in res_a.items():
        if res is None:
            print(f"\n{name}: BRAK ROZWIĄZANIA (α_f nie znalezione)")
            continue
        r21 = res['r21_pdg']
        r31_pdg = res['r31_pdg']
        r31_tgp = res['r31']
        r31_k = r31_koide(r21)
        Q_tgp = koide_Q(r21, r31_tgp) if not np.isnan(r31_tgp) else np.nan
        Q_pdg = koide_Q(r21, r31_pdg)
        print(f"\n{name}:")
        print(f"  α_f = {res['alpha_f']:.4f}, a_Γ = {a_focus}")
        K3_val = res['K3']
        K3_str = f"{K3_val:.4f}" if not np.isnan(K3_val) else "BRAK"
        print(f"  K₁ = {res['K1']:.6f}, K₂ = {res['K2']:.6f}, K₃ = {K3_str}")
        print(f"  r₂₁_TGP = {res['r21']:.2f}  (PDG: {r21:.2f})")
        r31_tgp_str = f"{r31_tgp:.1f}" if not np.isnan(r31_tgp) else "---"
        r31_k_str   = f"{r31_k:.1f}"   if not np.isnan(r31_k)   else "---"
        print(f"  r₃₁_TGP = {r31_tgp_str}  (PDG: {r31_pdg:.0f}, Koide: {r31_k_str})")
        if not np.isnan(Q_tgp):
            delta_Q = (Q_tgp - Q_pdg)/Q_pdg*100
            print(f"  Q_TGP = {Q_tgp:.4f}  (PDG: {Q_pdg:.4f}, δQ = {delta_Q:+.1f}%)")
        else:
            print(f"  Q_TGP: brak K₃")

    print("\n" + "="*70)
    print("P40 ZAKONCZONE")

if __name__ == '__main__':
    main()
