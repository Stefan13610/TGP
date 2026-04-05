"""
p21_psi_ratio.py
================
!! PRZESTARZALY (DEPRECATED) !!
   Skrypt uzywa siatki ROWNOMIERNEJ, ktora daje blad ~14% dla K*3.
   Uzywaj p21b_psi_ratio_log.py (siatka logarytmiczna -- poprawna).
   Zachowany tylko dla historii.
!! END DEPRECATED !!

CEL: Wyjaśnienie obserwacji ψ_core/ψ_new ≈ 1.93–1.97 = const dla WSZYSTKICH
     rozwiązań w rodzinie optymalnej (α, a_Γ, λ*).

HIPOTEZA (z ANALIZA_SPRZEZONY.md):
  ψ_core(M₃) / ψ_new ≈ e/√2 = 1.9221...

PLAN:
  Krok 1: Wyznacz K*₁, K*₂, K*₃ dla 4 members rodziny optymalnej
           (te same parametry które dają r₂₁=207, r₃₁=3477).
  Krok 2: Oblicz ψ_core(Mn) = 1 + K*n / a_Γ * exp(-a_Γ) dla każdego n
  Krok 3: Oblicz ψ_new = 1/sqrt(λ*) dla każdej family member
  Krok 4: Sprawdź czy ψ_core(M₃)/ψ_new = const ≈ e/√2
  Krok 5: Derivacja analityczna — z warunku samospójności g(K*₃)=0 + ψ_new=1/√λ*
  Krok 6: Skan "pseudo-rodzin" — (α, a_Γ, λ) z r₂₁=207 ale różne r₃₁:
           czy ratio jest stałe też tam?

PARAMETRY: Ta sama funkcja energy_yukawa co w p18.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAMMA     = 1.0
PHI0_BASE = 1.0

def energy_yukawa(K, Phi0, alpha, a_gam, lam, n_eval=3000):
    """Energia solitonu Yukawa w tle Phi0 (identyczna jak p18)."""
    msp   = np.sqrt(max(GAMMA * PHI0_BASE / Phi0, 1e-9))
    r_max = max(80.0 / msp, 20.0)
    r     = np.linspace(a_gam, r_max, n_eval)

    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / phi) * r**2, r)

    psi = phi / Phi0
    V1  = GAMMA/3 - GAMMA/4
    Ep  = 4*np.pi * Phi0**3 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4
         + lam/6*(psi - 1.0)**6 - V1) * r**2, r)

    return Ek + Ep


def g_func(M, Phi0, alpha, a_gam, lam):
    if M < 1e-12: return 0.0
    K = M / (4.0 * np.pi * Phi0)
    E = energy_yukawa(K, Phi0, alpha=alpha, a_gam=a_gam, lam=lam)
    return E - M


def find_M_stars(Phi0, alpha, a_gam, lam, N=400):
    """Znajdź wszystkie zera g(M;Phi0) – skan wieloskalowy + brentq."""
    M_arr = np.unique(np.concatenate([
        np.linspace(1e-4,   2.0,   int(N*0.15)),
        np.linspace(2.0,    100,   int(N*0.15)),
        np.linspace(100,    5000,  int(N*0.25)),
        np.linspace(5000,   50000, int(N*0.25)),
        np.linspace(50000, 5e5,    int(N*0.20)),
    ]))
    g_arr = np.array([g_func(M, Phi0, alpha, a_gam, lam) for M in M_arr])

    roots = []
    for i in range(len(M_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)):
            continue
        if gi * gj < 0:
            try:
                M_z = brentq(
                    lambda M: g_func(M, Phi0, alpha, a_gam, lam),
                    M_arr[i], M_arr[i+1],
                    xtol=1e-4, rtol=1e-6, maxiter=60)
                # Sprawdź czy to nie duplikat
                if not roots or abs(M_z - roots[-1]) > M_z * 0.01:
                    roots.append(M_z)
            except Exception:
                pass

    return sorted(roots)


# ============================================================
# Rodzina optymalna: (α, a_Γ, λ*) dające r₂₁=207, r₃₁=3477
# ============================================================
FAMILY = [
    (5.9148, 0.025, 2.883e-6),
    (6.8675, 0.030, 3.743e-6),
    (7.7449, 0.035, 4.618e-6),
    (8.5616, 0.040, 5.501e-6),
]

print("P21: Analiza stosunku ψ_core(M₃)/ψ_new w rodzinie optymalnej")
print("=" * 70)
print(f"  Hipoteza: ψ_core(M₃)/ψ_new ≈ e/√2 = {np.e/np.sqrt(2):.5f}")
print()

results = []

for alpha, a_gam, lam in FAMILY:
    print(f"\n  --- α={alpha}, a_Γ={a_gam}, λ*={lam:.3e} ---")
    PHI0 = 1.0
    msp  = np.sqrt(GAMMA * PHI0_BASE / PHI0)  # = 1.0

    # Znajdź M*₁, M*₂, M*₃
    roots = find_M_stars(PHI0, alpha, a_gam, lam)
    print(f"  Znaleziono {len(roots)} zer g(M):")
    for idx, M in enumerate(roots):
        K = M / (4*np.pi*PHI0)
        psi_c = 1.0 + K / a_gam * np.exp(-msp * a_gam)
        print(f"    M*_{idx+1} = {M:.4f},  K = {K:.5f},  ψ_core = {psi_c:.3f}")

    if len(roots) < 3:
        print("  UWAGA: nie znaleziono 3 zer — poszerzam skan...")
        # Spróbuj z większym zakresem i gęstszą siatką
        roots2 = find_M_stars(PHI0, alpha, a_gam, lam, N=800)
        if len(roots2) > len(roots):
            roots = roots2
            print(f"  Po poszerzeniu: {len(roots)} zer")
            for idx, M in enumerate(roots):
                K = M / (4*np.pi*PHI0)
                psi_c = 1.0 + K / a_gam * np.exp(-msp * a_gam)
                print(f"    M*_{idx+1} = {M:.4f},  K = {K:.5f},  ψ_core = {psi_c:.3f}")

    if len(roots) >= 3:
        M1, M2, M3 = roots[0], roots[1], roots[2]
        K1 = M1 / (4*np.pi*PHI0)
        K2 = M2 / (4*np.pi*PHI0)
        K3 = M3 / (4*np.pi*PHI0)
        psi_core_1 = 1.0 + K1 / a_gam * np.exp(-msp * a_gam)
        psi_core_2 = 1.0 + K2 / a_gam * np.exp(-msp * a_gam)
        psi_core_3 = 1.0 + K3 / a_gam * np.exp(-msp * a_gam)
        psi_new    = 1.0 / np.sqrt(lam)
        r21 = M2/M1
        r31 = M3/M1
        ratio3 = psi_core_3 / psi_new
        ratio2 = psi_core_2 / psi_new
        ratio1 = psi_core_1 / psi_new
        print(f"  r₂₁ = {r21:.2f},  r₃₁ = {r31:.1f}")
        print(f"  ψ_new = 1/√λ* = {psi_new:.2f}")
        print(f"  ψ_core(M₁)/ψ_new = {ratio1:.5f}")
        print(f"  ψ_core(M₂)/ψ_new = {ratio2:.5f}")
        print(f"  ψ_core(M₃)/ψ_new = {ratio3:.5f}  (hipoteza: {np.e/np.sqrt(2):.5f})")
        results.append({
            'alpha': alpha, 'a_gam': a_gam, 'lam': lam,
            'M1': M1, 'M2': M2, 'M3': M3,
            'K1': K1, 'K2': K2, 'K3': K3,
            'psi_core_1': psi_core_1, 'psi_core_2': psi_core_2,
            'psi_core_3': psi_core_3, 'psi_new': psi_new,
            'ratio3': ratio3, 'r21': r21, 'r31': r31,
        })
    elif len(roots) == 2:
        M1, M2 = roots[0], roots[1]
        K1 = M1/(4*np.pi); K2 = M2/(4*np.pi)
        psi_core_3 = np.nan; psi_new = 1.0/np.sqrt(lam)
        print(f"  Tylko 2 zera — brak M₃ (g wciąż ujemne po M₂?)")
        print(f"  r₂₁ = {M2/M1:.2f}")
        results.append({
            'alpha': alpha, 'a_gam': a_gam, 'lam': lam,
            'M1': M1, 'M2': M2, 'M3': np.nan,
            'K1': K1, 'K2': K2, 'K3': np.nan,
            'psi_core_3': np.nan, 'psi_new': psi_new, 'ratio3': np.nan,
            'r21': M2/M1, 'r31': np.nan,
        })
    else:
        print(f"  Znaleziono < 2 zera — problem numeryczny.")
        results.append({'alpha': alpha, 'a_gam': a_gam, 'lam': lam,
                        'ratio3': np.nan, 'r21': np.nan, 'r31': np.nan})


# ============================================================
# KROK 5: Próba analitycznego wyprowadzenia ratio
# ============================================================
print()
print("=" * 70)
print("KROK 5: Analityczne wyjaśnienie ψ_core(M₃)/ψ_new = e/√2 ?")
print("=" * 70)
print()
print("Dla Yukawa przy Phi0=1, m_sp=1, a_Γ małe:")
print("  ψ_core(M₃) ≈ K*₃ / a_Γ × exp(-a_Γ)")
print("  ψ_new = 1/√λ*")
print("  ratio = K*₃ × √λ* × exp(-a_Γ) / a_Γ")
print()

# Dla każdej family member: oblicz analityczne składniki
for r in results:
    if np.isnan(r.get('ratio3', np.nan)): continue
    K3 = r['K3']; a_gam = r['a_gam']; lam = r['lam']
    anal = K3 * np.sqrt(lam) * np.exp(-a_gam) / a_gam
    print(f"  α={r['alpha']}, a_Γ={a_gam}: K₃={K3:.4f}, "
          f"K₃×√λ×exp(-a_Γ)/a_Γ = {anal:.5f}"
          f"  (vs ratio = {r['ratio3']:.5f})")

print()
print("Jeśli K₃×√λ*/a_Γ × exp(-a_Γ) ≈ e/√2, to:")
print("  K₃ ≈ (e/√2) × a_Γ/√λ* × exp(a_Γ) = (e/√2) × a_Γ × ψ_new × exp(a_Γ)")
print()
print("Warunek g(K₃)=0 daje dodatkowe ograniczenie — to wyznacza K₃.")
print()

# Sprawdź "proste" skalowanie K₃ ~ a_Γ × ψ_new
for r in results:
    if np.isnan(r.get('ratio3', np.nan)): continue
    K3 = r['K3']; a_gam = r['a_gam']; psi_new = r['psi_new']
    ratio_kp = K3 / (a_gam * psi_new)
    print(f"  α={r['alpha']}: K₃/(a_Γ×ψ_new) = {K3:.4f}/({a_gam}×{psi_new:.1f}) = {ratio_kp:.5f}")

# ============================================================
# KROK 6: Sprawdź ratio dla "pseudo-rodziny"
# (ustalone α=8.56, a_Γ=0.040, λ* zmienne — dając r₃₁≠3477)
# ============================================================
print()
print("=" * 70)
print("KROK 6: Czy ratio jest stałe dla różnych λ* (pseudo-rodzina)?")
print("=" * 70)
alpha_ref = 8.5616
a_gam_ref = 0.040
lam_scan  = [1e-6, 2e-6, 3e-6, 5.501e-6, 8e-6, 1.2e-5]

pseudo_results = []
print(f"  Stałe: α={alpha_ref}, a_Γ={a_gam_ref}, Phi0=1")
print(f"  {'λ*':>10} {'ψ_new':>8} {'M₃':>10} {'K₃':>8} {'ratio3':>8} {'r₃₁':>7}")
print("  " + "-"*56)

for lam in lam_scan:
    roots = find_M_stars(1.0, alpha_ref, a_gam_ref, lam, N=500)
    psi_new = 1.0/np.sqrt(lam)
    if len(roots) >= 3:
        M3 = roots[2]; K3 = M3/(4*np.pi)
        psi_c3 = 1.0 + K3/a_gam_ref * np.exp(-a_gam_ref)
        ratio3 = psi_c3 / psi_new
        r31 = roots[2]/roots[0] if roots else np.nan
        print(f"  {lam:>10.2e} {psi_new:>8.1f} {M3:>10.2f} {K3:>8.4f} {ratio3:>8.4f} {r31:>7.0f}")
        pseudo_results.append((lam, psi_new, ratio3, r31))
    else:
        print(f"  {lam:>10.2e} {psi_new:>8.1f} {'<3 zera':>30}")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE p21: ψ_core(M₃)/ψ_new")
print("=" * 70)
print()
print(f"  {'α':>6} {'a_Γ':>6} {'λ*':>10} {'K*₃':>8} {'ψ_new':>7} "
      f"{'ratio':>7} {'r₂₁':>6} {'r₃₁':>7}")
print("  " + "-"*64)
for r in results:
    if np.isnan(r.get('ratio3', np.nan)):
        continue
    ratio3_str = f"{r['ratio3']:.4f}"
    r21_str    = f"{r['r21']:.1f}"
    r31_str    = f"{r['r31']:.0f}"
    print(f"  {r['alpha']:>6.4f} {r['a_gam']:>6.3f} {r['lam']:>10.3e} "
          f"{r['K3']:>8.4f} {r['psi_new']:>7.1f} "
          f"{ratio3_str:>7} {r21_str:>6} {r31_str:>7}")

ratios_valid = [r['ratio3'] for r in results if not np.isnan(r.get('ratio3', np.nan))]
if ratios_valid:
    mean_r = np.mean(ratios_valid)
    std_r  = np.std(ratios_valid)
    print()
    print(f"  Średni stosunek: {mean_r:.5f} ± {std_r:.5f}")
    print(f"  e/√2           : {np.e/np.sqrt(2):.5f}")
    print(f"  Odchylenie od e/√2: {abs(mean_r - np.e/np.sqrt(2))/np.e*np.sqrt(2)*100:.2f}%")
    print()
    if std_r / mean_r < 0.015:
        print("  ✓ RATIO JEST STAŁE w rodzinie (std/mean < 1.5%)")
    else:
        print("  ✗ Ratio zmienne w rodzinie (std/mean > 1.5%)")

# ============================================================
# WYKRES
# ============================================================
print("\nTworze wykres...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P21: ψ_core(M₃)/ψ_new — stałość w rodzinie optymalnej', fontsize=11)

ax = axes[0]
alphas_v = [r['alpha'] for r in results if not np.isnan(r.get('ratio3', np.nan))]
ratios_v = [r['ratio3'] for r in results if not np.isnan(r.get('ratio3', np.nan))]
if alphas_v:
    ax.plot(alphas_v, ratios_v, 'bo-', ms=8, lw=2, label='ψ_core(M₃)/ψ_new')
    ax.axhline(np.e/np.sqrt(2), color='r', lw=1.5, ls='--', label=f'e/√2={np.e/np.sqrt(2):.4f}')
ax.set_xlabel('α'); ax.set_ylabel('ψ_core/ψ_new')
ax.set_title('Stosunek vs α (rodzina optymalna)'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
# g(M) dla najlepszych parametrów
alpha_b, a_gam_b, lam_b = FAMILY[-1]
M_plot = np.concatenate([
    np.linspace(1e-3, 2, 30),
    np.logspace(0.3, 2.5, 40),
    np.logspace(2.5, 5, 40),
])
g_plot = [g_func(M, 1.0, alpha_b, a_gam_b, lam_b) for M in M_plot]
g_clip = np.clip(g_plot, -5e4, 5e4)
ax.semilogx(M_plot, g_clip, 'b-', lw=1.5)
ax.axhline(0, color='k', lw=1, ls='--')
# Zaznacz M*₁, M*₂, M*₃
r_last = [r for r in results if r['alpha']==alpha_b]
if r_last and not np.isnan(r_last[0].get('M3', np.nan)):
    for Mk, lab in [(r_last[0]['M1'],'M*₁'),(r_last[0]['M2'],'M*₂'),(r_last[0]['M3'],'M*₃')]:
        ax.axvline(Mk, color='green', lw=1, ls=':')
        ax.text(Mk, 1e3, lab, fontsize=7, ha='center')
ax.set_xlabel('M'); ax.set_ylabel('g(M)')
ax.set_title(f'g(M) dla α={alpha_b}, a_Γ={a_gam_b}', fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[2]
# Pseudo-rodzina: ratio vs λ*
if pseudo_results:
    lams_p, pns_p, rats_p, r31s_p = zip(*pseudo_results)
    ax.semilogx(lams_p, rats_p, 'rs-', ms=7, lw=2, label='pseudo-rodzina')
    ax.axhline(np.e/np.sqrt(2), color='k', lw=1.5, ls='--', label=f'e/√2')
ax.set_xlabel('λ*'); ax.set_ylabel('ψ_core(M₃)/ψ_new')
ax.set_title('Ratio vs λ* (pseudo-rodzina, α=8.56 fixed)', fontsize=9)
ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

# ============================================================
# ANALIZA ANALITYCZNA (przy założeniu, że ratio jest stałe)
# ============================================================
print()
print("=" * 70)
print("ANALIZA ANALITYCZNA: Warunek na ratio = C = const")
print("=" * 70)
print()
print("Z definicji:")
print("  ψ_core(M₃) = 1 + K*₃/a_Γ × exp(-a_Γ)  [m_sp=1, Phi0=1]")
print("  ψ_new = 1/√λ*")
print("  ratio = ψ_core/ψ_new = [1 + K*₃/a_Γ exp(-a_Γ)] × √λ*")
print()
print("Dla K*₃/a_Γ >> 1 (co jest prawdą dla M₃ ~ 10^4-10^5 >> a_Γ ~ 0.025-0.04):")
print("  ratio ≈ K*₃ × √λ* × exp(-a_Γ) / a_Γ")
print()
print("Z warunku samospójności g(K*₃)=0: E(K*₃) = 4π K*₃")
print("Dla dużych K₃ zdominowanych przez potencjał V_mod ~ -ψ⁴/4 + λψ⁶/6:")
print("  ψ(r) ~ K*₃/r × exp(-r) >> ψ_new = 1/√λ  w rdzeniu")
print("  Energia: E ≈ (4π K*₃²) × [kinetyk] + (4π) × [potenc. z λ]")
print()
print("Warunek E=4πK*₃ przy ψ_core ~ C/√λ gdzie C=ratio:")
print("  Kluczowy term energii: λ/6 × ψ^6 całkuje się do ~λ/6 × ψ_core^6 × a_Γ³")
print("  = λ/6 × (C/√λ)^6 × a_Γ³ = C^6/(6λ²) × a_Γ³")
print("  Musi być porównywalne z 4πK*₃ ≈ 4π × C × a_Γ/√λ × exp(a_Γ)")
print()
print("WNIOSEK: Ratio ≈ e/√2 wynika z balansu między członem λψ⁶")
print("  (stabilizującym) a krytycznym gradientem solitonu na granicy a_Γ.")
print("  Dokładna wartość wymaga analizy pełnego warunku g(K*₃)=0.")
