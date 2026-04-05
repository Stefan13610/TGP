"""
p21b_psi_ratio_log.py
=====================
KOREKTA p21: uzywamy siatki LOGARYTMICZNEJ (jak v2_bisekcja_alpha),
ktora daje poprawne K* = {K1=0.009820, K2=2.033, K3=34.14} dla optymalnych param.

Problem z p21: siatka rownierna nie oddaje ostrego pieku energii przy r=a_Gamma
dla duzych K (K3~30-34), dajac blad ~15% w K3.

PLAN:
  Krok 1: Zweryfikuj K* dla optymalnych parametrow (alfa=8.5616, a_Gamma=0.040, lam=5.501e-6)
           -> porownaj z K1=0.009820, K2=2.032726, K3=34.14423 (z v2_weryfikacja_finalna)
  Krok 2: Oblicz psi_core/psi_new dla calej rodziny optymalnej
           -> zweryfikuj ze ratio in [1.93, 1.97]
  Krok 3: Analiza psi_core/psi_new w pseudo-rodzinie (alfa=8.56, zmienne lam)
           -> czy ratio stale dla stalego a_Gamma?
  Krok 4: Wyprowadzenie analityczne ratio z warunku samospojnosci
  Krok 5: Sprawdz hipoteze e/sqrt(2) = 1.9221

KLUCZOWE ODKRYCIE z p21: siatka rownierna (p18/p21) daje K3 o 14% za male,
  zanizona psi_core, ratio=1.57-1.68 zamiast 1.93-1.97.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAMMA = 1.0
R_MAX = 60.0

# ============================================================
# ENERGIA -- identyczna jak v2_bisekcja_alpha (siatka LOG)
# ============================================================
def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_log(K, alpha, a_gam, lam, N=1500):
    """Energia solitonu -- siatka LOGARYTMICZNA od a_gam do R_MAX."""
    t   = np.linspace(0, 1, N)
    r   = a_gam * (R_MAX/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi= K*np.exp(-r)*(-r - 1.0)/r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep  = 4*np.pi*np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam, N=1500):
    return energy_log(K, alpha, a_gam, lam, N)/K - 4*np.pi

def find_K_stars(alpha, a_gam, lam, K_max=150.0, N_scan=800, N_en=1500):
    """Znajdz wszystkie zera g(K) = 0."""
    K_arr = np.unique(np.concatenate([
        np.linspace(1e-4, 0.5, 150),
        np.linspace(0.5,  10.0, 250),
        np.linspace(10.0, K_max, 400),
    ]))
    f_arr = np.array([g_func(K, alpha, a_gam, lam, N_en) for K in K_arr])
    roots = []
    for i in range(len(K_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                K_z = brentq(lambda K: g_func(K, alpha, a_gam, lam, N_en),
                             K_arr[i], K_arr[i+1], xtol=1e-7, maxiter=60)
                if not roots or abs(K_z - roots[-1]) > K_z*0.01:
                    roots.append(K_z)
            except Exception:
                pass
    return sorted(roots)

# ============================================================
# KROK 1: Weryfikacja dla najlepszego rozwiazania
# ============================================================
print("P21b: Analiza psi_core(M3)/psi_new -- SIATKA LOGARYTMICZNA")
print("=" * 70)
print()
print("KROK 1: Weryfikacja K* dla alfa=8.5616, a_Gamma=0.040, lam=5.501e-6")
print("-" * 70)

ALPHA_REF = 8.5616
A_GAM_REF = 0.040
LAM_REF   = 5.501e-6

roots_ref = find_K_stars(ALPHA_REF, A_GAM_REF, LAM_REF)
print(f"  Znaleziono {len(roots_ref)} zer K:")
for i, K in enumerate(roots_ref[:5]):
    psi_c = 1.0 + K/A_GAM_REF * np.exp(-A_GAM_REF)
    print(f"    K*_{i+1} = {K:.6f},  psi_core = {psi_c:.3f}")

if len(roots_ref) >= 3:
    K1, K2, K3 = roots_ref[0], roots_ref[1], roots_ref[2]
    r21 = K2/K1; r31 = K3/K1
    psi_new = 1.0/np.sqrt(LAM_REF)
    psi_c3 = 1.0 + K3/A_GAM_REF * np.exp(-A_GAM_REF)
    ratio3  = psi_c3/psi_new
    print()
    print(f"  K1 = {K1:.6f}  (oczekiwane: 0.009820)")
    print(f"  K2 = {K2:.6f}  (oczekiwane: 2.032726)")
    print(f"  K3 = {K3:.6f}  (oczekiwane: 34.14423)")
    print(f"  r21 = {r21:.3f}  (oczekiwane: ~207)")
    print(f"  r31 = {r31:.1f}  (oczekiwane: 3477)")
    print()
    print(f"  psi_new = 1/sqrt(lam*) = {psi_new:.2f}")
    print(f"  psi_core(K3) = {psi_c3:.2f}")
    print(f"  ratio = psi_core/psi_new = {ratio3:.4f}")
    print(f"  hipoteza e/sqrt(2) = {np.e/np.sqrt(2):.4f}")
    print(f"  odchylenie: {abs(ratio3-np.e/np.sqrt(2))/np.e*np.sqrt(2)*100:.3f}%")

# ============================================================
# KROK 2: Pelna rodzina optymalna (4 members)
# ============================================================
print()
print("=" * 70)
print("KROK 2: Rodzina optymalna -- psi_core/psi_new dla wszystkich")
print("=" * 70)

# Parametry z v2_bisekcja_alpha (tabela ANALIZA_SPRZEZONY linie 430-435)
FAMILY = [
    (5.9148, 0.025, 2.883e-6,  1158, 589),
    (6.8675, 0.030, 3.743e-6,  1009, 517),
    (7.7449, 0.035, 4.618e-6,   902, 465),
    (8.5616, 0.040, 5.501e-6,   821, 426),
]

results_fam = []
print(f"\n  {'alpha':>7} {'a_Gam':>6} {'lam':>10}  {'K1':>8} {'K2':>8} {'K3':>8}  "
      f"{'r21':>6} {'r31':>6}  {'psi_c3':>7} {'psi_new':>7} {'ratio':>7}")
print("  " + "-"*100)

for alpha, a_gam, lam, psi_c3_ref, psi_new_ref in FAMILY:
    roots = find_K_stars(alpha, a_gam, lam)
    psi_new = 1.0/np.sqrt(lam)
    if len(roots) >= 3:
        K1, K2, K3 = roots[0], roots[1], roots[2]
        r21 = K2/K1; r31 = K3/K1
        psi_c3 = 1.0 + K3/a_gam * np.exp(-a_gam)
        ratio3 = psi_c3/psi_new
        results_fam.append({
            'alpha': alpha, 'a_gam': a_gam, 'lam': lam,
            'K1': K1, 'K2': K2, 'K3': K3, 'r21': r21, 'r31': r31,
            'psi_c3': psi_c3, 'psi_new': psi_new, 'ratio3': ratio3,
            'psi_c3_ref': psi_c3_ref, 'psi_new_ref': psi_new_ref,
        })
        match_c = "OK" if abs(psi_c3-psi_c3_ref)/psi_c3_ref < 0.02 else "???"
        print(f"  {alpha:>7.4f} {a_gam:>6.3f} {lam:>10.3e}  "
              f"{K1:>8.5f} {K2:>8.4f} {K3:>8.4f}  "
              f"{r21:>6.1f} {r31:>6.0f}  "
              f"{psi_c3:>7.1f}({match_c}) {psi_new:>7.1f} {ratio3:>7.4f}")
    else:
        print(f"  {alpha:>7.4f} {a_gam:>6.3f} {lam:>10.3e}  <{len(roots)} zer>")

if results_fam:
    ratios = [r['ratio3'] for r in results_fam]
    print()
    print(f"  Ratios: {', '.join(f'{r:.4f}' for r in ratios)}")
    print(f"  Srednia: {np.mean(ratios):.4f},  std: {np.std(ratios):.4f}")
    print(f"  e/sqrt(2) = {np.e/np.sqrt(2):.4f}")
    print(f"  Zakres odchylen od e/sqrt(2): "
          f"{min(abs(r-np.e/np.sqrt(2)) for r in ratios):.4f} .. "
          f"{max(abs(r-np.e/np.sqrt(2)) for r in ratios):.4f}")

# ============================================================
# KROK 3: Pseudo-rodzina -- alfa=8.56 stale, zmienne lam
# ============================================================
print()
print("=" * 70)
print("KROK 3: Pseudo-rodzina (alfa=8.56, a_Gam=0.040, lam zmienny)")
print("=" * 70)

lam_scan = np.logspace(-6.0, -4.3, 12)
pseudo_res = []
print(f"  {'lam':>10} {'K3':>8} {'psi_new':>7} {'psi_c3':>8} {'ratio':>7} {'r31':>6}")
print("  " + "-"*55)
for lam in lam_scan:
    roots = find_K_stars(ALPHA_REF, A_GAM_REF, lam)
    if len(roots) >= 3:
        K1, K2, K3 = roots[0], roots[1], roots[2]
        psi_new = 1.0/np.sqrt(lam)
        psi_c3 = 1.0 + K3/A_GAM_REF * np.exp(-A_GAM_REF)
        ratio3 = psi_c3/psi_new
        r31 = K3/K1
        print(f"  {lam:>10.3e} {K3:>8.4f} {psi_new:>7.2f} {psi_c3:>8.2f} "
              f"{ratio3:>7.4f} {r31:>6.0f}")
        pseudo_res.append((lam, K3, psi_new, psi_c3, ratio3, r31))

# ============================================================
# KROK 4: Analityczne wyprowadzenie ratio
# ============================================================
print()
print("=" * 70)
print("KROK 4: Analiza analityczna psi_core(K3)/psi_new")
print("=" * 70)
print()
print("Definicje:")
print("  psi_core(K3) = 1 + K3/a_Gam * exp(-a_Gam)  [Yukawa, m_sp=1, Phi0=1]")
print("  psi_new = 1/sqrt(lam*)")
print("  ratio R = psi_core/psi_new = K3 * sqrt(lam*) * exp(-a_Gam) / a_Gam + sqrt(lam*)")
print()
print("Dla K3/a_Gam >> 1:")
print("  R ~ K3 * sqrt(lam*) / a_Gam  (dominujacy czlon)")
print()
print("Skalowanie K3 z lam* (z danych pseudo-rodziny):")
if pseudo_res:
    lams_p = np.array([x[0] for x in pseudo_res])
    K3s_p  = np.array([x[1] for x in pseudo_res])
    # Dopasowanie K3 ~ lam^beta
    log_l = np.log(lams_p); log_K = np.log(K3s_p)
    c = np.polyfit(log_l, log_K, 1)
    print(f"  log-log fit: K3 ~ lam^{c[0]:.4f} x exp({c[1]:.4f})")
    print(f"  Tj. K3 ~ lam^{c[0]:.4f}  (oczekiwane: -1/2 gdyby K3~1/sqrt(lam))")
    print()
    # Sprawdz K3 * sqrt(lam)
    prod = K3s_p * np.sqrt(lams_p)
    print(f"  K3 * sqrt(lam*): {[f'{p:.4f}' for p in prod]}")
    print(f"  Srednia: {np.mean(prod):.5f}, std: {np.std(prod):.5f}")
    if np.std(prod)/np.mean(prod) < 0.01:
        print(f"  -> K3 * sqrt(lam*) = const = {np.mean(prod):.5f}  (< 1% std)")
    print()
    # Wnioskuj ratio
    C_prod = np.mean(prod)
    R_anal = C_prod * np.exp(-A_GAM_REF) / A_GAM_REF
    print(f"  Analityczne ratio = K3*sqrt(lam)*exp(-a)/a = {R_anal:.5f}")
    print(f"  vs pomierzone ratio: {np.mean([x[4] for x in pseudo_res]):.5f}")

print()
print("KLUCZOWE PYTANIE: Skad pochodzi K3*sqrt(lam*) = const?")
print()
print("Warunek g(K3) = E(K3)/(4pi*K3) - 1 = 0 => E(K3) = 4pi*K3")
print()
print("Dla K3 duzego (K3/a_Gam >> 1, psi_core >> psi_new >> 1):")
print("Dominujacy wklad energii: czlon lam*phi^6/6 przy r ~ a_Gam")
print("  E_lambda ~ 4pi * integral_a_Gam^infty (lam/6 * phi^6) * r^2 dr")
print("  phi(r) ~ K3/r * exp(-r), phi^6 ~ K3^6/r^6 * exp(-6r)")
print("  Integral ~ K3^6 * integral (exp(-6r)/r^4) dr from a_Gam")
print("  ~ K3^6 * exp(-6*a_Gam) / (5*a_Gam^5)  [dominant near lower bound]")
print()
print("Warunek E_lambda ~ 4pi*K3:")
print("  lam * K3^6 * exp(-6a) / (5*a^5) ~ K3")
print("  K3^5 ~ 5*a^5 * exp(6a) / lam")
print("  K3 ~ (5/lam)^(1/5) * a * exp(6a/5)")
print("  K3 * sqrt(lam) ~ sqrt(lam) * (5/lam)^(1/5) * a * exp(6a/5)")
print("                   = 5^(1/5) * lam^(1/2-1/5) * a * exp(6a/5)")
print("                   = 5^(1/5) * lam^(3/10) * a * exp(6a/5)")
print()
print("Sprawdzenie: K3 ~ lam^(-1/5)?")
if pseudo_res:
    print(f"  Dopasowanie log-log: K3 ~ lam^{c[0]:.4f}  (teoria: -1/5 = -0.200)")
    print(f"  Odchylenie: {abs(c[0]+0.2)/0.2*100:.1f}%")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE p21b")
print("=" * 70)
print()
e_sqrt2 = np.e/np.sqrt(2)
if results_fam:
    print(f"  Rodzina optymalna psi_core(M3)/psi_new:")
    print(f"    Zakres: [{min(r['ratio3'] for r in results_fam):.4f}, "
          f"{max(r['ratio3'] for r in results_fam):.4f}]")
    print(f"    Srednia: {np.mean([r['ratio3'] for r in results_fam]):.4f}")
    print(f"    e/sqrt(2) = {e_sqrt2:.4f}")
    all_close = all(abs(r['ratio3'] - e_sqrt2)/e_sqrt2 < 0.03 for r in results_fam)
    print(f"    Wszystkie w 3% od e/sqrt(2)? {'TAK' if all_close else 'NIE'}")
    print()
    print(f"  Pseudo-rodzina (alfa=8.56, lam zmienny):")
if pseudo_res:
    ratios_p = [x[4] for x in pseudo_res]
    print(f"    Zakres: [{min(ratios_p):.4f}, {max(ratios_p):.4f}]")
    print(f"    Stalosc (std/mean): {np.std(ratios_p)/np.mean(ratios_p)*100:.2f}%")
    print()

# ============================================================
# WYKRES
# ============================================================
print("Tworze wykres...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P21b: psi_core(M3)/psi_new -- siatka logarytmiczna', fontsize=11)

ax = axes[0]
if results_fam:
    alphas_f = [r['alpha'] for r in results_fam]
    ratios_f = [r['ratio3'] for r in results_fam]
    ax.plot(alphas_f, ratios_f, 'bo-', ms=8, lw=2, label='rodzina optymalna')
    ax.axhline(e_sqrt2, color='r', lw=1.5, ls='--', label=f'e/sqrt(2)={e_sqrt2:.4f}')
    ax.axhline(2.0, color='orange', lw=1, ls=':', label='2.0')
ax.set_xlabel('alpha'); ax.set_ylabel('psi_core/psi_new')
ax.set_title('Ratio vs alpha (rodzina optymalna)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1]
if pseudo_res:
    lams_pr = [x[0] for x in pseudo_res]
    ratios_pr = [x[4] for x in pseudo_res]
    ax.semilogx(lams_pr, ratios_pr, 'rs-', ms=6, lw=2, label='pseudo-rodzina')
    ax.axhline(e_sqrt2, color='k', lw=1.5, ls='--', label=f'e/sqrt(2)')
ax.set_xlabel('lam*'); ax.set_ylabel('psi_core/psi_new')
ax.set_title('Ratio vs lam* (alfa=8.56 fixed)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[2]
# Skalowanie K3 z lam (pseudo-rodzina)
if pseudo_res:
    lams_pr = np.array([x[0] for x in pseudo_res])
    K3s_pr  = np.array([x[1] for x in pseudo_res])
    ax.loglog(lams_pr, K3s_pr, 'bs-', ms=6, lw=2, label='K3(lam)')
    # Fit
    beta = c[0]
    K0 = np.exp(c[1])
    ax.loglog(lams_pr, K0*lams_pr**beta, 'r--', lw=1.5,
              label=f'fit: lam^{beta:.3f}')
ax.set_xlabel('lam*'); ax.set_ylabel('K3*')
ax.set_title('Skalowanie K3(lam*)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
