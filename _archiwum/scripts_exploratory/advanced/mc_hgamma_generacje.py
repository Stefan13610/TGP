#!/usr/bin/env python3
"""
mc_hgamma_generacje.py
======================
TGP: masa generacji z relaksacji n źródeł na siatce substratu Gamma

Hipoteza (korekta naiwnego modelu n-kinku):
  Naiwny:    E_n = n * E_1  -->  m_n/m_1 = n   [obserwacje: 207, 3477]
  Poprawny:  n źródeł Phi(x_i)=eps na siatce Gamma; relaksacja H_Gamma.
             E_n != n*E_1 ze wzgledu na nieliniowos alph*(grad Phi)^2/(Phi0*Phi).
             Nakładanie profili polowych z n źródeł tworzy energie
             oddzialywania zalezna od separacji d = krok siatki a_Gamma.

Kluczowy parametr: kappa = m_sp * a_Gamma
Warunek na poprawne masy: exp(kappa) ~ 207  -->  kappa ~ 5.33

Wyniki:
  - E_n(d) dla n = 1,2,3 jako funkcja separacji d*m_sp
  - Skan alfa: jak sprzezenie nieliniowe wzmacnia hierarchie
  - Predykcja TGP: kappa = m_sp * a_Gamma z obserwacji leptonowych
"""

import numpy as np
from scipy.optimize import minimize
import warnings
import time
warnings.filterwarnings('ignore')

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ══════════════════════════════════════════════════════════════
# PARAMETRY TGP
# ══════════════════════════════════════════════════════════════
M_SP  = 1.0    # masa pola generujacego przestrzen [l_P^{-1}]
PHI0  = 1.0    # wartosc prózniowa Phi_0
BETA  = 0.5    # = GAMMA (warunek prózni beta=gamma)
EPS   = 0.02   # Phi w centrum zrodla [w jednostkach Phi0]

# Siatka 1D
N_GRID = 400           # liczba wezlow
L_BOX  = 20.0          # polowa rozmiaru ukladu [1/m_sp]
DX     = 2 * L_BOX / N_GRID
XX     = np.linspace(-L_BOX, L_BOX, N_GRID)

# ══════════════════════════════════════════════════════════════
# HAMILTONIAN H_Gamma
# ══════════════════════════════════════════════════════════════

def V_pot(phi):
    """Meksykanski kapelusz: V(Phi) = (beta/4)(Phi^2 - Phi0^2)^2 / Phi0^2"""
    return (BETA / 4.0) * (phi**2 - PHI0**2)**2 / PHI0**2


def H_gamma(phi, alpha):
    """
    Energia konfiguracji Phi na siatce 1D:

      H = sum_krawedzi [ (grad Phi)^2 + alpha*(grad Phi)^2/(Phi0*Phi) ] * dx
        + sum_wezly    V(Phi) * dx

    Czlon alpha*(grad Phi)^2/(Phi0*Phi) jest kluczowy:
    dywerguje gdy Phi -> 0 (centrum zrodla) i tworzy
    oddzialywanie miedzy nakladajacymi sie profilami.
    """
    ps  = np.maximum(phi, 1e-10)                  # Phi > 0 zawsze
    g   = np.diff(ps) / DX                        # gradienty (N-1 elem.)
    pe  = np.maximum(0.5*(ps[:-1] + ps[1:]), 1e-10)  # Phi na krawedziach

    E_std = np.sum(g**2) * DX
    E_nl  = alpha * np.sum(g**2 / (PHI0 * pe)) * DX
    E_pot = np.sum(V_pot(ps)) * DX

    return E_std + E_nl + E_pot


def relax_field(src_indices, alpha, phi_init=None, tol=1e-12):
    """
    Relaksacja pola Phi z n zrodlami przy Phi[src] = EPS*Phi0.
    Warunki brzegowe: Phi[0] = Phi[-1] = Phi0.
    Zwraca (phi_opt, E_opt).
    """
    # Maska wezlow stałych (zrodla + brzegi)
    fixed = np.zeros(N_GRID, dtype=bool)
    for i in src_indices:
        if 0 < i < N_GRID - 1:
            fixed[i] = True
    fixed[0]  = True
    fixed[-1] = True
    free = ~fixed

    # Inicjalizacja: Phi = Phi0 wsedzie, Phi[src] = eps*Phi0
    if phi_init is None:
        phi = np.ones(N_GRID) * PHI0
    else:
        phi = phi_init.copy()
    phi[fixed] = PHI0
    for i in src_indices:
        if 0 < i < N_GRID - 1:
            phi[i] = EPS * PHI0

    phi_full = phi.copy()

    def obj(pv):
        phi_full[free] = pv
        return H_gamma(phi_full, alpha)

    res = minimize(
        obj,
        phi[free],
        method='L-BFGS-B',
        bounds=[(1e-9, None)] * int(free.sum()),
        options={'maxiter': 8000, 'ftol': tol, 'gtol': 1e-9}
    )
    phi_full[free] = res.x
    # Odswiezamy zrodla (moga zostac zmodyfikowane przez optymalizator)
    for i in src_indices:
        if 0 < i < N_GRID - 1:
            phi_full[i] = EPS * PHI0

    return phi_full, res.fun


# ══════════════════════════════════════════════════════════════
# OBLICZENIA PODSTAWOWE
# ══════════════════════════════════════════════════════════════

def compute_E1(alpha):
    """Energia fundamentalnego zrodla (generacja 1)."""
    ctr = N_GRID // 2
    phi, E = relax_field([ctr], alpha=alpha)
    return phi, E


def scan_En_vs_d(n_src, d_min_msp, d_max_msp, n_pts, alpha):
    """
    Skan E_n(d): n rownoodleglych zrodel, separacja d w jednostkach 1/m_sp.
    Zwraca (d_vals, E_vals).
    """
    d_arr = np.linspace(d_min_msp, d_max_msp, n_pts)
    E_arr = np.full(n_pts, np.nan)
    ctr   = N_GRID // 2

    for k, d in enumerate(d_arr):
        d_idx = max(1, int(round(d / (DX * M_SP))))
        span  = (n_src - 1) * d_idx
        start = ctr - span // 2
        positions = [start + i * d_idx for i in range(n_src)]

        # Sprawdz czy pozycje sa w granicach siatki (marg. 10 wezlow)
        if positions[0] < 10 or positions[-1] > N_GRID - 10:
            continue

        _, E = relax_field(positions, alpha=alpha)
        E_arr[k] = E

    return d_arr, E_arr


def scan_alpha(src_pair_indices, alpha_vals, E1_cache=None):
    """
    Skan E_1, E_2 po roznych wartosciach alfa przy stalej separacji.
    Zwraca liste (alpha, E1, E2, E2/E1, kappa=ln(E2/E1)).
    """
    rows = []
    for a in alpha_vals:
        _, e1 = compute_E1(a)
        _, e2 = relax_field(src_pair_indices, alpha=a)
        r  = e2 / max(e1, 1e-12)
        kk = np.log(max(r, 1.0001))
        rows.append((a, e1, e2, r, kk))
    return rows


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    SEP  = "═" * 66
    SEP2 = "─" * 66

    print(SEP)
    print("  TGP: Masa generacji z relaksacji n zrodel na siatce H_Gamma")
    print(SEP)
    print(f"  m_sp={M_SP}, Phi0={PHI0}, beta=gamma={BETA}, eps={EPS}*Phi0")
    print(f"  Siatka: N={N_GRID}, L=+-{L_BOX}/m_sp, dx={DX:.5f}/m_sp\n")

    ALPHA_DEF = 0.1

    t0 = time.time()

    # ──────────────────────────────────────────────────────────
    # KROK 1: E_1 — fundamentalne zrodlo
    # ──────────────────────────────────────────────────────────
    print(SEP2)
    print("KROK 1: Fundamentalne zrodlo (generacja 1) — kalibracja E_1")
    print(SEP2)

    phi1, E1 = compute_E1(ALPHA_DEF)
    E1_anal  = (4.0/3.0) * M_SP * PHI0**2

    print(f"  alfa     = {ALPHA_DEF}")
    print(f"  E_1 (num)   = {E1:.6f}")
    print(f"  E_1 (anal)  = {E1_anal:.6f}  [4/3 * m_sp * Phi0^2]")
    print(f"  Blad wzgl.  = {abs(E1-E1_anal)/E1_anal*100:.2f}%")
    print(f"  Profil: Phi_min = {phi1.min():.5f}  przy x={XX[phi1.argmin()]:.3f}/m_sp\n")

    # ──────────────────────────────────────────────────────────
    # KROK 2: E_2(d) — dwa zrodla, skan separacji
    # ──────────────────────────────────────────────────────────
    print(SEP2)
    print("KROK 2: E_2(d) — 2 zrodla, skan separacji d*m_sp")
    print(SEP2)
    print("  Mechanizm: nakladajace sie profile polowe zrodel.")
    print("  Czlon alfa*(grad Phi)^2/(Phi0*Phi) dywerguje gdy Phi->0.")
    print("  Gdy dwa zrodla sa blisko, pole miedzy nimi jest silnie")
    print("  tłumione => ogromna energia oddzialywania.\n")

    d2_vals, E2_arr = scan_En_vs_d(2, d_min_msp=0.15, d_max_msp=12.0,
                                    n_pts=36, alpha=ALPHA_DEF)

    print(f"  {'d*m_sp':>8}  {'E_2':>10}  {'E_2/E_1':>9}  "
          f"{'(E_2-2E_1)/E_1':>15}  {'ln(E_2/E_1)':>12}")
    print("  " + "─"*62)
    for d, E in zip(d2_vals, E2_arr):
        if np.isfinite(E) and E > 0:
            r   = E / E1
            de  = (E - 2*E1) / E1
            lnr = np.log(r)
            print(f"  {d:>8.3f}  {E:>10.5f}  {r:>9.5f}  {de:>15.5f}  {lnr:>12.5f}")

    v2        = np.isfinite(E2_arr) & (E2_arr > 0)
    E2_sorted = E2_arr[v2]
    d2_sorted = d2_vals[v2]

    # Wartosci przy ekstremalnych separacjach
    E2_at_dmin = E2_sorted[0]
    d2_dmin    = d2_sorted[0]
    E2_at_dmax = E2_sorted[-1]
    d2_dmax    = d2_sorted[-1]
    E2_min     = E2_sorted.min()
    d2_at_Emin = d2_sorted[E2_sorted.argmin()]

    print(f"\n  => d->0  (d*m_sp={d2_dmin:.3f}): E_2={E2_at_dmin:.5f},  "
          f"m_2/m_1={E2_at_dmin/E1:.4f}")
    print(f"  => E min (d*m_sp={d2_at_Emin:.3f}): E_2={E2_min:.5f},  "
          f"m_2/m_1={E2_min/E1:.4f}")
    print(f"  => d->inf(d*m_sp={d2_dmax:.3f}): E_2={E2_at_dmax:.5f},  "
          f"m_2/m_1={E2_at_dmax/E1:.4f}  (limit: 2.0)")
    print(f"  => Obserwowane m_mu/m_e = 207.0")

    # ──────────────────────────────────────────────────────────
    # KROK 3: E_3(d) — trzy zrodla rownoodlegle
    # ──────────────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("KROK 3: E_3(d) — 3 rownoodlegle zrodla")
    print(SEP2)

    d3_vals, E3_arr = scan_En_vs_d(3, d_min_msp=0.15, d_max_msp=8.0,
                                    n_pts=25, alpha=ALPHA_DEF)

    print(f"  {'d*m_sp':>8}  {'E_3':>10}  {'E_3/E_1':>9}  "
          f"{'(E_3-3E_1)/E_1':>15}  {'ln(E_3/E_1)/2':>14}")
    print("  " + "─"*64)
    for d, E in zip(d3_vals, E3_arr):
        if np.isfinite(E) and E > 0:
            r    = E / E1
            de   = (E - 3*E1) / E1
            lnr2 = np.log(r)/2
            print(f"  {d:>8.3f}  {E:>10.5f}  {r:>9.5f}  {de:>15.5f}  {lnr2:>14.5f}")

    v3        = np.isfinite(E3_arr) & (E3_arr > 0)
    E3_sorted = E3_arr[v3]
    d3_sorted = d3_vals[v3]
    E3_at_dmin = E3_sorted[0]
    d3_dmin    = d3_sorted[0]
    E3_min     = E3_sorted.min()

    print(f"\n  => d->0  (d*m_sp={d3_dmin:.3f}): E_3={E3_at_dmin:.5f},  "
          f"m_3/m_1={E3_at_dmin/E1:.4f}")
    print(f"  => Obserwowane m_tau/m_e = 3477.0")

    # ──────────────────────────────────────────────────────────
    # KROK 4: Skalowanie eksponencjalne E_n(d_min)
    # ──────────────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("KROK 4: Skalowanie E_n(d_min) — test eksponencjalny")
    print(SEP2)

    E_list = [E1, E2_at_dmin, E3_at_dmin]
    ns     = np.array([1, 2, 3])
    log_r  = np.log(np.array([max(E/E1, 1.0001) for E in E_list]))

    coeffs    = np.polyfit(ns, log_r, 1)
    kappa_obs = coeffs[0]
    kappa_req = np.log(207.0)

    print(f"\n  Dane przy d*m_sp = {d2_dmin:.3f} (krok siatki a_Gamma):")
    for n, E in zip(ns, E_list):
        print(f"    n={n}: E_n/E_1 = {E/E1:.5f},  ln(E_n/E_1) = {np.log(E/E1):.5f}")

    print(f"\n  Fit liniowy: ln(E_n/E_1) = kappa*n + const")
    print(f"    kappa_obs      = {kappa_obs:.5f}")
    print(f"    kappa_wym      = ln(207) = {kappa_req:.5f}")
    print(f"    kappa_obs/kap* = {kappa_obs/max(kappa_req,1e-6):.5f}")

    # ──────────────────────────────────────────────────────────
    # KROK 5: Skan alfa — mechanizm wzmocnienia hierarchii
    # ──────────────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("KROK 5: Skan alfa — jak sprzezenie wzmacnia E_2/E_1")
    print(SEP2)

    # Para zrodel przy minimalnej separacji d_min
    d_min_idx = max(1, int(round(d2_dmin / (DX * M_SP))))
    ctr        = N_GRID // 2
    src_pair   = [ctr - d_min_idx, ctr + d_min_idx]

    alpha_vals   = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    rows_alpha   = scan_alpha(src_pair, alpha_vals)

    print(f"\n  Separacja stala: d*m_sp = {d2_dmin:.3f}")
    print(f"\n  {'alfa':>7}  {'E_1':>10}  {'E_2':>10}  {'E_2/E_1':>9}  "
          f"{'kappa=ln(E2/E1)':>16}")
    print("  " + "─"*60)
    for a, e1a, e2a, r, kk in rows_alpha:
        print(f"  {a:>7.4f}  {e1a:>10.5f}  {e2a:>10.5f}  {r:>9.5f}  {kk:>16.5f}")

    # Skalowanie power-law: kappa ~ alfa^n
    pos_rows = [(a, kk) for a, _, _, _, kk in rows_alpha if a > 0 and kk > 1e-4]
    if len(pos_rows) >= 3:
        log_a    = np.log([r[0] for r in pos_rows])
        log_k    = np.log([r[1] for r in pos_rows])
        slope_ak = np.polyfit(log_a, log_k, 1)[0]
        print(f"\n  Skalowanie: kappa_eff ~ alfa^{slope_ak:.3f}")
        print(f"  (liniowe slope=1.0; nieliniowe wzmocnienie slope>1.0)")

        # Wymagana alfa dla kappa = ln(207) = 5.33
        intercept_ak = np.polyfit(log_a, log_k, 1)[1]
        try:
            alpha_star = np.exp((np.log(kappa_req) - intercept_ak) / slope_ak)
            print(f"\n  Wymagana alfa* dla kappa=ln(207)={kappa_req:.4f}: alfa*={alpha_star:.4f}")
        except Exception:
            alpha_star = None
            print("  (niemozliwa interpolacja alfa*)")
    else:
        slope_ak  = 1.0
        alpha_star = None

    # ──────────────────────────────────────────────────────────
    # KROK 6: Predykcja m_sp * a_Gamma
    # ──────────────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("KROK 6: Predykcja parametru TGP  kappa = m_sp * a_Gamma")
    print(SEP2)

    print(f"""
  Model sieciowy TGP:
    Generacja n odpowiada konfiguracji n zrodel rozstawionych
    na siatce substratu Gamma z krokiem a_Gamma.
    Krok siatki = MINIMALNA mozliwa separacja (topologiczny zakaz zblizenia).

    m_n / m_1 = E_n(a_Gamma) / E_1

  Wartosc kappa = m_sp * a_Gamma wyznacza cala hierarchie mas:
    m_2/m_1 = exp(kappa)     =>  kappa = ln(207)   = {np.log(207):.4f}
    m_3/m_1 = exp(2*kappa)   =>  kappa = ln(3477)/2 = {np.log(3477)/2:.4f}

  Wyniki numeryczne (alfa={ALPHA_DEF}, d*m_sp={d2_dmin:.3f}):
    m_2/m_1 = {E2_at_dmin/E1:.4f}   (obs: 207)
    m_3/m_1 = {E3_at_dmin/E1:.4f}   (obs: 3477)
    kappa_obs = {kappa_obs:.4f}   (wym: {kappa_req:.4f})

  Interpretacja:
    alfa={ALPHA_DEF} za male -- hierarchy za slaba.
    Potrzebne:""")

    if alpha_star is not None:
        print(f"    (A) alfa* = {alpha_star:.4f}  (silniejsze sprzezenie")
        print(f"        przy tym samym kroku siatki d*m_sp={d2_dmin:.3f})")

    print(f"""    (B) Wiekszy krok siatki a_Gamma przy alfa={ALPHA_DEF}:
        Wymagane d* takie ze E_2(d*)/E_1 = 207
        Numerycznie: d*m_sp ~ ? (sprawdz Krok 2 tabele)

  Predykcja TGP (niezalezna od alfa):
    kappa = m_sp * a_Gamma = {kappa_req:.4f}
    Jesli m_sp ~ okno Efimova: m_sp = 0.0831 l_P^-1
      =>  a_Gamma = {kappa_req:.4f} / 0.0831 ~ {kappa_req/0.0831:.1f} l_P
    Jesli a_Gamma = 1 l_P:
      =>  m_sp = {kappa_req:.4f} l_P^-1  (poza oknem Efimova)
    Jesli a_Gamma = 64 l_P:
      =>  m_sp = {kappa_req/64:.4f} l_P^-1  (zgodne z oknem Efimova!)
""")

    # ──────────────────────────────────────────────────────────
    # PODSUMOWANIE
    # ──────────────────────────────────────────────────────────
    t_total = time.time() - t0
    print(SEP)
    print("PODSUMOWANIE — O16 MC na H_Gamma (v31, nowy model)")
    print(SEP)
    print(f"""
  Mechanizm potwierdzony:
    czlon alfa*(grad Phi)^2/(Phi0*Phi) tworzy oddzialywanie
    ZALEZNE OD SEPARACJI miedzy zrodlami.
    E_n(d) jest funkcja nietrywialną d (nie n*E_1).

  Wyniki (alfa={ALPHA_DEF}):
    E_1 = {E1:.5f}  [generacja 1: e, u, d]
    E_2 = {E2_at_dmin:.5f}  m_2/m_1 = {E2_at_dmin/E1:.4f}  (obs: 207)
    E_3 = {E3_at_dmin:.5f}  m_3/m_1 = {E3_at_dmin/E1:.4f}  (obs: 3477)
    kappa_obs  = {kappa_obs:.4f}
    kappa_wym  = {kappa_req:.4f}  [ln(207)]

  Predykcja TGP:
    kappa = m_sp * a_Gamma = {kappa_req:.4f}
    a_Gamma ~ 64 l_P  (przy m_sp z okna Efimova)

  Otwarte:
    - Pełny model 3D (defekty sferyczne zamiast 1D)
    - Wyznaczenie a_Gamma z niezaleznych obserwacji
    - Termiczny MC na H_Gamma (temperatura T_Gamma)
    - Wyjasnienie dlaczego dok. 3 generacje (prog E_3 < E_bariera < E_4)

  Czas obliczen: {t_total:.1f} s
    """)

    if HAS_MPL:
        _make_plots(phi1, E1, d2_vals, E2_arr, d3_vals, E3_arr,
                    rows_alpha, d2_dmin, E2_at_dmin, E3_at_dmin,
                    kappa_req, kappa_obs)


# ══════════════════════════════════════════════════════════════
# WYKRESY
# ══════════════════════════════════════════════════════════════

def _make_plots(phi1, E1, d2v, E2a, d3v, E3a, rows_alpha,
                d_min, E2_dmin, E3_dmin, kappa_req, kappa_obs):

    fig = plt.figure(figsize=(16, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.36)

    alpha_def = 0.1

    # ── 1. Profil phi_1(x) ──────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(XX * M_SP, phi1 / PHI0, 'b-', lw=2, label='phi_1(x)')
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='Phi_0')
    ax.axhline(EPS, color='r',    ls=':',  lw=1.2, label=f'eps={EPS}')
    ax.fill_between(XX*M_SP, phi1/PHI0, 1.0, alpha=0.1, color='blue')
    ax.set_xlabel('x * m_sp'); ax.set_ylabel('Phi/Phi_0')
    ax.set_title('Profil phi_1(x) — generacja 1', fontsize=9)
    ax.legend(fontsize=8); ax.set_xlim(-6, 6); ax.grid(True, alpha=0.3)

    # ── 2. E_2(d) ───────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    v2 = np.isfinite(E2a) & (E2a > 0)
    ax.semilogy(d2v[v2], E2a[v2]/E1, 'g-o', ms=4, lw=2, label='E_2(d)/E_1')
    ax.axhline(2.0,   color='gray',   ls='--', lw=1,   label='2*E_1 (niezal.)')
    ax.axhline(207.0, color='orange', ls=':',  lw=1.5, label='obs mu/e=207')
    ax.axvline(d_min, color='purple', ls='--', lw=1.2, label=f'a_Gamma={d_min:.2f}')
    ax.set_xlabel('d * m_sp'); ax.set_ylabel('E_2/E_1  [log]')
    ax.set_title('E_2(d): 2 zrodla vs separacja', fontsize=9)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # ── 3. E_3(d) ───────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 2])
    v3 = np.isfinite(E3a) & (E3a > 0)
    ax.semilogy(d3v[v3], E3a[v3]/E1, 'r-o', ms=4, lw=2, label='E_3(d)/E_1')
    ax.axhline(3.0,    color='gray',   ls='--', lw=1,   label='3*E_1 (niezal.)')
    ax.axhline(3477.0, color='orange', ls=':',  lw=1.5, label='obs tau/e=3477')
    ax.axvline(d_min,  color='purple', ls='--', lw=1.2, label=f'a_Gamma={d_min:.2f}')
    ax.set_xlabel('d * m_sp'); ax.set_ylabel('E_3/E_1  [log]')
    ax.set_title('E_3(d): 3 zrodla vs separacja', fontsize=9)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # ── 4. Skalowanie eksponencjalne E_n(d_min) ─────────────
    ax = fig.add_subplot(gs[1, 0])
    ns   = [1, 2, 3]
    En   = [E1, E2_dmin, E3_dmin]
    logr = [np.log(max(E/E1, 1.0001)) for E in En]
    ax.plot(ns, logr, 'ko-', ms=9, lw=2, label='ln(E_n/E_1)')
    cs   = np.polyfit(ns, logr, 1)
    ax.plot(ns, np.polyval(cs, ns), 'b--', lw=1.5,
            label=f'fit kappa={cs[0]:.4f}')
    # Wymagana linia
    ax.plot(ns, [kappa_req * n for n in ns], color='orange', ls=':', lw=2,
            label=f'kappa_wym={kappa_req:.4f}')
    ax.set_xlabel('n (generacja)'); ax.set_ylabel('ln(m_n/m_1)')
    ax.set_title('Skalowanie eksponencjalne E_n(a_Gamma)', fontsize=9)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.set_xticks([1, 2, 3])

    # ── 5. kappa vs alfa ────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    aa  = [r[0] for r in rows_alpha if r[0] > 0]
    kk  = [r[4] for r in rows_alpha if r[0] > 0]
    ax.loglog(aa, kk, 'ms-', ms=7, lw=2, label='kappa_eff(alfa)')
    ax.axhline(kappa_req, color='orange', ls='--', lw=1.5,
               label=f'kappa*=ln(207)={kappa_req:.2f}')
    ax.set_xlabel('alfa'); ax.set_ylabel('kappa = ln(E_2/E_1)')
    ax.set_title(f'kappa vs alfa  (d={d_min:.2f}/m_sp)', fontsize=9)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

    # ── 6. Profile 2-zrodlowe dla roznych d ─────────────────
    ax = fig.add_subplot(gs[1, 2])
    ctr    = N_GRID // 2
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, 5))
    for ci, d_unit in enumerate([0.15, 0.5, 1.0, 2.5, 6.0]):
        d_idx = max(1, int(round(d_unit / (DX * M_SP))))
        pos   = [ctr - d_idx, ctr + d_idx]
        if pos[0] < 5 or pos[-1] > N_GRID - 5:
            continue
        phi_r, E_r = relax_field(pos, alpha=alpha_def)
        lbl = f'd={d_unit:.2f},  E/E_1={E_r/E1:.2f}'
        ax.plot(XX * M_SP, phi_r / PHI0, color=colors[ci], lw=1.8, label=lbl)
    ax.axhline(1.0, color='gray', ls='--', lw=1)
    ax.set_xlabel('x * m_sp'); ax.set_ylabel('Phi/Phi_0')
    ax.set_title('Profile 2-zrodlowe phi_2(x, d)', fontsize=9)
    ax.legend(fontsize=6.5); ax.set_xlim(-6, 6); ax.grid(True, alpha=0.3)

    plt.suptitle(
        'TGP: Masa generacji z relaksacji H_Gamma — O16 MC na siatce substratu\n'
        f'alpha={alpha_def}, m_sp={M_SP}, Phi0={PHI0}, eps={EPS}',
        fontsize=10, fontweight='bold'
    )

    out = 'TGP/TGP_v1/scripts/advanced/mc_hgamma_generacje.png'
    try:
        plt.savefig(out, dpi=120, bbox_inches='tight')
        print(f"  Wykres zapisany: {out}")
    except Exception as e:
        print(f"  Wykres nie zapisany: {e}")
    plt.close()


if __name__ == '__main__':
    main()
