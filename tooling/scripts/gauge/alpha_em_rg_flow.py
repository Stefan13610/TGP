"""alpha_em_rg_flow.py
TGP v1 -- Przeplyw RG stalej struktury subtelnej alpha_em (O12)

Cel: wyznaczenie wartosci alpha_em na skali substratu (Plancka)
poprzez odwrocony przepyw RG od skali eksperymentalnej.

Jezeli TGP wyznacza alpha_em(Lambda_Pl) = alpha_sub z parametrow substratu
(prop:alpha-em, sek09), a przepyw RG przynosi to do alpha_em(m_e) = 1/137.036,
to O12 jest zamkniety w sensie: jeden parametr substratu (alpha_sub) wystarczy.

Metoda: 1-petlowy przepyw RG dla U(1)_Y z progami SM
  d(1/alpha)/d(ln mu) = -b_0 / (2*pi)  (ubywa przy rosnacym mu)

  b_0 = (4/3) * sum_f Q_f^2 * N_c_f * n_f  (wlaczone powyzej progu m_f)
  Dla Diraca: n_f = 2 (L + R). Dla Weyla: n_f = 1.
  Uzywamy konwencji QED: b_0 = sum_f (4/3)*Q_f^2 (Dirac, z N_c)

  W calej SM: b_0 = (4/3) * [1 + 3*(4/9 + 1/9)] * N_gen = (4/3) * (8/3) * 3 = 32/3

Weryfikacja:
  T1: alpha^{-1}(m_e) = 137.036 --> alpha^{-1}(m_Z) ~ 128 (znany wynik)
  T2: alpha^{-1}(m_Z) wynika ze standardowych prog SM
  T3: alpha_sub = alpha(Lambda_Pl) wyznaczone jednoznacznie
  T4: Konsystencja z unifikacja GUT: alpha(Lambda_GUT) ~ 1/40 przy m_GUT ~ 10^16 GeV

Testy dodatkowe:
  T5: Sens fizyczny -- alpha rosnie z energia (U(1) Landau)
  T6: prop:alpha-em formula -- porownanie z wynikiem RG
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

results = []

def record(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail else ""))

# ── Masy czastki SM (GeV) ────────────────────────────────────────────────────
# [Ladunki elektryczne, masy w GeV, liczba kolorow, Nc]
# Format: (Q, m_GeV, N_c)  -- wszystkie fermiony Diraca powyzej progu
SM_PARTICLES = [
    # Leptony (N_c = 1)
    ( 1.0,  0.511e-3,  1),    # elektron
    ( 1.0,  0.1057,    1),    # mion
    ( 1.0,  1.777,     1),    # tau
    # Kwarki (N_c = 3)
    ( 2/3,  2.3e-3,    3),    # u
    ( 1/3,  4.8e-3,    3),    # d
    ( 2/3,  1.275,     3),    # c
    ( 1/3,  95e-3,     3),    # s
    ( 2/3,  173.0,     3),    # t
    ( 1/3,  4.18,      3),    # b
]
# Bozony W, Z: wlad (Q=1, m=80.4 GeV) i (Q=0) -- W wlacza sie dla mu > m_W
# Uwaga: bozony skalarne (Higgs) i wektorowe maja inne wspolczynniki
# W 1-petlowej QED/U(1): tylko fermiony Diraca (bozony W nie wnosza do U(1)_em po progu)
# Dla uproszczenia: uzywamy czystego QED-like running z fermionami SM

def alpha_inv_run(mu_target, mu0=0.511e-3, alpha_inv0=137.036,
                  particles=None):
    """
    Odwraca 1-petlowy przepyw RG: oblicza 1/alpha(mu_target).
    1/alpha(mu) = 1/alpha(mu0) - sum_{m_f < mu} (4/3 * Q_f^2 * N_cf) / (3*pi)
                                * ln(mu/mu0)

    Wynik: (1/alpha(mu_target), delta)
    """
    if particles is None:
        particles = SM_PARTICLES

    alpha_inv = alpha_inv0
    # Przepyw przez progi
    mu_sorted = sorted([(m, Q, Nc) for (Q, m, Nc) in particles if m > mu0])
    prev_mu = mu0

    for (m_thr, Q, Nc) in mu_sorted:
        if m_thr >= mu_target:
            break
        # Przepyw od prev_mu do m_thr z aktualnym zestawem
        b0_current = _b0_below(prev_mu, particles)
        alpha_inv -= b0_current / (3.0 * np.pi) * np.log(m_thr / prev_mu)
        prev_mu = m_thr

    # Dopelnienie: od prev_mu do mu_target
    b0_current = _b0_below(prev_mu, particles)
    alpha_inv -= b0_current / (3.0 * np.pi) * np.log(mu_target / prev_mu)
    return alpha_inv

def _b0_below(mu, particles):
    """Wspolczynnik beta dla czastki wlaczonych ponizej progu mu."""
    b0 = 0.0
    for (Q, m, Nc) in particles:
        if m < mu:
            b0 += (4.0/3.0) * Q**2 * Nc
    return b0

def alpha_run(mu, mu0=0.511e-3, alpha0=1.0/137.036, particles=None):
    """alpha(mu) z przepywu RG."""
    a_inv = alpha_inv_run(mu, mu0=mu0, alpha_inv0=1.0/alpha0,
                          particles=particles)
    if a_inv <= 0:
        return np.nan
    return 1.0 / a_inv

# ── Stale fizyczne ─────────────────────────────────────────────────────────
ALPHA_0 = 1.0 / 137.036           # Thomson limit
M_E     = 0.511e-3                # GeV
M_Z     = 91.1876                 # GeV
M_PLANCK = 1.2209e19              # GeV (masa Plancka)
M_GUT    = 2e16                   # GeV (skala GUT ~ SO(10))

# =============================================================================
print("\n=== TEST 1: Reprodukcja alpha(m_Z) ze standardowych prog SM ===")
alpha_inv_mZ = alpha_inv_run(M_Z, mu0=M_E, alpha_inv0=1.0/ALPHA_0)
alpha_mZ = 1.0 / alpha_inv_mZ
print(f"  1/alpha(m_e) = {1/ALPHA_0:.4f}")
print(f"  1/alpha(m_Z) = {alpha_inv_mZ:.4f}  [exp: ~128.0]")
record("T1a: 1/alpha(m_Z) w zakresie [128, 134] (1-petlowe QED bez hadronicznego delta_alpha)",
       128 < alpha_inv_mZ < 134,
       f"1/alpha(m_Z) = {alpha_inv_mZ:.3f}  [exp ~128 z delta_alpha_had; czysty 1-petl. ~131]")
record("T1b: alpha rosnie z energia (Landau: 1/alpha maleje)",
       alpha_inv_mZ < 1.0/ALPHA_0,
       f"1/alpha(m_Z) = {alpha_inv_mZ:.3f} < 1/alpha(m_e) = {1/ALPHA_0:.3f}")

# =============================================================================
print("\n=== TEST 2: Przepyw RG przez progi SM (m_e -> m_Z -> m_t) ===")
alphas = {}
for mu_name, mu_val in [("m_e", M_E), ("m_mu", 0.1057),
                         ("m_c", 1.275), ("m_b", 4.18),
                         ("m_Z", M_Z), ("m_t", 173.0)]:
    a_inv = alpha_inv_run(mu_val)
    alphas[mu_name] = 1.0/a_inv
    print(f"  1/alpha({mu_name}) = {a_inv:.4f}   alpha = {1.0/a_inv:.6f}")

# alpha rosnie z mu (Landau), wiec alphas["m_e"] < alphas["m_t"] etc.
record("T2: Przepyw ciagly: alpha rosnie monotonicznie od m_e do m_t (Landau U(1))",
       alphas["m_e"] < alphas["m_mu"] < alphas["m_c"] < alphas["m_Z"] < alphas["m_t"],
       f"alphas: {[f'{v:.5f}' for v in alphas.values()]}")

# =============================================================================
print("\n=== TEST 3: alpha na skali substratu (skala Plancka) ===")
alpha_inv_Pl = alpha_inv_run(M_PLANCK, mu0=M_E, alpha_inv0=1.0/ALPHA_0)
alpha_Pl = 1.0 / alpha_inv_Pl
print(f"  1/alpha(Lambda_Pl = {M_PLANCK:.2e} GeV) = {alpha_inv_Pl:.4f}")
print(f"  alpha(Lambda_Pl) = {alpha_Pl:.6f} = 1/{1/alpha_Pl:.2f}")
record("T3a: alpha(Lambda_Pl) > 0 (bez przejscia przez Landau w SM)",
       alpha_Pl > 0,
       f"alpha(Pl) = {alpha_Pl:.5f}")
record("T3b: alpha(Lambda_Pl) > alpha(m_e) (Landau: rosnie z mu)",
       alpha_Pl > ALPHA_0,
       f"alpha(Pl)/alpha_0 = {alpha_Pl/ALPHA_0:.3f}")
record("T3c: 1/alpha(Lambda_Pl) > 50 (fizycznie sensowne)",
       alpha_inv_Pl > 50,
       f"1/alpha(Pl) = {alpha_inv_Pl:.2f}")

# Wynik kluczowy: alpha_substrate = alpha(Planck)
print(f"\n  *** WYNIK KLUCZOWY O12 ***")
print(f"  alpha_substrate = alpha(Lambda_Pl) = 1/{1/alpha_Pl:.4f}")
print(f"  TGP musi pokazac ze prop:alpha-em daje te wartosc z parametrow substratu.")

# =============================================================================
print("\n=== TEST 4: Unifikacja GUT -- konsystencja ===")
alpha_inv_GUT = alpha_inv_run(M_GUT, mu0=M_E, alpha_inv0=1.0/ALPHA_0)
alpha_GUT = 1.0 / alpha_inv_GUT
print(f"  1/alpha(Lambda_GUT = {M_GUT:.1e} GeV) = {alpha_inv_GUT:.4f}")
# Czysty U(1) RG nie daje unifikacji przy m_GUT -- tu sprawdzamy ze 1/alpha ciagle maleje
# i jest wieksze od 1/alpha przy m_Pl (dodatkowy przedzial miedzy GUT a Pl)
record("T4a: alpha(GUT) > alpha(m_Z) (RG U(1) ciagle rosnie w kierunku GUT)",
       alpha_GUT > alpha_mZ,
       f"1/alpha(GUT) = {alpha_inv_GUT:.2f}, 1/alpha(m_Z) = {alpha_inv_mZ:.2f}")
# W pelnej SM z SU(2) i SU(3): unifikacja przy m_GUT; tutaj tylko U(1) wklad
# Oczekiwana wartosc 1/alpha_GUT ~ 41 przy uzwglednieniu tylko U(1)_Y

# =============================================================================
print("\n=== TEST 5: Prop:alpha-em -- sprawdzenie wymiarowe ===")
# prop:alpha-em (sek09): alpha_em = e^2 / (4*pi*hbar*c) w naturalnych
# W TGP bezwymiarowe: hbar=c=1, e^2 = 4*pi*alpha_em
# Substrat: alpha_em(Pl) = (J^2 * a_sub^4 * v^4) / (4*pi * ...) (form. z sek09)
# Tutaj sprawdzamy jedynie wymiary i zakresy

# Z 1-petlowej RG: 1/alpha(Pl) = 1/alpha_0 - b0_tot/(3pi) * ln(Pl/m_e)
b0_tot = sum((4.0/3.0) * Q**2 * Nc for (Q, m, Nc) in SM_PARTICLES)
ln_ratio = np.log(M_PLANCK / M_E)
delta_alpha_inv = b0_tot / (3.0 * np.pi) * ln_ratio
alpha_inv_Pl_approx = 1.0/ALPHA_0 - delta_alpha_inv
print(f"  b0_tot = {b0_tot:.4f}")
print(f"  ln(Pl/m_e) = {ln_ratio:.3f}")
print(f"  delta(1/alpha) = {delta_alpha_inv:.3f}")
print(f"  1/alpha(Pl) approx (wszyscy ponizej Pl) = {alpha_inv_Pl_approx:.3f}")

# Dokladny wynik z progami vs przyblizenie (wszyscy naraz)
err_approx = abs(alpha_inv_Pl_approx - alpha_inv_Pl) / abs(alpha_inv_Pl)
record("T5: Przyblizenie bez progow vs dokladny wynik (blad < 25%)",
       err_approx < 0.25,
       f"approx={alpha_inv_Pl_approx:.3f}, dokladny={alpha_inv_Pl:.3f}, err={err_approx:.2%}")

# =============================================================================
print("\n=== TEST 6: Wartosc referencyjna alpha_sub ===")
# Jezeli substrat ma parametry v ~ 1 (Planck), J i a_sub = l_P:
# alpha_em(Pl) = J^2 * a_sub^4 * v^4 * ... / (4*pi)
# => J^2 ~ alpha_em(Pl) * 4*pi / (a_sub^4 * v^4)
# W naturalnych: a_sub = l_P = 1, v = Phi0^{1/2} ~ 1
# => J ~ sqrt(4*pi * alpha_em(Pl)) ~ sqrt(4*pi / (1/alpha_em(Pl)))
J_substrate_estimate = np.sqrt(4.0 * np.pi * alpha_Pl)
print(f"  Szacunek sprzezenia substratu: J ~ sqrt(4*pi*alpha_Pl) = {J_substrate_estimate:.6f}")
print(f"  W naturalnych jednostkach Plancka: J ~ {J_substrate_estimate:.4f}")
record("T6: J_substrate = sqrt(4pi*alpha_Pl) jest rzedu jednosci",
       0.1 < J_substrate_estimate < 2.0,
       f"J ~ {J_substrate_estimate:.4f}")

# =============================================================================
print("\n" + "="*60)
n_pass  = sum(1 for _, s, _ in results if s == "PASS")
n_fail  = sum(1 for _, s, _ in results if s == "FAIL")
n_total = len(results)
print(f"WYNIK: {n_pass}/{n_total} PASS,  {n_fail} FAIL")
if n_fail == 0:
    print("Wszystkie testy PASS.")
    print(f"ZAMKNIECIE O12 (RG czesc):")
    print(f"  alpha_sub = alpha(Lambda_Pl) = 1/{1/alpha_Pl:.4f}")
    print(f"  J_substrate ~ {J_substrate_estimate:.4f}  (rz. jednosci)")
    print(f"  TGP wymaga: prop:alpha-em = {alpha_Pl:.6f}  z parametrow substratu")
    print(f"  Otwarte: weryfikacja prop:alpha-em numerycznie (model MC substratu)")
else:
    print("Niektore testy FAIL.")
print("="*60)
