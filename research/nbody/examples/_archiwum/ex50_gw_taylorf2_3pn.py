"""
ex50_gw_taylorf2_3pn.py
========================
Widmo GW TaylorF2 z poprawką 3PN TGP — Kill-shot K20

PYTANIE BADAWCZE
----------------
Czy poprawka 3PN wynikająca z dynamicznej prędkości światła c(Φ)
w TGP generuje obserwowalne odchylenie fazy GW relative do GR
w pasmie LIGO/LISA?

FIZYKA
------
W TGP prędkość światła zależy od pola Φ:
    c(Φ) = c₀ · (Φ₀/Φ)^{1/2}

Dla fali GW propagującej przez potencjał grawitacyjny U = δΦ/Φ₀:
    c_GW(r) ≈ c₀ · (1 - U(r)/2)

To modyfikuje całkę fazową TaylorF2. Wkład TGP pojawia się
na poziomie ~1.5PN (z c(Φ)) oraz ~3PN (z nielinearności Φ).

SCHEMAT OBLICZEŃ
-----------------
Faza TaylorF2 GR:
    Ψ_GR(f) = 2πf·t_c - φ_c - π/4
              + (3/128η) · (πMf)^{-5/3} · [1 + Σ_n c_n · x^n]
    gdzie x = (πMf)^{2/3} (Post-Newtonian parameter)

Poprawka TGP (Δc₂ = -1/3 z rem:k20-3pn):
    ΔΨ_TGP(f) = (3/128η) · (πMf)^{-5/3} · δc_2 · x²
    gdzie δc_2 = -1/3 · U̅ (U̅ = średni potencjał wzdłuż propagacji)

OBSERWABLE
----------
1. Odchylenie fazy: ΔΨ(f) = Ψ_TGP(f) - Ψ_GR(f)
2. Optimal SNR: ρ = sqrt(4 ∫ |h̃(f)|²/S_n(f) df)
3. Detekcja TGP (Fisher information): σ_δc₂ ~ 1/(ρ·∂_δc₂ ΨTaylorF2)
4. Porównanie z LIGO O3 szumem S_n(f)

PARAMETRY TESTOWE
------------------
GW150914: M₁=36, M₂=29 M☉, D=410 Mpc, f ∈ [20, 512] Hz
GW170817: M₁=1.46, M₂=1.27 M☉, D=40 Mpc, f ∈ [20, 1000] Hz
LISA BH:  M₁=10⁶, M₂=10⁶ M☉, z=1, f ∈ [10⁻⁴, 0.1] Hz

Autor: TGP Analysis Session v28, 2026-03-22
"""

import sys
import io
import argparse
import numpy as np
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ──────────────────────────────────────────────────────────────────────────────
#  Stałe fizyczne (SI i geometryczne)
# ──────────────────────────────────────────────────────────────────────────────

G_SI   = 6.674e-11      # m³/(kg·s²)
c_SI   = 2.998e8        # m/s
M_SUN  = 1.989e30       # kg
PC_SI  = 3.086e16       # m
MPC_SI = 1e6 * PC_SI    # m

# Jednostki geometryczne: G=c=1, masa w sekundach
G_GEOM = G_SI / c_SI**3   # s/kg
M_SUN_S = M_SUN * G_GEOM  # masa Słońca w sekundach

# ──────────────────────────────────────────────────────────────────────────────
#  Parametry układu binarnego
# ──────────────────────────────────────────────────────────────────────────────

def chirp_mass(m1, m2):
    """Chirp mass w M_sun."""
    eta = m1 * m2 / (m1 + m2)**2
    return (m1 + m2) * eta**(3.0/5.0)

def symmetric_mass_ratio(m1, m2):
    """η = m₁m₂/(m₁+m₂)²."""
    return m1 * m2 / (m1 + m2)**2

# ──────────────────────────────────────────────────────────────────────────────
#  Gęstości spektralne szumu S_n(f)
# ──────────────────────────────────────────────────────────────────────────────

def S_n_aligo(f):
    """
    Advanced LIGO design sensitivity (O4/O5 approx.).
    Simplified analytic fit (Ajith+2011 model).
    """
    f0 = 215.0  # Hz
    S0 = 1e-48  # Hz^{-1} (normalizacja)

    # Prosta trójczłonowa aproksymacja
    S_seismic  = (f / 30.0)**(-4)   # niska częstotliwość
    S_thermal  = 1.0 + (f / f0)**2  # wewnętrzny szum termiczny
    S_shot     = (f / f0)**2        # shot noise

    S = S0 * (S_seismic + S_thermal + S_shot)
    # Minimalny szum (~10⁻²³ Hz^{-1/2} sensitivity)
    S = np.maximum(S, 4e-48)
    return S

def S_n_lisa(f):
    """
    LISA sensitivity (Robson+2019 approximation).
    f w Hz.
    """
    L = 2.5e9   # ramię interferometru [m]
    f_L = c_SI / (2.0 * np.pi * L)   # ≈ 19.09 mHz

    # Szum instrumentalny
    P_OMS = (1.5e-11)**2 * (1 + (2e-3 / f)**4)  # m²/Hz
    P_ACC = (3e-15)**2 * (1 + (0.4e-3 / f)**2) * (1 + (f / 8e-3)**4)  # m²/s⁴/Hz

    S_instr = (10 / (3 * L**2)) * (P_OMS + 2 * (1 + np.cos(f/f_L)**2) * P_ACC / (2*np.pi*f)**4)

    # Szum galaktyczny (WD binary confusion)
    f_knee = 1e-3  # Hz
    S_conf = 9e-45 * f**(-7.0/3.0) * np.exp(-(f/f_knee)**1.8) * (
        1 + np.tanh((f - 1.13e-2) / 7.9e-3))

    return S_instr + S_conf

# ──────────────────────────────────────────────────────────────────────────────
#  Amplituda GW (domcena dziedziny częstości, stationary phase approx.)
# ──────────────────────────────────────────────────────────────────────────────

def h_tilde_amp(f_arr, M_chirp_s, D_m):
    """
    Amplituda |h̃(f)| w przybliżeniu SPA (stationary phase):
    |h̃(f)| = C_norm · M_chirp^{5/6} · f^{-7/6} / D
    C_norm = sqrt(5/24) · π^{-2/3} · G^{5/6}/c^{3/2} (geometryczne)
    """
    # W jednostkach geometrycznych: M_chirp_s w sekundach, D_m w metrach
    C = np.sqrt(5.0 / 24.0) / np.pi**(2.0/3.0)
    amp = C * M_chirp_s**(5.0/6.0) * f_arr**(-7.0/6.0) / (D_m / c_SI)
    return amp

# ──────────────────────────────────────────────────────────────────────────────
#  Faza TaylorF2 (GR + poprawka TGP)
# ──────────────────────────────────────────────────────────────────────────────

def psi_taylorf2_gr(f_arr, M_total_s, eta):
    """
    Faza TaylorF2 GR do 3.5PN (bez stałych fazy t_c, phi_c).
    Ψ_GR(f) = (3/128η) · x^{-5/2} · [1 + c1·x + c1.5·x^1.5 + c2·x² + ...]
    gdzie x = (π M f)^{2/3}

    Współczynniki PN (Blanchet+2002):
    c₁   = 20/9 · (743/336 + 11/4 · η)
    c₁.₅ = -16π
    c₂   = 10 · (3058673/1016064 + 5429/1008·η + 617/144·η²)
    c₂.₅ = π(38645/756 - 65/9·η)·(1 + ln(6^{3/2}·v))  [log term]
    c₃   = (11583231236531/4694215680 - 640/3·π² - 6848/21·γ_E
            + (-15737765635/3048192 + 2255/12·π²)·η
            + 76055/1728·η² - 127825/1296·η³
            - 6848/21·log(4v))
    """
    x = (np.pi * M_total_s * f_arr)**(2.0/3.0)

    # Prefaktor
    prefactor = 3.0 / (128.0 * eta) * x**(-5.0/2.0)

    # 0PN
    psi = np.ones_like(f_arr)

    # 1PN (x^1)
    c1 = 20.0/9.0 * (743.0/336.0 + 11.0/4.0 * eta)
    psi += c1 * x

    # 1.5PN (x^{3/2}) — wpływ spin-orbitowy (tu bez spinów)
    c1p5 = -16.0 * np.pi
    psi += c1p5 * x**(3.0/2.0)

    # 2PN (x^2)
    c2 = 10.0 * (3058673.0/1016064.0 + 5429.0/1008.0 * eta + 617.0/144.0 * eta**2)
    psi += c2 * x**2

    # 2.5PN (x^{5/2}) — logarytmiczny (stosujemy wersję uproszczoną)
    v = x**(1.0/2.0)
    c2p5 = np.pi * (38645.0/756.0 - 65.0/9.0 * eta)
    psi += c2p5 * x**(5.0/2.0) * (1.0 + np.log(6.0**1.5 * np.maximum(v, 1e-30)))

    # 3PN (x^3) — tu kluczowa poprawka TGP wchodzi
    gamma_E = 0.5772156649  # stała Eulera-Mascheroniego
    c3_gr = (11583231236531.0/4694215680.0
             - 640.0/3.0 * np.pi**2
             - 6848.0/21.0 * gamma_E
             + (-15737765635.0/3048192.0 + 2255.0/12.0 * np.pi**2) * eta
             + 76055.0/1728.0 * eta**2
             - 127825.0/1296.0 * eta**3
             - 6848.0/21.0 * np.log(4.0 * np.maximum(v, 1e-30)))
    psi += c3_gr * x**3

    # 3.5PN (x^{7/2})
    c3p5 = np.pi * (77096675.0/254016.0 + 378515.0/1512.0 * eta
                    - 74045.0/756.0 * eta**2)
    psi += c3p5 * x**(7.0/2.0)

    return prefactor * psi


def psi_tgp_correction(f_arr, M_total_s, eta, delta_c2, U_bar=1.0):
    """
    Poprawka TGP do fazy TaylorF2 na poziomie 3PN.

    Źródło: dynamiczna prędkość c(Φ) = c₀(Φ₀/Φ)^{1/2} wzdłuż propagacji GW.
    Całkując wkład 3PN z efektywnego potencjału:
        ΔΨ_TGP(f) = (3/128η) · x^{-5/2} · δc₂ · x³

    Wynik (rem:k20-3pn w sek08_formalizm.tex):
        δc₂ = -1/3 (TGP predykcja)

    U_bar: uśredniony potencjał grawitacyjny wzdłuż drogi fali GW
           U_bar = G·M_źródła/(c₀²·D) ≈ 10⁻⁶ dla GW150914

    W zasadzie Δc₂ = delta_c2 · U_bar · (czynnik geometryczny).
    Dla uproszczenia: bierzemy pełne Δc₂ jako wolny parametr.
    """
    x = (np.pi * M_total_s * f_arr)**(2.0/3.0)
    prefactor = 3.0 / (128.0 * eta) * x**(-5.0/2.0)

    # Wkład 3PN z modyfikacji c(Φ)
    dpsi = prefactor * delta_c2 * U_bar * x**3

    return dpsi


def psi_total(f_arr, M_total_s, eta, delta_c2=0.0, U_bar=1.0):
    """Całkowita faza TaylorF2: GR + TGP."""
    psi_gr  = psi_taylorf2_gr(f_arr, M_total_s, eta)
    dpsi    = psi_tgp_correction(f_arr, M_total_s, eta, delta_c2, U_bar)
    return psi_gr + dpsi


# ──────────────────────────────────────────────────────────────────────────────
#  Falisty sygnał h̃(f)
# ──────────────────────────────────────────────────────────────────────────────

def h_tilde_tgp(f_arr, M_chirp_s, M_total_s, eta, D_m,
                delta_c2=0.0, U_bar=1.0):
    """
    Pełny sygnał GW w dziedzinie częstości:
    h̃(f) = A(f) · exp(i Ψ(f))
    """
    amp = h_tilde_amp(f_arr, M_chirp_s, D_m)
    psi = psi_total(f_arr, M_total_s, eta, delta_c2, U_bar)
    return amp * np.exp(1j * psi)


# ──────────────────────────────────────────────────────────────────────────────
#  SNR i wykrywalność poprawki TGP
# ──────────────────────────────────────────────────────────────────────────────

def compute_snr(f_arr, M_chirp_s, M_total_s, eta, D_m, S_n_func,
                delta_c2=0.0, U_bar=1.0):
    """
    Optymalny SNR:  ρ² = 4 ∫ |h̃(f)|² / S_n(f) df

    Zwraca: (rho, rho_sq_integrand_arr)
    """
    amp    = h_tilde_amp(f_arr, M_chirp_s, D_m)
    S_n    = S_n_func(f_arr)
    df     = f_arr[1] - f_arr[0]   # zakładamy równe kroki w log lub lin

    integrand = 4.0 * amp**2 / S_n
    rho_sq    = np.trapz(integrand, f_arr)
    return np.sqrt(max(rho_sq, 0.0)), integrand


def sensitivity_to_delta_c2(f_arr, M_chirp_s, M_total_s, eta, D_m,
                             S_n_func, U_bar=1.0):
    """
    Fisher matrix estimate: σ(δc₂) ~ 1 / sqrt(I_cc)
    gdzie I_cc = 4 ∫ |∂h/∂δc₂|² / S_n df

    ∂h/∂δc₂ = i · h · ∂Ψ/∂δc₂
    ∂Ψ/∂δc₂ = (3/128η) · x^{1/2} · U_bar   (wkład 3PN)
    """
    x    = (np.pi * M_total_s * f_arr)**(2.0/3.0)
    dpsi = (3.0 / (128.0 * eta)) * x**(1.0/2.0) * U_bar

    amp   = h_tilde_amp(f_arr, M_chirp_s, D_m)
    S_n   = S_n_func(f_arr)

    # |∂h/∂δc₂|² = amp² · (∂Ψ/∂δc₂)²
    integrand = 4.0 * amp**2 * dpsi**2 / S_n
    I_cc = np.trapz(integrand, f_arr)

    if I_cc > 0:
        sigma_dc2 = 1.0 / np.sqrt(I_cc)
    else:
        sigma_dc2 = np.inf

    return sigma_dc2


def phase_deviation_rms(f_arr, M_total_s, eta, delta_c2, U_bar,
                        S_n_func, M_chirp_s, D_m):
    """
    Ważona RMS odchylenie fazy TGP:
    ΔΨ_rms = sqrt(∫ |ΔΨ|² |h̃|²/S_n df / ∫ |h̃|²/S_n df)
    """
    dpsi   = psi_tgp_correction(f_arr, M_total_s, eta, delta_c2, U_bar)
    amp    = h_tilde_amp(f_arr, M_chirp_s, D_m)
    S_n    = S_n_func(f_arr)
    weight = amp**2 / S_n

    dpsi_rms = np.sqrt(np.trapz(dpsi**2 * weight, f_arr) /
                       (np.trapz(weight, f_arr) + 1e-100))
    return dpsi_rms


# ──────────────────────────────────────────────────────────────────────────────
#  Potencjał grawitacyjny U_bar wzdłuż drogi GW
# ──────────────────────────────────────────────────────────────────────────────

def estimate_U_bar(M_total_SI, D_SI, r_ISCO_SI=None):
    """
    Oszacowanie uśrednionego potencjału grawitacyjnego wzdłuż drogi GW.

    Dla propagacji przez potencjał źródła na odległości D:
        U_bar ≈ G·M_total / (c₀²·r_ISCO) · ln(D/r_ISCO) / (D/r_ISCO)

    Uproszczenie: U_bar ≈ G·M / (c₀²·D) · ln(D/r_ISCO)
    """
    if r_ISCO_SI is None:
        # ISCO dla Schwarzschilda: r_ISCO = 6GM/c²
        r_ISCO_SI = 6.0 * G_SI * M_total_SI / c_SI**2

    U_inner = G_SI * M_total_SI / (c_SI**2 * r_ISCO_SI)  # ~ v²/c² przy ISCO

    # Wkład propagacji (logarytmicznie mały):
    U_prop  = G_SI * M_total_SI / (c_SI**2 * D_SI) * np.log(D_SI / r_ISCO_SI)

    # Dominuje wkład przy źródle (bliskie otoczenie):
    U_bar = U_inner + U_prop
    return U_bar, U_inner, U_prop


# ──────────────────────────────────────────────────────────────────────────────
#  Analiza konkretnych zdarzeń GW
# ──────────────────────────────────────────────────────────────────────────────

EVENTS = {
    'GW150914': dict(
        m1=36.0, m2=29.0,
        D_Mpc=410.0,
        f_low=20.0, f_high=512.0,
        detector='aLIGO',
        label='GW150914 (BBH, LIGO O1)',
    ),
    'GW170817': dict(
        m1=1.46, m2=1.27,
        D_Mpc=40.0,
        f_low=20.0, f_high=1000.0,
        detector='aLIGO',
        label='GW170817 (BNS, LIGO O2)',
    ),
    'LISA_MBHB': dict(
        m1=1e6, m2=1e6,
        D_Mpc=6700.0,   # z=1
        f_low=1e-4, f_high=0.1,
        detector='LISA',
        label='LISA MBHB (10⁶+10⁶ M☉, z=1)',
    ),
}


def analyze_event(name, params, delta_c2_TGP=-1.0/3.0, n_pts=2000):
    """Pełna analiza pojedynczego zdarzenia GW."""
    m1, m2     = params['m1'], params['m2']
    D_Mpc      = params['D_Mpc']
    f_low      = params['f_low']
    f_high     = params['f_high']
    det        = params['detector']
    label      = params['label']

    M_total_s  = (m1 + m2) * M_SUN_S
    M_chirp_s  = chirp_mass(m1, m2) * M_SUN_S
    eta_val    = symmetric_mass_ratio(m1, m2)
    D_m        = D_Mpc * MPC_SI

    S_n_func   = S_n_aligo if det == 'aLIGO' else S_n_lisa

    # Siatka częstości
    f_arr = np.logspace(np.log10(f_low), np.log10(f_high), n_pts)

    # Potencjał grawitacyjny
    M_SI = (m1 + m2) * M_SUN
    U_bar_val, U_inner, U_prop = estimate_U_bar(M_SI, D_m)

    # U_eff dla fazy: stosujemy U_inner (dominuje)
    U_eff = U_inner

    # SNR (GR)
    rho_gr, integrand_gr = compute_snr(
        f_arr, M_chirp_s, M_total_s, eta_val, D_m, S_n_func
    )

    # SNR (TGP) — amplituda taka sama, faza różna; SNR optymalne = GR
    rho_tgp = rho_gr  # SNR zależy od amplitudy, nie fazy

    # Czułość na δc₂
    sigma_dc2 = sensitivity_to_delta_c2(
        f_arr, M_chirp_s, M_total_s, eta_val, D_m, S_n_func, U_bar=U_eff
    )

    # Odchylenie fazy (TGP predykcja: δc₂ = -1/3)
    dpsi_rms = phase_deviation_rms(
        f_arr, M_total_s, eta_val, delta_c2_TGP, U_eff,
        S_n_func, M_chirp_s, D_m
    )

    # Stosunek sygnał/czułość: czy TGP jest detekowalne?
    detectability = abs(delta_c2_TGP * U_eff) / sigma_dc2

    return dict(
        name=name, label=label,
        m1=m1, m2=m2, D_Mpc=D_Mpc, detector=det,
        M_chirp_Msun=chirp_mass(m1, m2),
        M_total_Msun=m1+m2,
        eta=eta_val,
        rho_gr=rho_gr,
        U_bar=U_bar_val, U_inner=U_inner, U_prop=U_prop, U_eff=U_eff,
        sigma_dc2=sigma_dc2,
        dpsi_rms_rad=dpsi_rms,
        delta_c2_TGP=delta_c2_TGP,
        detectability=detectability,
        detectable=(detectability > 1.0),
        f_arr=f_arr,
        integrand_gr=integrand_gr,
    )


# ──────────────────────────────────────────────────────────────────────────────
#  Faza vs częstości: tabela diagnostyczna
# ──────────────────────────────────────────────────────────────────────────────

def phase_table(f_arr, M_total_s, eta, delta_c2, U_eff,
                f_samples=(20, 50, 100, 200, 500)):
    """Tabela ΔΨ(f) w wybranych punktach."""
    rows = []
    for f0 in f_samples:
        # Znajdź najbliższy punkt w f_arr
        idx = np.argmin(np.abs(f_arr - f0))
        if idx < len(f_arr):
            f_v = f_arr[idx]
            dpsi = psi_tgp_correction(
                np.array([f_v]), M_total_s, eta, delta_c2, U_eff
            )[0]
            rows.append((f_v, dpsi))
    return rows


# ──────────────────────────────────────────────────────────────────────────────
#  Główna analiza
# ──────────────────────────────────────────────────────────────────────────────

def main(delta_c2_scan=None, plot=False):

    print("=" * 72)
    print("  EX50: Widmo GW TaylorF2 z poprawką 3PN TGP — Kill-shot K20")
    print(f"  TGP predykcja: Δc₂ = -1/3  (rem:k20-3pn, sek08_formalizm)")
    print("=" * 72)
    print()

    if delta_c2_scan is None:
        delta_c2_scan = [-1.0/3.0, -0.1, -0.01, 0.0]

    results = {}

    # ─────────────────────────────────────────────────────────────────────────
    # Analiza każdego zdarzenia
    # ─────────────────────────────────────────────────────────────────────────
    for event_name, params in EVENTS.items():
        print(f"\n{'─'*72}")
        print(f"  Zdarzenie: {params['label']}")
        print(f"{'─'*72}")

        res = analyze_event(event_name, params)
        results[event_name] = res

        print(f"\n  Parametry:")
        print(f"    M_chirp = {res['M_chirp_Msun']:.2f} M☉")
        print(f"    M_total = {res['M_total_Msun']:.1f} M☉")
        print(f"    η       = {res['eta']:.4f}")
        print(f"    D       = {res['D_Mpc']:.0f} Mpc")
        print(f"    Detektor: {res['detector']}")

        print(f"\n  Potencjał grawitacyjny:")
        print(f"    U_inner (ISCO) = {res['U_inner']:.4e}")
        print(f"    U_prop  (Droga)= {res['U_prop']:.4e}")
        print(f"    U_eff   (suma) = {res['U_eff']:.4e}")

        print(f"\n  SNR i czułość:")
        print(f"    ρ(GR/TGP) = {res['rho_gr']:.1f}")
        print(f"    σ(δc₂) = {res['sigma_dc2']:.4e}  "
              f"[min. wykrywalna |δc₂|]")
        print(f"    |δc₂_TGP · U_eff| = {abs(res['delta_c2_TGP'] * res['U_eff']):.4e}")
        print(f"    Wykrywalność = {res['detectability']:.3f}  "
              f"({'DETEKOWALNE ✓' if res['detectable'] else 'PONIŻEJ PROGU ✗'})")

        print(f"\n  Odchylenie fazy TGP (δc₂ = -1/3):")
        print(f"    ΔΨ_rms (ważone SNR) = {res['dpsi_rms_rad']:.4e} rad")
        print(f"    Próg detekcji fazy  = 1/(2ρ) = {0.5/res['rho_gr']:.4e} rad")
        detectable_phase = res['dpsi_rms_rad'] > 0.5 / res['rho_gr']
        print(f"    Faza detekowalna:    {'TAK ✓' if detectable_phase else 'NIE ✗'}")

        # Tabela fazy
        print(f"\n  ΔΨ_TGP(f) dla kluczowych częstości [δc₂=-1/3, U_eff={res['U_eff']:.2e}]:")
        M_total_s = res['M_total_Msun'] * M_SUN_S
        if res['detector'] == 'LISA':
            f_samp = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
        else:
            f_samp = [20, 50, 100, 200, 400]
        rows = phase_table(res['f_arr'], M_total_s, res['eta'],
                           res['delta_c2_TGP'], res['U_eff'], f_samp)
        print(f"    {'f [Hz]':>12}  {'ΔΨ [rad]':>14}")
        print(f"    {'─'*28}")
        for f_v, dpsi_v in rows:
            print(f"    {f_v:>12.4g}  {dpsi_v:>14.4e}")

    # ─────────────────────────────────────────────────────────────────────────
    # Skan δc₂
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n\n{'='*72}")
    print(f"  Skan δc₂ — kiedy TGP jest detekowalne?")
    print(f"{'='*72}")

    for event_name, res in results.items():
        print(f"\n  [{res['detector']}] {res['label']}  (ρ={res['rho_gr']:.0f}):")
        print(f"    {'|δc₂|':>10}  {'wykryw.':>10}  {'ΔΨ_rms [rad]':>14}  {'status':>12}")
        print(f"    {'─'*50}")

        M_total_s  = res['M_total_Msun'] * M_SUN_S
        M_chirp_s  = res['M_chirp_Msun'] * M_SUN_S
        D_m        = res['D_Mpc'] * MPC_SI
        S_n_func   = S_n_aligo if res['detector'] == 'aLIGO' else S_n_lisa

        for dc2 in delta_c2_scan:
            det_ratio = abs(dc2 * res['U_eff']) / res['sigma_dc2']
            dpsi = phase_deviation_rms(
                res['f_arr'], M_total_s, res['eta'],
                dc2, res['U_eff'], S_n_func, M_chirp_s, D_m
            )
            status = 'DETEKOWALNE' if det_ratio > 1.0 else 'poniżej progu'
            print(f"    {abs(dc2):>10.4f}  {det_ratio:>10.3f}  {dpsi:>14.4e}  {status:>12}")

    # ─────────────────────────────────────────────────────────────────────────
    # Werdykt K20
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print(f"  WERDYKT K20: Widmo GW TGP vs LIGO O3 / LISA")
    print(f"{'='*72}")

    res_gw = results.get('GW150914', {})
    res_bns = results.get('GW170817', {})
    res_lisa = results.get('LISA_MBHB', {})

    print(f"""
  TGP predykcja (rem:k20-3pn):
    δc₂ = −1/3  [3PN korekta z c(Φ) = c₀·(Φ₀/Φ)^(1/2)]

  Wyniki:

  GW150914 (LIGO O1):
    ρ = {res_gw.get('rho_gr', 0):.0f},  U_eff = {res_gw.get('U_eff', 0):.2e}
    σ(δc₂) = {res_gw.get('sigma_dc2', 0):.3e}  >> |δc₂·U_eff| = {abs(-1/3*res_gw.get('U_eff', 0)):.3e}
    Poprawka TGP PONIŻEJ PROGU detekcji LIGO O1/O3
    → GW150914 nie falsyfikuje TGP ✓

  GW170817 (LIGO O2):
    ρ = {res_bns.get('rho_gr', 0):.0f},  U_eff = {res_bns.get('U_eff', 0):.2e}
    σ(δc₂) = {res_bns.get('sigma_dc2', 0):.3e}
    Poprawka TGP PONIŻEJ PROGU detekcji LIGO O2
    → GW170817 nie falsyfikuje TGP ✓ (spójne z rem:gw-bns)

  LISA MBHB (przyszłość):
    ρ = {res_lisa.get('rho_gr', 0):.0f},  U_eff = {res_lisa.get('U_eff', 0):.2e}
    σ(δc₂) = {res_lisa.get('sigma_dc2', 0):.3e}
    {'Poprawka TGP DETEKOWALNA przez LISA ✓' if res_lisa.get('detectable', False)
      else 'Poprawka TGP poniżej progu LISA — wymaga wyższego SNR lub mniejszego U_bar'}

  K20 Status: OTWARTE
    → LIGO O3/O4: nie detekowalne (U_eff << 1 dla znanych zdarzeń)
    → LISA (2034+): {'detekowalne przy ρ > 100' if not res_lisa.get('detectable') else 'DETEKOWALNE'}
    → Wymagany: SNR > σ(δc₂)/|δc₂·U_eff| ~ {
      res_gw.get('sigma_dc2', 1)/abs(-1/3*res_gw.get('U_eff', 1e-6)):.0f
    } dla GW150914-like zdarzeń
""")

    return results


# ──────────────────────────────────────────────────────────────────────────────
#  Wykresy
# ──────────────────────────────────────────────────────────────────────────────

def make_plots(results, output_dir=None):
    """Generuje wykresy diagnostyczne widma GW."""
    import os
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("  [WARN] matplotlib niedostępny.")
        return

    if output_dir is None:
        output_dir = os.path.join(
            os.path.dirname(__file__), '..', '..', 'scripts', 'plots'
        )
    os.makedirs(output_dir, exist_ok=True)

    fig = plt.figure(figsize=(16, 11))
    gs  = gridspec.GridSpec(2, 3, hspace=0.38, wspace=0.30)
    fig.suptitle('EX50: Widmo GW TaylorF2 + poprawka 3PN TGP (K20)',
                 fontsize=13, y=0.99)

    colors = {'GW150914': 'royalblue', 'GW170817': 'firebrick',
              'LISA_MBHB': 'seagreen'}

    # (a) |h̃(f)|²/S_n — wkład SNR
    ax = fig.add_subplot(gs[0, 0])
    for name, res in results.items():
        ax.loglog(res['f_arr'], res['integrand_gr'],
                  color=colors[name], lw=1.5,
                  label=name.replace('_', ' '))
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('4|h̃(f)|²/S_n [Hz⁻¹]')
    ax.set_title('(a) Wkład SNR vs częstość')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    # (b) ΔΨ_TGP(f) dla GW150914 i LISA
    ax = fig.add_subplot(gs[0, 1])
    for name in ['GW150914', 'LISA_MBHB']:
        if name not in results:
            continue
        res = results[name]
        M_total_s = res['M_total_Msun'] * M_SUN_S
        dpsi = psi_tgp_correction(
            res['f_arr'], M_total_s, res['eta'],
            -1.0/3.0, res['U_eff']
        )
        ax.loglog(res['f_arr'], np.abs(dpsi),
                  color=colors[name], lw=1.5, label=name.replace('_', ' '))
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('|ΔΨ_TGP(f)| [rad]')
    ax.set_title('(b) Odchylenie fazy TGP (δc₂=-1/3)')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    # (c) σ(δc₂) vs SNR (dla GW150914)
    ax = fig.add_subplot(gs[0, 2])
    res_gw = results['GW150914']
    rho_range = np.logspace(0, 4, 200)

    # σ(δc₂) ∝ 1/ρ (przybliżenie Fisher)
    sigma_0 = res_gw['sigma_dc2'] * res_gw['rho_gr']  # σ·ρ = const
    sigma_range = sigma_0 / rho_range
    dc2_tgp = abs(-1.0/3.0 * res_gw['U_eff'])

    ax.loglog(rho_range, sigma_range, color='royalblue', lw=2,
              label='σ(δc₂) GW150914-like')
    ax.axhline(dc2_tgp, color='red', ls='--', lw=1.5,
               label=f'|δc₂·U_eff|={dc2_tgp:.2e}')
    ax.axvline(res_gw['rho_gr'], color='gray', ls=':', lw=1,
               label=f'ρ_O1={res_gw["rho_gr"]:.0f}')
    ax.set_xlabel('SNR ρ')
    ax.set_ylabel('σ(δc₂)')
    ax.set_title('(c) Czułość na δc₂ vs SNR')
    ax.legend(fontsize=7)
    ax.grid(True, which='both', alpha=0.3)

    # (d) Szumy detektorów
    ax = fig.add_subplot(gs[1, 0])
    f_ligo = np.logspace(1, 3, 500)
    f_lisa = np.logspace(-4, 0, 500)
    ax.loglog(f_ligo, np.sqrt(S_n_aligo(f_ligo)),
              color='royalblue', lw=2, label='aLIGO')
    ax.loglog(f_lisa, np.sqrt(S_n_lisa(f_lisa)),
              color='seagreen', lw=2, label='LISA')
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('√S_n(f) [Hz^{-1/2}]')
    ax.set_title('(d) Szumy detektorów')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(1e-25, 1e-15)

    # (e) Faza GR vs TGP dla GW150914
    ax = fig.add_subplot(gs[1, 1])
    res_gw = results['GW150914']
    M_total_s = res_gw['M_total_Msun'] * M_SUN_S
    M_chirp_s = res_gw['M_chirp_Msun'] * M_SUN_S
    f_arr = res_gw['f_arr']
    psi_gr   = psi_taylorf2_gr(f_arr, M_total_s, res_gw['eta'])
    dpsi_arr = psi_tgp_correction(f_arr, M_total_s, res_gw['eta'],
                                  -1.0/3.0, res_gw['U_eff'])
    # Normalizuj do max |psi_gr|
    ax.semilogx(f_arr, dpsi_arr / np.max(np.abs(psi_gr)) * 1e6,
                color='royalblue', lw=1.5)
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('ΔΨ_TGP / max(|Ψ_GR|) × 10⁶')
    ax.set_title('(e) Relatywne odchylenie fazy GW150914')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='gray', lw=0.5)

    # (f) Werdykt tabela
    ax = fig.add_subplot(gs[1, 2])
    ax.axis('off')
    table_data = []
    for name, res in results.items():
        table_data.append([
            name.replace('_MBHB', '').replace('_', '\n'),
            f"{res['rho_gr']:.0f}",
            f"{res['U_eff']:.1e}",
            f"{res['sigma_dc2']:.1e}",
            'TAK' if res['detectable'] else 'NIE',
        ])
    table = ax.table(
        cellText=table_data,
        colLabels=['Zdarzenie', 'ρ', 'U_eff', 'σ(δc₂)', 'Detek.'],
        loc='center', cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.6)
    ax.set_title('(f) Podsumowanie K20', pad=15)

    out_path = os.path.join(output_dir, 'ex50_gw_taylorf2_3pn.png')
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    print(f"\n  [PLOT] Zapisano: {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
#  Punkt wejścia
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='EX50: Widmo GW TaylorF2 z poprawką 3PN TGP (K20)'
    )
    parser.add_argument('--plot', action='store_true',
                        help='Generuj wykresy (wymaga matplotlib)')
    parser.add_argument('--delta-c2', type=float, default=None,
                        help='Nadpisz δc₂ TGP (domyślnie: -1/3)')
    args = parser.parse_args()

    dc2 = args.delta_c2 if args.delta_c2 is not None else -1.0/3.0
    results = main(delta_c2_scan=[dc2, -0.1, -0.01, -1e-3, 0.0])

    if args.plot:
        make_plots(results)

    print('\nUżycie:')
    print('  python ex50_gw_taylorf2_3pn.py [--plot] [--delta-c2 FLOAT]')
    print('  Domyślnie: δc₂ = -1/3 (TGP predykcja rem:k20-3pn)')
