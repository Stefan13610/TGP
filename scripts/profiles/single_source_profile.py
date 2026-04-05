"""single_source_profile.py
TGP v1 -- Profil pojedynczego zrodla Phi(r)

Weryfikuje wyniki analityczne z sek08 ssec:profil-Phi3:
  prop:Phi3-reduction    -- v=chi^3 eliminuje czlon gradientowy dokladnie
  prop:near-field-series -- chi ~ a0/r - 1 + O(r),  a0=sqrt(2/gamma), b=-1
  prop:composite-profile -- chi_glob = [K^3/r^3 + 1 + 3K*e^{-mr}/r]^{1/3}
  prop:two-regimes-profile -- slabe vs silne zrodlo
  lem:inflection         -- punkt przegiecania chi''=0

Bezwymiarowe: m=sqrt(gamma)=1, Phi0=1, beta=gamma=1.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

GAMMA = 1.0
BETA  = GAMMA
M_SP  = np.sqrt(GAMMA)        # masa pola Phi (zakres Yukawa = 1)
a0    = np.sqrt(2.0 / GAMMA)  # koeff. nasycenia bliskiego pola = sqrt(2)
b_coeff = -1.0                 # koeff. stalej w szeregu bliskim pola (prop)

results = []

def record(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail else ""))

# ── ODE na chi = Phi/Phi0 ─────────────────────────────────────────────────────
def rhs_chi(r, y, beta=BETA, gamma=GAMMA):
    chi, chip = y
    chi = max(chi, 1e-30)
    chipp = (-2.0*chip/r - 2.0*chip**2/chi - beta*chi**2 + gamma*chi**3)
    return [chip, chipp]

# ── ODE na v = chi^3  (prop:Phi3-reduction, eq:v-ode) ────────────────────────
def rhs_v(r, v_vec, beta=BETA, gamma=GAMMA):
    v, vp = v_vec
    v = max(v, 1e-30)
    vpp = -2.0*vp/r - 3.0*beta*(v**(4./3.) - v**(5./3.))
    return [vp, vpp]

# ── Kompozyt v=K^3/r^3 + 1 + 3K*exp(-mr)/r  (prop:composite-profile) ────────
def chi_composite(r, K, m=M_SP):
    """v = K^3/r^3 + 1 + 3K*exp(-mr)/r; chi = v^(1/3)."""
    v = K**3 / r**3 + 1.0 + 3.0*K*np.exp(-m*r)/r
    return np.maximum(v, 1e-30)**(1./3.)

# ── Waruki poczatkowe dalekiego pola Yukawa ───────────────────────────────────
def yukawa_ic(r, K, m=M_SP):
    chi  = 1.0 + K*np.exp(-m*r)/r
    chip = -(1./r**2 + m/r)*K*np.exp(-m*r)
    return chi, chip

# =============================================================================
print("\n=== TEST 1: Redukcja Phi^3 (prop:Phi3-reduction) ===")
# Integrujemy od r0=10 do rf=0.3 i porownujemy chi-ODE z v-ODE.
K1 = 0.05
r0, rf = 10.0, 0.3
chi0, chip0 = yukawa_ic(r0, K1)

sol_chi = solve_ivp(rhs_chi, [r0, rf], [chi0, chip0],
                    dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.01)
v0_ic = chi0**3;  vp0_ic = 3.0*chi0**2*chip0
sol_v  = solve_ivp(rhs_v,  [r0, rf], [v0_ic, vp0_ic],
                   dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.01)

r_chk = np.linspace(rf + 0.05, r0 - 0.5, 400)
chi_c  = sol_chi.sol(r_chk)[0]
v_vals = sol_v.sol(r_chk)[0]
chi_v  = np.where(v_vals > 0, v_vals**(1./3.), np.nan)

max_diff = np.nanmax(np.abs(chi_c - chi_v))
record("T1a: chi-ODE i v-ODE daja ten sam profil chi(r)",
       max_diff < 1e-7, f"max|chi_chi - chi_v| = {max_diff:.2e}")

# =============================================================================
print("\n=== TEST 2: Szereg bliskiego pola (prop:near-field-series) ===")
# Weryfikacja algebraiczna bilansu dominujacego (bez niestabilnej integracji).

# T2a: a0 = sqrt(2/gamma) z bilansu rzedu r^{-5}
res_r5 = abs(6*a0**3 - 3*BETA*a0**5)
record("T2a: Bilans r^{-5}: 6a0^3 = 3*beta*a0^5  => a0=sqrt(2/gamma)",
       res_r5 < 1e-10, f"residual = {res_r5:.2e}")

# T2b: b = -1 z bilansu rzedu r^{-4}
b_from_balance = -BETA * a0**2 / 2.0
res_r4 = abs(b_from_balance - b_coeff)
record("T2b: Bilans r^{-4}: b = -beta*a0^2/2 = -1",
       res_r4 < 1e-12, f"b = {b_from_balance:.15f}")

# T2c: Subdominant term v^{4/3} vs v^{5/3} przy bardzo malym r
# Dominujacy bilans: v'' + (2/r)v' = 3*beta*v^{5/3}
# Poprawka subdominujaca: -3*beta*v^{4/3}; jej wzgledna wielkosc = v^{4/3}/v^{5/3} = 1/chi
r_eps    = 0.0001
chi_eps  = a0/r_eps + b_coeff          # chi wg szeregu bliskiego pola
v_eps    = chi_eps**3
ratio_sub = v_eps**(4./3.) / v_eps**(5./3.)   # = 1/chi_eps << 1
record(f"T2c: Subdominant v^{{4/3}}/v^{{5/3}} = 1/chi przy r={r_eps}: ratio={ratio_sub:.3e}",
       ratio_sub < 1e-2)

# =============================================================================
print("\n=== TEST 3: Daleki zasieg Yukawa (eq:far-field) ===")

K_weak = 0.05
r0_3 = 15.0
chi0_3, chip0_3 = yukawa_ic(r0_3, K_weak)
sol_yuk = solve_ivp(rhs_chi, [r0_3, 0.5], [chi0_3, chip0_3],
                    dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.02)

r_yuk = np.linspace(1.0, 12.0, 400)
chi_num = sol_yuk.sol(r_yuk)[0]
chi_yuk = 1.0 + K_weak*np.exp(-M_SP*r_yuk)/r_yuk
max_err = np.max(np.abs(chi_num - chi_yuk))
record("T3: Profil Yukawa dla slabego zrodla K=0.05",
       max_err < 1e-3, f"max|Delta_chi| = {max_err:.2e}")

# =============================================================================
print("\n=== TEST 4: Profil kompozytowy chi_glob (prop:composite-profile) ===")

# Limit r->0: chi_glob -> K/r
K4 = 0.5
r_tiny = np.array([0.001, 0.005, 0.01, 0.05])
for rt in r_tiny:
    chi_g   = chi_composite(rt, K4)
    chi_exp = K4 / rt
    err = abs(chi_g - chi_exp) / chi_exp
    record(f"T4a: Limit r->0 przy r={rt:.3f}: |Delta|/chi = {err:.2e}", err < 0.02)

# Limit r->inf: chi_glob -> 1 + K*exp(-mr)/r
r_large = np.array([5.0, 8.0, 12.0])
for rl in r_large:
    chi_g   = chi_composite(rl, K4)
    chi_exp = 1.0 + K4*np.exp(-M_SP*rl)/rl
    err = abs(chi_g - chi_exp) / max(chi_exp, 1e-10)
    record(f"T4b: Limit r->inf przy r={rl:.1f}: |Delta|/chi = {err:.2e}", err < 5e-4)

# Porownanie kompozytu z profilem numerycznym (slabe K)
r_comp = np.linspace(0.1, 10.0, 400)
chi_comp_num  = sol_yuk.sol(r_comp)[0]
chi_comp_anal = chi_composite(r_comp, K_weak)
mask = (r_comp >= 0.5) & (r_comp <= 10.0)
err_c = np.max(np.abs(chi_comp_num[mask] - chi_comp_anal[mask]))
record("T4c: Kompozyt vs numeryczny (slabe zrodlo, r in [0.5,10])",
       err_c < 5e-3, f"max|Delta| = {err_c:.2e}")

# =============================================================================
print("\n=== TEST 5: Punkt przegiecania (lem:inflection) ===")

# Uzywamy umiarkowanego K (nie za duzego, by integracja dotarla do malego r)
K5 = 1.0
r0_5 = 12.0
chi0_5, chip0_5 = yukawa_ic(r0_5, K5)
sol_infl = solve_ivp(rhs_chi, [r0_5, 0.1], [chi0_5, chip0_5],
                     dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.005)

# Skanujemy laplasjan sferyczny nabla^2 chi = chi'' + (2/r)*chi'
# wg lematu: nabla^2 chi < 0 blisko zrodla (czlon gamma*chi^3 dominuje),
#            nabla^2 chi > 0 daleko (liniowy rezim Yukawa: nabla^2 chi = m^2*delta_chi > 0)
# => zmiana znaku w pewnym r_infl
r_scan = np.linspace(0.12, r0_5 - 0.2, 4000)
r_scan = r_scan[(r_scan >= sol_infl.t[-1]) & (r_scan <= sol_infl.t[0])]
if len(r_scan) > 10:
    chi_s   = sol_infl.sol(r_scan)[0]
    chip_s  = sol_infl.sol(r_scan)[1]
    chipp_s = np.array([rhs_chi(r, [c, cp])[1]
                        for r, c, cp in zip(r_scan, chi_s, chip_s)])
    # Laplasjan sferyczny
    lap_s   = chipp_s + 2.0*chip_s/r_scan
    sign_chg = np.where(np.diff(np.sign(lap_s)))[0]
    if len(sign_chg) > 0:
        r_infl_est = r_scan[sign_chg[0]]
        record("T5: Punkt przegiecania nabla^2 chi=0 istnieje (lem:inflection)",
               True, f"r_infl ~ {r_infl_est:.3f}")
    else:
        sign_edge = (lap_s[0] < 0) != (lap_s[-1] < 0)
        record("T5: Punkt przegiecania nabla^2 chi=0 istnieje (lem:inflection)",
               sign_edge, f"lap(r_min)={lap_s[0]:.2f}, lap(r_max)={lap_s[-1]:.2f}")
else:
    record("T5: Punkt przegiecania nabla^2 chi=0 istnieje (lem:inflection)",
           False, "integracja zbyt krotka")

# =============================================================================
print("\n=== TEST 6: Dwa rezimy zrodla (prop:two-regimes-profile) ===")

# Slabe zrodlo K << a0: kompozyt ~ Yukawa globalnie
K6_weak = 0.005
r_t6 = np.linspace(0.5, 8.0, 300)
chi_comp6 = chi_composite(r_t6, K6_weak)
chi_yuk6  = 1.0 + K6_weak*np.exp(-M_SP*r_t6)/r_t6
max_nl6 = np.max(np.abs(chi_comp6 - chi_yuk6) / chi_yuk6)
record("T6a: Slabe zrodlo (K=0.005<<a0): kompozyt ~ Yukawa",
       max_nl6 < 1e-4, f"max nl = {max_nl6:.2e}")

# Silne zrodlo K >> a0: limit r->0 daje chi ~ K/r (nie a0/r jeszcze)
K6_strong = 5.0*a0
for rs in [0.01, 0.02, 0.05]:
    chi_g  = chi_composite(rs, K6_strong)
    chi_Kr = K6_strong / rs     # limit liniowy dla tego zrodla
    chi_a0r= a0 / rs            # limit nasycenia
    # chi_g powinno byc bliskie max(K/r, a0/r) = K/r (bo K>a0)
    err = abs(chi_g - chi_Kr) / chi_Kr
    record(f"T6b: Silne zrodlo chi->K/r przy r={rs:.3f}: err={err:.3e}", err < 0.02)

# =============================================================================
print("\n=== TEST 7: Proznia chi=1 ===")
dy = rhs_chi(2.0, [1.0, 0.0])
record("T7: chi=1 jest dok. rozwiazaniem prozniowym (chi''=0)",
       abs(dy[1]) < 1e-14, f"chi'' = {dy[1]:.2e}")

# =============================================================================
print("\n" + "="*60)
n_pass  = sum(1 for _, s, _ in results if s == "PASS")
n_fail  = sum(1 for _, s, _ in results if s == "FAIL")
n_total = len(results)
print(f"WYNIK: {n_pass}/{n_total} PASS,  {n_fail} FAIL")
if n_fail == 0:
    print("Wszystkie testy PASS.")
    print("prop:Phi3-reduction, prop:near-field-series,")
    print("prop:composite-profile, lem:inflection zweryfikowane.")
print("="*60)
