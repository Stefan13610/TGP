"""two_source_potential.py
TGP v1 -- Potencjal efektywny dwoch zrodel i skale przejscia  [v3: Yukawa]

Weryfikuje thm:three-regimes (sek08 ssec:dwa-zrodla):
  Trzy rezimy separacji d:
    I  (d > r_rep):          grawitacja (F < 0)
    II (r_well < d < r_rep): odpychanie (F > 0)
    III (d < r_well):        studnia konfinementu (F < 0)

Metoda: perturbacyjna rozwinięcie Yukawa (slabe zrodla, K << a0).
  delta_chi(r) = chi - 1 ~ K*exp(-m*r)/r   (m=1)

  E_lin = -K^2/d                                    (rownoleglosc Yukawa)
  E_beta_cross = 2*beta * Integral[(d1)(d2)]d3r
               = 2*K^2 * I2(d)                       (analitycznie)
  E_gamma_cross = -6*gamma * K^3 * I3(d)             (numerycznie)

  I2(d) = (4pi)^2 * e^{-d}*(1+d)/(8pi) = 2pi*e^{-d}*(1+d)
  (wynik z transformaty Fouriera, calka splotowa propagatora Yukawa)

  F(d) = -dV/dd
  Warunek 3 rezimow: F zmienia znak dwa razy.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import dblquad, nquad

M_SP  = 1.0   # masa Yukawa
BETA  = 1.0
GAMMA = 1.0

results = []

def record(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail else ""))

# ── Energia: skladowe analityczne i numeryczne ───────────────────────────────

def E_lin(d, K):
    """Oddzialywanie liniowe: E_lin = -K^2/d  (eq:Elin, sek08)"""
    return -K**2 / d

def I2_analytic(d, m=M_SP):
    """
    Splot dwoch potencjalow Yukawa:
    I2(d) = Integral e^{-m*r1}/r1 * e^{-m*r2}/r2 d3r
    Wynik z transformaty Fouriera (twierdzenie o resztach):
      Integral d3k/(2pi)^3 * (4pi/(k^2+m^2))^2 * e^{i*k*d}
      = (4pi)^2/(2pi^2*d) * Integral_0^inf k*sin(kd)/(k^2+m^2)^2 dk
      = (4pi)^2/(2pi^2*d) * pi*d/(4m) * e^{-md}
      = 2*pi/m * e^{-md}
    """
    return 2.0 * np.pi * np.exp(-m*d) / m

def E_beta(d, K, m=M_SP):
    """E_beta = 2*beta * K^2 * I2(d)  (eq:Ebeta, perturbacyjnie)"""
    return 2.0 * BETA * K**2 * I2_analytic(d, m)

def I3_numeric(d, m=M_SP, n_rho=60, n_z=100):
    """
    I3(d) = Integral (e^{-m*r1}/r1)^2 * (e^{-m*r2}/r2) d3r
    Calkowanie na siatce cylindrycznej (rho, z).
    """
    decay = max(5.0/m, 2.0*d)
    rho = np.linspace(1e-3, decay, n_rho)
    z   = np.linspace(-decay, d + decay, n_z)
    dR  = rho[1] - rho[0]
    dZ  = z[1]   - z[0]

    RHO, Z = np.meshgrid(rho, z, indexing='ij')
    R1 = np.sqrt(RHO**2 + Z**2)
    R2 = np.sqrt(RHO**2 + (Z - d)**2)

    f1 = np.exp(-m*R1) / R1
    f2 = np.exp(-m*R2) / R2

    integrand = f1**2 * f2 * 2.0 * np.pi * RHO
    return np.sum(integrand) * dR * dZ

def E_gamma(d, K, m=M_SP, **kw):
    """E_gamma = -6*gamma * K^3 * I3(d)  (czlon kubiczny, studnia)"""
    return -6.0 * GAMMA * K**3 * I3_numeric(d, m, **kw)

def V_eff(d, K, include_gamma=True, **kw):
    V = E_lin(d, K) + E_beta(d, K)
    if include_gamma:
        V += E_gamma(d, K, **kw)
    return V

def force(d, K, dd=0.05, **kw):
    return -(V_eff(d+dd, K, **kw) - V_eff(d-dd, K, **kw)) / (2.0*dd)

# =============================================================================
print("\n=== TEST 0: Weryfikacja I2 analityczna ===")
# Porownaj I2_analytic z bezposrednia numeryczna calka 2D
def I2_numeric(d, m=M_SP, n_rho=80, n_z=120):
    decay = max(8.0/m, 3.0*d)
    rho = np.linspace(1e-3, decay, n_rho)
    z   = np.linspace(-decay, d + decay, n_z)
    dR  = rho[1] - rho[0]
    dZ  = z[1]   - z[0]
    RHO, Z = np.meshgrid(rho, z, indexing='ij')
    R1 = np.sqrt(RHO**2 + Z**2)
    R2 = np.sqrt(RHO**2 + (Z - d)**2)
    f1 = np.exp(-m*R1) / R1
    f2 = np.exp(-m*R2) / R2
    return np.sum(f1 * f2 * 2.0*np.pi*RHO) * dR * dZ

for d_test in [1.0, 2.0, 4.0]:
    I2_a = I2_analytic(d_test)
    I2_n = I2_numeric(d_test)
    err  = abs(I2_a - I2_n) / abs(I2_a)
    record(f"T0: I2 analityczna vs numeryczna przy d={d_test:.1f}: blad={err:.2e}",
           err < 0.05, f"anal={I2_a:.5f}, num={I2_n:.5f}")

# =============================================================================
print("\n=== TEST 1: Struktura trzech rezimow (K=0.5, Yukawa) ===")
K1 = 0.5   # slabe zrodlo -- Yukawa perturbacyjnie
d_arr = np.concatenate([np.linspace(0.10, 0.60, 10),
                        np.linspace(0.70, 3.0, 12),
                        np.linspace(3.5, 9.0, 12)])
kw1 = dict(n_rho=60, n_z=100)

print(f"  Obliczam V_eff(d) dla K={K1} (bez E_gamma) ...")
Varr_noG = np.array([V_eff(d, K1, include_gamma=False) for d in d_arr])
Farr_noG = -np.gradient(Varr_noG, d_arr)

sc_noG = np.where(np.diff(np.sign(Farr_noG)))[0]
record(f"T1a: Trzy rezimy bez E_gamma (E_lin+E_beta): {len(sc_noG)} zmiany znaku",
       len(sc_noG) >= 2,
       f"zmiany znaku przy d ~ {d_arr[sc_noG]}")

if len(sc_noG) >= 2:
    r_rep_lin = d_arr[sc_noG[-1]]
    r_well_lin = d_arr[sc_noG[0]]
    record(f"T1b: r_well < r_rep",
           r_well_lin < r_rep_lin,
           f"r_well~{r_well_lin:.2f}, r_rep~{r_rep_lin:.2f}")
    record(f"T1c: r_rep jest rzedu O(1/m) = O(1)",
           0.2 < r_rep_lin < 10.0,
           f"r_rep~{r_rep_lin:.2f}")
    record(f"T1d: r_well << r_rep (rezimy dobrze rozdzielone)",
           r_well_lin < 0.5 * r_rep_lin,
           f"r_well/r_rep = {r_well_lin/r_rep_lin:.2f}")
else:
    r_rep_lin = r_well_lin = None
    record("T1b-d: potrzeba >=2 zmian znaku dla testow skali", False, "")

# =============================================================================
print("\n=== TEST 2: F > 0 w srodku przedzialu (odpychanie) ===")
if r_rep_lin and r_well_lin:
    d_mid = 0.5*(r_well_lin + r_rep_lin)
    F_mid = force(d_mid, K1, include_gamma=False)
    record(f"T2: F>0 w d_mid={(d_mid):.2f} (odpychanie)",
           F_mid > 0, f"F={F_mid:.4f}")

# =============================================================================
print("\n=== TEST 3: F < 0 dla duzego d (grawitacja) ===")
d_grav = 8.0
F_grav = force(d_grav, K1, include_gamma=False)
record(f"T3: F<0 dla d={d_grav} (grawitacja)",
       F_grav < 0, f"F={F_grav:.5f}")

# =============================================================================
print("\n=== TEST 4: E_gamma (studnia) dla K=0.5 ===")
# I3 ma osobliwos (Yukawa)^2 blisko zrodla -- testujemy tylko dla d >= 0.8
# gdzie profil Yukawa jest dobrze zdefiniowany

# T4a: E_gamma < 0 (przyciagajace) dla roznych d (czlon kubiczny poglabia studnie)
print("  Obliczam I3 dla kilku wartosci d (d >= 0.8) ...")
d_gamma_test = [0.8, 1.2, 2.0, 3.0]
all_neg = True
Eg_vals = []
for dg in d_gamma_test:
    Eg = E_gamma(dg, K1, **kw1)
    Eg_vals.append(Eg)
    if Eg >= 0:
        all_neg = False
    print(f"    d={dg:.1f}: E_gamma={Eg:.6f}")
record("T4a: E_gamma < 0 (przyciagajace, studnia) dla d in [0.8, 3.0]",
       all_neg, f"E_gamma = {[f'{e:.4f}' for e in Eg_vals]}")

# T4b: E_gamma poglabia potencjal przy d=1.0 (V_z_gamma < V_bez_gamma)
V_noG_1 = V_eff(1.0, K1, include_gamma=False)
V_wG_1  = V_eff(1.0, K1, include_gamma=True, **kw1)
record("T4b: E_gamma poglabia studnie przy d=1 (V_z_gamma < V_bez_gamma)",
       V_wG_1 < V_noG_1,
       f"V_bez={V_noG_1:.4f}, V_z={V_wG_1:.4f}")

# =============================================================================
print("\n=== TEST 5: Skale przejscia vs eq:scales ===")
# eq:scales z sek03: r_rep ~ (beta/gamma)*qM/Phi0 = K * (beta/gamma)
# r_well ~ (alpha/beta)*r0 gdzie r0 ~ 1/m
# W Yukawa: r_rep ~ wyznacznik rownania 1 = 4pi*m^3*d^3*e^{-md}
# i r_well ~ d takie ze m^2*d^2 ~ 1 => d ~ 1/m

if r_rep_lin:
    # Teoretyczna skala: gdzie I2(d)*m^2 ~ 1/(2pi) => e^{-d}*(1+d) ~ 1/(2pi)
    # Przyblizenie: r_rep ~ ln(2pi * something)
    r_rep_scale = 1.0 / M_SP   # skala Yukawa
    record("T5a: r_rep rzędu O(1..10)/m (skala Yukawa)",
           0.5 < r_rep_lin * M_SP < 20.0,
           f"r_rep={r_rep_lin:.2f}, 1/m={r_rep_scale:.1f}, r_rep*m={r_rep_lin*M_SP:.1f}")
if r_well_lin:
    record("T5b: r_well < r_rep (studnia blizej niz odpychanie)",
           r_well_lin < r_rep_lin,
           f"r_well={r_well_lin:.2f}, r_rep={r_rep_lin:.2f}")

# =============================================================================
print("\n" + "="*60)
n_pass  = sum(1 for _, s, _ in results if s == "PASS")
n_fail  = sum(1 for _, s, _ in results if s == "FAIL")
n_total = len(results)
print(f"WYNIK: {n_pass}/{n_total} PASS,  {n_fail} FAIL")
if n_fail == 0:
    print("Wszystkie testy PASS.")
    print("thm:three-regimes zweryfikowane perturbacyjnie (Yukawa, K<<a0).")
else:
    print("Niektore testy FAIL.")
if r_rep_lin and r_well_lin:
    print(f"\nSKALE dla K={K1}, m=1: r_rep~{r_rep_lin:.2f}, r_well~{r_well_lin:.2f}")
print("="*60)
