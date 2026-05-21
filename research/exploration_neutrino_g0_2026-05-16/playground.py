"""
NEUTRINO PLAYGROUND — TGP universal formula extrapolation to low-g_0
=====================================================================

NIE jest to formalny cykl. Brak Phase 0 gate, brak BINDING contract.
YOLO exploration: szybka ekstrapolacja universal formuły m=c·A²·g_0^(e²/2)
z lepton calibration na neutrino mass range.

PRE-DISCLAIMERS (uczciwie):
  1. A_tail(g_0) power-law z 3-point lepton fit — ekstrapolacja 8× w dół
  2. Cross-section heurystyka σ ~ A_tail⁴ — naiwna, nie rigorous Feynman
  3. ω_C_ratio² jako suppression — heuristic argument, nie derivation
  4. Wartości m_ν użyte: oscylacyjne best-fit + cosmology bound
  5. Power-law przedłużone w region g_0 < 0.5 jest pure extrapolation
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import sympy as sp
from sympy import Float, log, exp, sqrt, E, symbols, nsolve

# =========================================================================
# Lepton calibration (LIVE — z why_n3 PHASE2_n_alpha_derivation 2026-05-01)
# =========================================================================
LEPTON_DATA = {
    'e':  (Float('0.86941'), Float('0.11003')),
    'mu': (Float('1.40673'), Float('0.65041')),
    'tau':(Float('1.75505'), Float('1.66645')),
}
g0_e, A_e = LEPTON_DATA['e']
e2_half = E**2 / 2  # ≈ 3.6945

# Power-law fit z poprzedniego cyklu: A_tail = exp(b)·g_0^a, R²=0.999
slope_a = Float('3.8404')
intercept_b = Float('-1.6866')

def A_tail(g0):
    return exp(intercept_b) * g0**slope_a

def m_over_me(g0):
    return (A_tail(g0)/A_e)**2 * (g0/g0_e)**(e2_half)

def psi_of_g(g0):
    """ψ = 0.3814·g + 0.6186 (Faza 1 R3→M9.1'' reparametryzacja)"""
    return Float('0.3814')*g0 + Float('0.6186')

def A_metric(psi):
    """A(ψ) = (4-3ψ)/ψ (M9.1'' time-component)"""
    return (4 - 3*psi) / psi

def solve_g0_for_mass(target_ratio, guess=0.3):
    """Find g_0 such that m/m_e = target_ratio"""
    g = symbols('g', positive=True)
    try:
        sol = nsolve(m_over_me(g) - target_ratio, g, guess, prec=20)
        return float(sol)
    except Exception:
        return None

m_e_MeV = Float('0.5109989461')

# =========================================================================
# Neutrino mass scenarios — PDG oscillation + cosmology
# =========================================================================
# Δm²_21 = 7.42e-5 eV²; Δm²_32 = 2.51e-3 eV² (PDG 2024 normal hierarchy)
# Cosmology Σm_ν < 0.12 eV (Planck 2018)
# Pick lightest m_1 = 0.005 eV (minimum sum hierarchy compatible)

m1_eV = 0.005
m2_eV = float(sqrt(Float(m1_eV)**2 + Float('7.42e-5')).evalf())
m3_eV = float(sqrt(Float(m2_eV)**2 + Float('2.51e-3')).evalf())

# Also: degenerate hierarchy candidate (m_1 ~ 0.05 eV)
m1_deg = 0.05
m2_deg = float(sqrt(Float(m1_deg)**2 + Float('7.42e-5')).evalf())
m3_deg = float(sqrt(Float(m2_deg)**2 + Float('2.51e-3')).evalf())

# =========================================================================
# Sekcja 1: Lepton calibration ground truth
# =========================================================================
print("="*88)
print("NEUTRINO PLAYGROUND — TGP universal formula w low-g_0 range")
print("="*88)
print()
print("LEPTON CALIBRATION (LIVE z why_n3 PHASE2 2026-05-01):")
print(f"  {'particle':10s} {'g_0':>10s} {'A_tail':>10s} {'psi':>10s} {'A_metric':>10s} {'c_loc/c_0':>10s}")
print("  " + "-"*72)
for label, (g, A) in LEPTON_DATA.items():
    psi_sym = psi_of_g(g)
    psi = float(psi_sym.evalf())
    Am = float(A_metric(psi_sym).evalf())
    c_loc = float(sqrt(A_metric(psi_sym)).evalf())
    print(f"  {label:10s} {float(g):>10.4f} {float(A):>10.4f} {psi:>10.4f} {Am:>10.4f} {c_loc:>10.4f}")
print()
print(f"  Power-law fit: A_tail = {float(exp(intercept_b)):.4f} · g_0^{float(slope_a):.4f}  (R²=0.9988)")
print(f"  Exponent e²/2 = {float(e2_half):.4f}")
print(f"  R3 Lorentzian horizon: psi = 4/3 = 1.3333")
print()

# =========================================================================
# Sekcja 2: Neutrino predictions (light hierarchy + degenerate + sterile)
# =========================================================================
print("="*88)
print("NEUTRINO SCENARIOS — TGP universal formula prediction")
print("="*88)
print()

scenarios = [
    ('nu_1 (light)',  m1_eV),
    ('nu_2 (light)',  m2_eV),
    ('nu_3 (light)',  m3_eV),
    ('Sum(light)',    m1_eV + m2_eV + m3_eV),
    ('nu_1 (deg)',    m1_deg),
    ('nu_2 (deg)',    m2_deg),
    ('nu_3 (deg)',    m3_deg),
    ('Sum(deg)',      m1_deg + m2_deg + m3_deg),
    ('sterile (2eV)', 2.0),
]

print(f"  {'particle':16s} {'m [eV]':>10s} {'g_0_nu':>9s} {'A_tail':>11s} {'A_t/A_e':>11s} {'psi_nu':>9s} {'A(psi)':>9s} {'c_loc/c0':>9s}")
print("  " + "-"*102)

g0_results = {}
for name, m_eV in scenarios:
    target = (m_eV * 1e-6) / float(m_e_MeV)
    # Try multiple initial guesses
    g0 = None
    for guess in [0.1, 0.2, 0.3, 0.5, 0.7]:
        g0 = solve_g0_for_mass(target, guess)
        if g0 is not None and 0 < g0 < 1.0:
            break
    if g0 is None or g0 <= 0:
        print(f"  {name:16s} {m_eV:>10.4g} {'-- no solution --':>9s}")
        continue

    g0_results[name] = g0
    g0_sym = Float(g0)
    A_t = float(A_tail(g0_sym).evalf())
    psi_sym = psi_of_g(g0_sym)
    psi = float(psi_sym.evalf())
    Am = float(A_metric(psi_sym).evalf())
    c_loc = float(sqrt(Float(Am)).evalf()) if Am > 0 else float('nan')
    A_ratio = A_t / float(A_e)
    print(f"  {name:16s} {m_eV:>10.4g} {g0:>9.5f} {A_t:>11.3e} {A_ratio:>11.3e} {psi:>9.5f} {Am:>9.4f} {c_loc:>9.4f}")

print()

# =========================================================================
# Sekcja 3: Suppression mechanism analysis
# =========================================================================
print("="*88)
print("SUPPRESSION MECHANISMS — heurystyczne candidate (sigma_nu vs sigma_e)")
print("="*88)
print()
print("  Mechanism A: sigma ~ A_tail^4  (Yukawa coupling y ~ A_tail²)")
print("  Mechanism B: (omega_C_nu/omega_C_e)^2  (Compton frequency mismatch)")
print("  Reference:  sigma_nu_e_real / sigma_e_e_Thomson ~ 10^(-20)")
print()
print(f"  {'particle':16s} {'m [eV]':>10s} {'(A_t/A_e)^4':>14s} {'(omega_ratio)^2':>17s} {'product':>14s} {'log10(prod)':>12s}")
print("  " + "-"*94)

target_required = 1e-20  # observed σ suppression

for name, m_eV in scenarios:
    if name not in g0_results:
        continue
    target = (m_eV * 1e-6) / float(m_e_MeV)
    g0 = g0_results[name]
    A_t = float(A_tail(Float(g0)).evalf())
    A_ratio = A_t / float(A_e)
    A4 = A_ratio**4
    omega_sq = target**2
    product = A4 * omega_sq
    log10p = float(log(Float(product), 10).evalf()) if product > 0 else float('nan')
    print(f"  {name:16s} {m_eV:>10.4g} {A4:>14.3e} {omega_sq:>17.3e} {product:>14.3e} {log10p:>12.2f}")

print()
print(f"  Target log10(suppression) ~ -20  (sigma_nu-e / sigma_e-e Thomson)")
print()

# =========================================================================
# Sekcja 4: phi-ladder check (golden ratio between generations)
# =========================================================================
print("="*88)
print("phi-DRABINKA TEST (g_0^mu / g_0^e = phi = 1.6180; sprawdzić dla ν?)")
print("="*88)
print()

phi_val = float((1 + sqrt(Float(5)))/2)

# Lepton drabinka
g_e_mu = float(LEPTON_DATA['mu'][0] / LEPTON_DATA['e'][0])
g_mu_tau = float(LEPTON_DATA['tau'][0] / LEPTON_DATA['mu'][0])
print(f"  Lepton charged drabinka:")
print(f"    g_0^mu / g_0^e   = {g_e_mu:.5f}  (phi = {phi_val:.5f}, diff = {(g_e_mu - phi_val)*100:+.3f})")
print(f"    g_0^tau / g_0^mu = {g_mu_tau:.5f}  (phi expected? diff vs phi = {(g_mu_tau - phi_val)*100:+.3f})")
print()
print(f"  Neutrino candidates:")
for scenario_set in [('light hierarchy', ['nu_1 (light)', 'nu_2 (light)', 'nu_3 (light)']),
                      ('degenerate hierarchy', ['nu_1 (deg)', 'nu_2 (deg)', 'nu_3 (deg)'])]:
    label, keys = scenario_set
    if not all(k in g0_results for k in keys):
        continue
    g1, g2, g3 = [g0_results[k] for k in keys]
    r_21 = g2/g1
    r_32 = g3/g2
    print(f"    {label}:")
    print(f"      g_0_nu2/g_0_nu1 = {r_21:.5f}  (phi diff = {(r_21 - phi_val)*100:+.3f})")
    print(f"      g_0_nu3/g_0_nu2 = {r_32:.5f}  (phi diff = {(r_32 - phi_val)*100:+.3f})")

print()

# =========================================================================
# Sekcja 5: Cross-relationship: czy m_nu są na drabince z m_e?
# =========================================================================
print("="*88)
print("CROSS-LADDER: g_0_nu vs g_0_e?  (is neutrino on continuation of lepton ladder?)")
print("="*88)
print()
print(f"  g_0_e = 0.8694  (electron anchor)")
print(f"  If phi-ladder continues DOWN: g_0_e / phi = {float(g0_e)/phi_val:.5f}")
print(f"                                g_0_e / phi^2 = {float(g0_e)/phi_val**2:.5f}")
print(f"                                g_0_e / phi^3 = {float(g0_e)/phi_val**3:.5f}")
print(f"                                g_0_e / phi^4 = {float(g0_e)/phi_val**4:.5f}")
print()

# What mass would TGP predict for each ladder step?
print(f"  Predicted m_nu if g_0_nu is on ladder step:")
print(f"  {'step':10s} {'g_0_predicted':>15s} {'m_nu predicted [eV]':>22s}")
print("  " + "-"*52)
for n in [1, 2, 3, 4, 5, 6]:
    g_step = float(g0_e) / phi_val**n
    target = m_over_me(Float(g_step))
    target_f = float(target.evalf())
    m_pred_eV = target_f * float(m_e_MeV) * 1e6
    print(f"  e/phi^{n:1d}    {g_step:>15.5f} {m_pred_eV:>22.4e}")
print()

# =========================================================================
# Sekcja 6: Lorentzian horizon + R3 barrier check
# =========================================================================
print("="*88)
print("LORENTZIAN HORIZON CHECK (psi < 4/3 = 1.333; physical particle constraint)")
print("="*88)
print()
print("  Particles previously verified IN Lorentzian domain (why_n3 Phase 1):")
print(f"    e:  psi = 0.9502  (deep below horizon)")
print(f"    mu: psi = 1.1551  (closer)")
print(f"    tau: psi = 1.2880  (very close to 4/3)")
print(f"    R3 barrier g_0_crit = 1.874 -> psi = 4/3 (horizon)")
print()
print("  Neutrino predictions:")
for name in ['nu_1 (light)', 'nu_2 (light)', 'nu_3 (light)', 'sterile (2eV)']:
    if name not in g0_results:
        continue
    g0 = g0_results[name]
    psi = float(psi_of_g(Float(g0)).evalf())
    status = "OK" if psi < 4/3 else "HORIZON-VIOLATION"
    distance_to_psi_vacuum = abs(psi - 1.0)
    print(f"    {name:18s}: g_0={g0:.4f}, psi={psi:.4f}  [{status}]  |psi-1|={distance_to_psi_vacuum:.4f}")
print()

# =========================================================================
# Sekcja 7: Compton frequency w lokalnym frame
# =========================================================================
print("="*88)
print("COMPTON FREQUENCY w LOKALNYM frame (omega_C = m·c²/hbar)")
print("="*88)
print()
print(f"  omega_C_electron (vacuum) ≈ 7.76e20 Hz (reference)")
print()
print(f"  {'particle':16s} {'m [eV]':>10s} {'omega_C_vac [Hz]':>18s} {'omega_C_local [Hz]':>20s} {'ratio_local/vac':>16s}")
print("  " + "-"*86)
omega_C_e_vac = 7.76e20
hbar = 6.582e-16  # eV·s
c = 3e8  # m/s
for name, m_eV in scenarios:
    if name not in g0_results:
        continue
    g0 = g0_results[name]
    psi_sym = psi_of_g(Float(g0))
    psi = float(psi_sym.evalf())
    Am = float(A_metric(psi_sym).evalf())
    if Am <= 0:
        continue
    omega_C_vac_Hz = m_eV / hbar
    # Local: c_local = c·sqrt(A), m_local "as seen by substrate" — debatable
    # Naive: ω_C_local = m·c_local²/ħ = ω_C_vac · A(ψ)
    omega_C_local_Hz = omega_C_vac_Hz * Am
    ratio = Am  # = omega_local / omega_vac
    print(f"  {name:16s} {m_eV:>10.4g} {omega_C_vac_Hz:>18.3e} {omega_C_local_Hz:>20.3e} {ratio:>16.4f}")
print()
print("  Interpretacja: dla psi < 1 (low-g_0 neutrino), A(psi) > 1 → omega_C_local > omega_C_vac")
print("  Neutrino w swoim lokalnym frame 'tika szybciej' niż w vacuum frame")
print()

print("="*88)
print("END PLAYGROUND")
print("="*88)
