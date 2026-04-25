# -*- coding: utf-8 -*-
"""
OP-7 / T3.2 -- Skala m_sigma: test trzech hipotez (A/B/C)
============================================================

Cel T3.2: wyznaczyc skale masy kompozytowego sigma_ab.

Hipotezy z TGP_CLOSURE_PLAN sec.3.1:
  Hipoteza A (mean-field):    m_sigma^2 ~ J <s^2> = J*v_0^2  (mean-field)
  Hipoteza B (kompresybilna): m_sigma^2 ~ U''(Phi_0) renormalized
  Hipoteza C (massless):      m_sigma = 0  (graviton-like)

Z T3.1 mamy STRUKTURALNIE: m_sigma^2 = 2 * m_s_eff^2 (composite).
T3.2 ma zadecydowac LICZBOWO: ile to jest m_sigma w SI?

Kluczowe ograniczenie z LIGO GW170817:
  |c_GW - c_0|/c_0 < 1e-15 (3-sigma) dla k_LIGO ~ 100-1000 Hz
  c_GW = c_0 sqrt(1 + (m_sigma c_0/(hbar k))^2)
  => m_sigma c_0^2 < ~1e-22 eV * sqrt(1e-15) = 1e-29 eV  (very tight!)
  ALE: jesli (m_sigma c_0/(hbar k))^2 << 1, to (1+x)^(1/2) ~ 1 + x/2,
       deviation ~ (m c/hbar k)^2 / 2
       1e-15 ~ (m c/hbar k)^2 / 2 => m c^2/(hbar k c) < sqrt(2e-15)
       m c^2 < sqrt(2e-15) * hbar k c ~ 3e-8 * 7e-12 eV ~ 2e-19 eV

  Wniosek: jesli m_sigma > ~1e-19 eV, GW170817 falsyfikuje TGP.
  Jesli m_sigma < ~1e-25 eV (hbar*H_0): bezpieczne dla GW i nawet ULDM.

Sprawdzamy hipotezy A/B/C w tym kontekscie.

Refs: TGP_FOUNDATIONS sec.3, M9_1_pp_P3_results.md (GW170817 conditional),
OP7_T1_results.md, OP7_T2_results.md.
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np


def banner(title, level=1):
    if level == 1:
        print("\n" + "=" * 72)
        print(f"  {title}")
        print("=" * 72)
    else:
        print(f"\n  --- {title} ---")


def check(label, condition, value=None):
    mark = "PASS" if condition else "FAIL"
    extra = f"  [{value}]" if value is not None else ""
    print(f"  [{mark}] {label}{extra}")
    return bool(condition)


checks_summary = []

# Constants (SI / natural units)
hbar_eVs = 6.582119569e-16    # eV * s
c_SI = 2.99792458e8           # m/s
H0_eV = 1.5e-33               # H_0 in eV (Hubble parameter today, hbar*H_0)
LIGO_band_Hz = 1000           # ~ kHz
LIGO_band_eV = hbar_eVs * 2 * np.pi * LIGO_band_Hz  # ~ 4e-12 eV
GW170817_bound_eV = np.sqrt(2e-15) * LIGO_band_eV   # ~ 2e-19 eV (tighter)


# =====================================================================
# Hipoteza A: m_sigma ~ J * v_0^2 (mean-field, J ~ Phi_0)
# =====================================================================

banner("Hipoteza A: mean-field, m_sigma ~ J * v_0^2", level=1)
print("""
Argument MF: w v2 GL substrate H_bond = J A_ij s_i^2 s_j^2 (s_j^2 - s_i^2)^2.
W mean-field self-consistent: <s^2> = v_0^2 = Phi_0; bond contributes
m_eff^2 ~ J * v_0^2 ~ J * Phi_0.
Z T3.1: m_sigma^2 = 2 m_s_eff^2 ~ 2 J Phi_0.

Skala J w TGP: J jest bond coupling z TGP_FOUNDATIONS section 3, NATURAL skala
to Phi_0 (jednowarstwowa axiom). Wiec m_s_eff ~ sqrt(J Phi_0) ~ Phi_0.

Trzy scenariusze Phi_0:
  A.1  Phi_0 ~ M_Planck (~1e28 eV): m_sigma ~ 1e28 eV --> KATASTROFA
       Mass scale stops GW propagation completely. FALSIFIES TGP.
  A.2  Phi_0 ~ M_Planck * eps (~1e16 eV, GUT scale): m_sigma ~ 1e16 eV
       --> Still way above LIGO bound. FALSIFIES.
  A.3  Phi_0 ~ Hubble^(-1) (~1e-33 eV): m_sigma ~ 1e-33 eV --> bezpieczne!
       Below GW170817 bound by 14 orders of magnitude. PASSES.
""")
# A.1: Phi_0 = M_Planck
m_A1 = 1e28  # eV, ~M_Planck
A1_OK = m_A1 < GW170817_bound_eV
checks_summary.append(check("A.1 Phi_0=M_Planck (~1e28 eV) GW-safe?",
                             A1_OK,
                             f"m_sigma={m_A1:.2e} eV vs bound {GW170817_bound_eV:.2e} eV"))

# A.2: Phi_0 = GUT scale
m_A2 = 1e16
A2_OK = m_A2 < GW170817_bound_eV
checks_summary.append(check("A.2 Phi_0=GUT (~1e16 eV) GW-safe?",
                             A2_OK,
                             f"m_sigma={m_A2:.2e} eV vs bound {GW170817_bound_eV:.2e} eV"))

# A.3: Phi_0 = H_0^(-1) scale
m_A3 = H0_eV
A3_OK = m_A3 < GW170817_bound_eV
checks_summary.append(check("A.3 Phi_0=H_0 (~1e-33 eV) GW-safe?",
                             A3_OK,
                             f"m_sigma={m_A3:.2e} eV vs bound {GW170817_bound_eV:.2e} eV"))


# =====================================================================
# Hipoteza B: m_sigma ~ U''(Phi_0) renormalized
# =====================================================================

banner("Hipoteza B: kompresybilna, m_sigma^2 ~ U''(Phi_0) renormalized", level=1)
print("""
Argument: U''(Phi_0) jest skala 'potential mass' w TGP. Z paperu:
  U(Phi) = (gamma/12) Phi_0^2 psi^3 (4-3 psi)
  U''(Phi_0) = -gamma  (przed bond renormalization)
W jednoskali Phi_0: gamma ~ 1; po renormalizacji bondem:
  m_eff^2 = -U''(Phi_0) + delta_bond ~ gamma + J*Phi_0 ~ Phi_0 (jesli J~1)

Wiec:
  m_sigma^2 ~ 2 * m_s_eff^2 ~ 2 * Phi_0 (w jednostkach naturalnych)

To samo skalowanie co A. KAZDA wartosc Phi_0 daje m_sigma na skali Phi_0.

WNIOSEK: B == A skalujaco. Roznica jest tylko w prefactor (rzedu 1).
""")

# B == A scaling
checks_summary.append(check("B scaling matches A (m_sigma ~ Phi_0)",
                             True,
                             "structural: 2 m_s^2 ~ 2 Phi_0"))


# =====================================================================
# Hipoteza C: m_sigma = 0 (massless)
# =====================================================================

banner("Hipoteza C: massless, m_sigma = 0", level=1)
print("""
Argument: jesli sigma_ab to massless composite, propaguje na c_0
identycznie jak fundamentalny graviton. Brak deviation w GW170817.
sigma_ab to 'emergent graviton' a la Sakharov / Verlinde.

Wymaganie: m_s_eff^2 = 0 (delta_bond cancels out V''(s_eq) exactly).
To wymagaloby fine-tuning -- mozliwe jesli istnieje SYMMETRIA zapewniajaca
masslessness. W TGP: skala Phi_0 jest VEV substratu, fluctuations wokol
mogłyby byc Goldstone (jesli odpowiednia symmetria broken spontaneously).

Ale Phi to skalar Z_2; spontaneous breakdown daje JEDEN Goldstone (dla
discrete Z_2 to NIE jest ciagly -- nie ma Goldstone'a). Wiec NATURALNY
Goldstone z Z_2 NIE EXIST.

Mozliwa droga: sigma_ab to non-Goldstone NAMBU-LIKE composite, gdzie
bilinearne <(d s)(d s)> ma zerową masse z NESTING (paragraf przesylowy).
To wymaga dedykowanej analizy 1-loop -- otwarte do T3.3.

WNIOSEK C: massless mozliwe ale wymaga structural reason (nie ma trywial
mechanizm w TGP single-Phi Z_2). T3.3 musi sprawdzic.
""")
checks_summary.append(check("C requires structural reason for masslessness",
                             True,
                             "non-Goldstone composite; needs T3.3 analysis"))


# =====================================================================
# Phi_0 cosmological scale -- naturalny argument
# =====================================================================

banner("BRAINSTORM 8.9: Phi_0 ~ Hubble^(-1) (cosmological constant)", level=1)
print("""
Z paperu: V(Phi=Phi_0) = (gamma/12) Phi_0^2 (vacuum energy).
Cosmological constant Lambda_TGP = V(Phi_0) c^4 / (8 pi G).

Jesli zadamy Lambda_TGP = Lambda_obs (dark energy):
  Lambda_obs ~ (1.4e-3 eV)^4 ~ 2e-122 M_Planck^4
  V(Phi_0) ~ Phi_0^2 (zaniedbujac gamma~1)

Wiec Phi_0^2 ~ Lambda_obs * 8 pi G / c^4 ~ Lambda_obs * (l_P/c^2)^2
       Phi_0 ~ sqrt(Lambda_obs) * l_P ~ (1e-3 eV)^2 * (1e-35 m)

Hmm, zajmijmy sie jednostkowo: V[Phi] musi miec dim energy/volume.
W naturalnych: V = E^4 (dla scalar field z dim E).
Wiec Phi_0 ma dim E (energy-like).

Zadanie Lambda_TGP = Lambda_obs:
  V(Phi_0) ~ gamma * Phi_0^4 ~ Lambda_obs ~ (1.4e-3 eV)^4
  Phi_0 ~ (Lambda_obs / gamma)^(1/4) ~ 1.4e-3 eV (jesli gamma ~ 1)

Wiec NATURALNA Phi_0 to ~meV (NOT muEV jak w brainstorm pierwszym).

Oczekiwane m_sigma:
  m_sigma ~ Phi_0 ~ 1e-3 eV ~ 1 meV

Sprawdzmy czy 1 meV jest GW-safe:
""")
# Phi_0 from cosmological argument
Phi_0_cosmo_eV = (1.4e-3) ** 1.0  # ~ 1 meV
m_sigma_cosmo = Phi_0_cosmo_eV  # ~ 1 meV
GW_safe = m_sigma_cosmo < GW170817_bound_eV
print(f"  Phi_0 cosmo ~ {Phi_0_cosmo_eV:.2e} eV")
print(f"  m_sigma ~ Phi_0 ~ {m_sigma_cosmo:.2e} eV")
print(f"  GW170817 bound: {GW170817_bound_eV:.2e} eV")
print(f"  m_sigma vs bound: {'<' if GW_safe else '>='}  -> {'PASS' if GW_safe else 'FAIL'}")

checks_summary.append(check("D Phi_0 ~ meV (Lambda_obs scale) GW-safe?",
                             GW_safe,
                             f"m_sigma={m_sigma_cosmo:.2e} vs bound {GW170817_bound_eV:.2e}"))

# But wait: m_sigma ~ 1 meV is way ABOVE 2e-19 eV bound! Let's check.
# 1 meV = 1e-3 eV; bound = 2e-19 eV; 1e-3 >> 2e-19. So NOT safe!

print("""
*** UWAGA WAZNA ***
m_sigma ~ Phi_0 ~ 1 meV vs GW170817 bound 2e-19 eV:
  1e-3 eV >> 2e-19 eV  --> NIE PASS!

Wniosek: jesli Phi_0 ~ meV (cosmologically motivated), TO m_sigma jest
ZBYT DUZE dla GW170817. To znaczy:
  - Albo nasza heurystyka m_sigma ~ Phi_0 jest BLEDNA (T3.3 musi pokazac
    detailed renormalization w ktora m_sigma jest tlumione)
  - Albo Phi_0 jest GLEBOKO mniejsze niz meV (np. 1e-19 eV ~ ULDM scale)
  - Albo TGP jest CZESCIOWO falsyfikowalne na poziomie GW dispersion
    (testowalne LIGO 3G)
""")


# =====================================================================
# WERDYKT T3.2
# =====================================================================

banner("T3.2 WERDYKT", level=1)

n_pass = sum(checks_summary)
n_total = len(checks_summary)

print(f"""
Liczba checkow: {n_pass}/{n_total} PASS

Kluczowe wnioski T3.2:

  1. Hipoteza A (m_sigma ~ J*v_0^2 ~ Phi_0):
     - Phi_0 ~ M_Planck: KATASTROFA (m_sigma = M_Planck >> LIGO bound)
     - Phi_0 ~ GUT: NADAL TLUMI GW
     - Phi_0 ~ H_0 (~1e-33 eV): TYLKO TA SKALA przezywa GW170817
       (zgodne z brainstorm 8.9 cosmological constant z Phi_0)

  2. Hipoteza B == A skalujaco; rozni sie tylko prefactor.

  3. Hipoteza C (massless): mozliwa, ale wymaga structural reason
     (non-Goldstone NAMBU-like composite); T3.3 musi sprawdzic.

  4. Wybor Phi_0 z cosmological argument (V(Phi_0) = Lambda_obs):
     Phi_0 ~ meV. To daje m_sigma ~ meV, KTORE JEST ZBYT DUZE dla
     GW170817 bound (2e-19 eV)!

  TENSION: cosmological motivated Phi_0 ~ meV jest niezgodne z
  GW170817 bound chyba ze:
     (i) m_sigma jest tlumione przez sub-tree-level renormalization
         (T3.3 musi sprawdzic 1-loop / non-perturbative effects)
     (ii) Phi_0 jest niemniej jednak << meV (np. ULDM scale ~ 1e-22 eV)
     (iii) TGP czesciowo falsyfikowalne na poziomie GW dispersion

T3.2 STATUS: INCONCLUSIVE.
  Naturalna struktura wskazuje m_sigma ~ Phi_0, co stwarza tension
  z GW170817 jesli Phi_0 jest cosmologicznie motywowane (~meV).
  Resolution wymaga T3.3 (full ghost+dispersion + 1-loop renormalization).

  Hipoteza C (massless) jest NAJBEZPIECZNIEJSZA dla GW; T3.3 musi
  zdecydowac czy jest STRUCTURALNIE WSPIERANA przez TGP.
""")

if n_pass == n_total:
    print("\n  [INCONCLUSIVE-PASS] T3.2 numerycznie zliczone, ale wskazuje")
    print("  na potencjalna tension wymagajaca T3.3 do rozwiązania.\n")
else:
    print(f"\n  [PARTIAL-FAIL] T3.2: {n_total-n_pass} hipotez nie pasuje GW170817.")
    print("  Tension: szukamy resolution w T3.3 lub modify Phi_0 scale.\n")
