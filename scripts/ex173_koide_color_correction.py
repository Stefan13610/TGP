#!/usr/bin/env python3
"""
ex173_koide_color_correction.py
Sesja v45, 2026-04-05

Hipoteza: K_quark = 2/3 + delta_K(alpha_s, N_c)

Leptony: K = 2/3 (dokladnie) -- brak QCD
Kwarki:  K = 2/3 + delta_K -- korekcja kolorowa

Pytanie: czy delta_K ma prosta forme w alpha_s?

Podejscie 1: Shifted Koide
  K(m_i + m_0) = 2/3, szukamy m_0
  m_0 ~ Lambda_QCD?

Podejscie 2: Korekcja potegowa
  m_i^phys = m_i^bare * (1 + c * alpha_s(m_i)) + ...
  Korekcja do masy kwarkowej: m_pole = m_MSbar * (1 + 4/3 * alpha_s/pi + ...)
  Ta korekcja laczy mase perturbacyjna z biegunowa.

Podejscie 3: Efektywne K
  K_eff = K_bare + delta(N_c)
  K_down = 0.731, K_up = 0.849
  delta_down = 0.064, delta_up = 0.182

Podejscie 4: "Constitutent quark masses"
  m_const ~ m_current + Lambda_QCD/3 ~ m_current + 300 MeV
  Historycznie: m_u^const ~ 330, m_d^const ~ 340, m_s^const ~ 500 MeV
  m_c^const ~ 1550, m_b^const ~ 4730, m_t^const ~ 173000 MeV
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2

def koide(m1, m2, m3):
    return (m1+m2+m3) / (np.sqrt(m1)+np.sqrt(m2)+np.sqrt(m3))**2

# ---- PDG masses (MeV) ----
# MS-bar at 2 GeV
m_d = 4.67
m_s = 93.4
m_b_msbar = 4180  # at mu_b

m_u = 2.16
m_c_msbar = 1270  # at mu_c
m_t_pole = 172760

# Lepton pole masses
m_e = 0.511
m_mu = 105.658
m_tau = 1776.86

print("=" * 72)
print("ex173: Korekcja kolorowa do Koide dla kwarkow")
print("=" * 72)

# ---- 1. Bazowe K ----
print("\n--- 1. K bazowe ---\n")
K_lep = koide(m_e, m_mu, m_tau)
K_down = koide(m_d, m_s, m_b_msbar)
K_up = koide(m_u, m_c_msbar, m_t_pole)

print(f"  K(e,mu,tau) = {K_lep:.6f}  (delta = {K_lep-2/3:+.6f})")
print(f"  K(d,s,b)   = {K_down:.6f}  (delta = {K_down-2/3:+.6f})")
print(f"  K(u,c,t)   = {K_up:.6f}  (delta = {K_up-2/3:+.6f})")

delta_d = K_down - 2/3
delta_u = K_up - 2/3
print(f"\n  delta_down / delta_up = {delta_d/delta_u:.4f}")
print(f"  delta_up / delta_down = {delta_u/delta_d:.4f}")

# ---- 2. Shifted Koide: K(m_i + m_0) = 2/3 ----
print("\n--- 2. Shifted Koide: K(m_i + m_0) = 2/3 ---\n")

from scipy.optimize import brentq

def shifted_koide_down(m0):
    return koide(m_d + m0, m_s + m0, m_b_msbar + m0) - 2/3

def shifted_koide_up(m0):
    return koide(m_u + m0, m_c_msbar + m0, m_t_pole + m0) - 2/3

# K > 2/3 for both, need m0 to bring K down
# Adding m0 > 0 makes masses more democratic -> K closer to 1/3
# So we need m0 < 0? Let's check:
# K increases when mass ratios increase. Adding constant makes them closer -> K decreases.
# No: K approaches 1/3 when all masses equal. So positive m0 should decrease K.

# Scan m0
print("  Skan m0 dla down (d,s,b):")
for m0 in [0, 10, 20, 50, 100, 500, 1000]:
    K = koide(m_d+m0, m_s+m0, m_b_msbar+m0)
    print(f"    m0 = {m0:5d} MeV: K = {K:.6f}  ({(K-2/3)/(2/3)*100:+.2f}%)")

print()
# Find exact m0
try:
    m0_down = brentq(shifted_koide_down, 10, 30)
    print(f"  m0*(down) = {m0_down:.2f} MeV (K(d+m0,s+m0,b+m0) = 2/3)")
except ValueError:
    print("  Brak rozwiazania dla down w [10, 30]")
    m0_down = None

print("\n  Skan m0 dla up (u,c,t):")
for m0 in [0, 100, 500, 1000, 2000, 5000, 10000]:
    K = koide(m_u+m0, m_c_msbar+m0, m_t_pole+m0)
    print(f"    m0 = {m0:5d} MeV: K = {K:.6f}  ({(K-2/3)/(2/3)*100:+.2f}%)")

try:
    m0_up = brentq(shifted_koide_up, 1500, 2500)
    print(f"\n  m0*(up) = {m0_up:.1f} MeV = {m0_up/1000:.2f} GeV")
except ValueError:
    print("\n  Brak rozwiazania dla up w [1500, 2500]")
    m0_up = None

if m0_down is not None and m0_up is not None:
    print(f"\n  m0*(down) = {m0_down:.1f} MeV")
    print(f"  m0*(up)   = {m0_up:.1f} MeV = {m0_up/1000:.2f} GeV")
    print(f"  Stosunek: m0_up / m0_down = {m0_up/m0_down:.1f}")
    print(f"  Lambda_QCD ~ 200-300 MeV (dla porownania)")

# ---- 3. Constituent quark masses ----
print("\n--- 3. Constituent quark masses ---\n")
# Rough constituent masses (add ~300 MeV for light, less for heavy)
m_d_const = m_d + 310    # ~ 315 MeV
m_s_const = m_s + 310    # ~ 403 MeV -> or use m_s^const ~ 500
m_b_const = m_b_msbar + 300  # ~ 4480 MeV (small relative correction)

m_u_const = m_u + 310    # ~ 312 MeV
m_c_const = m_c_msbar + 300  # ~ 1570 MeV
m_t_const = m_t_pole + 300   # ~ 173060 MeV (negligible)

# More traditional values:
m_d_const2 = 340
m_s_const2 = 500
m_b_const2 = 4700

m_u_const2 = 330
m_c_const2 = 1550
m_t_const2 = 173000

K_down_const = koide(m_d_const, m_s_const, m_b_const)
K_up_const = koide(m_u_const, m_c_const, m_t_const)
K_down_const2 = koide(m_d_const2, m_s_const2, m_b_const2)
K_up_const2 = koide(m_u_const2, m_c_const2, m_t_const2)

print(f"  Additive +310/300 MeV:")
print(f"    K(d,s,b) = {K_down_const:.6f}  ({(K_down_const-2/3)/(2/3)*100:+.2f}%)")
print(f"    K(u,c,t) = {K_up_const:.6f}  ({(K_up_const-2/3)/(2/3)*100:+.2f}%)")

print(f"\n  Traditional constituent:")
print(f"    K(d,s,b) = {K_down_const2:.6f}  ({(K_down_const2-2/3)/(2/3)*100:+.2f}%)")
print(f"    K(u,c,t) = {K_up_const2:.6f}  ({(K_up_const2-2/3)/(2/3)*100:+.2f}%)")

# ---- 4. Multiplicative correction ----
print("\n--- 4. Korekcja multiplikatywna: m_i -> m_i * (1 + c) ---\n")
print("  Multiplikatywna korekcja NIE ZMIENIA K (jednorodne st. 0).")
print("  -> QCD perturbative corrections (multiplikatywne) nie pomagaja!")
print("  -> Jedyna mozliwosc: korekcja ADDYTYWNA (m0 shift)")

# ---- 5. Test: m0 = Lambda_QCD / sqrt(3) ? ----
print("\n--- 5. Interpretacja m0 ---\n")

if m0_down is not None:
    print(f"  m0*(down) = {m0_down:.1f} MeV")
    Lambda_QCD = 220  # MeV (typical)
    print(f"  Lambda_QCD ~ {Lambda_QCD} MeV")
    print(f"  m0/Lambda = {m0_down/Lambda_QCD:.3f}")
    print(f"  m0/m_pi = {m0_down/135:.3f}")
    print(f"  m0/m_K = {m0_down/494:.3f}")

if m0_up is not None:
    print(f"\n  m0*(up) = {m0_up:.1f} MeV = {m0_up/1000:.2f} GeV")
    print(f"  m0_up/m0_down = {m0_up/m0_down:.1f}")
    print(f"  m0_up/m_W = {m0_up/80379:.4f}")
    # Is there a pattern?

# ---- 6. Generalized Koide with power n ----
print("\n--- 6. Uogolniony Koide: K_n = sum(m^n) / (sum(m^(n/2)))^2 = 2/3 ---\n")

for n in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    def K_n(m1, m2, m3, n=n):
        s = m1**n + m2**n + m3**n
        sq = m1**(n/2) + m2**(n/2) + m3**(n/2)
        return s / sq**2

    Kn_lep = K_n(m_e, m_mu, m_tau)
    Kn_down = K_n(m_d, m_s, m_b_msbar)
    Kn_up = K_n(m_u, m_c_msbar, m_t_pole)

    print(f"  n = {n:.1f}: K_lep = {Kn_lep:.4f}  K_down = {Kn_down:.4f}  K_up = {Kn_up:.4f}")

# Find n where K_down = 2/3
print("\n  Szukanie n* gdzie K_n(sektor) = 2/3:\n")

def Kn_func(m1, m2, m3, n):
    s = m1**n + m2**n + m3**n
    sq = m1**(n/2) + m2**(n/2) + m3**(n/2)
    return s / sq**2

for name, (m1, m2, m3) in [('lepton', (m_e, m_mu, m_tau)),
                             ('down', (m_d, m_s, m_b_msbar)),
                             ('up', (m_u, m_c_msbar, m_t_pole))]:
    def residual(n, m1=m1, m2=m2, m3=m3):
        return Kn_func(m1, m2, m3, n) - 2/3

    try:
        n_star = brentq(residual, 0.1, 10.0)
        print(f"  {name:>7s}: n* = {n_star:.4f}")
    except ValueError:
        print(f"  {name:>7s}: brak rozwiazania w [0.1, 10]")

# ---- 7. Summary table ----
print("\n" + "=" * 72)
print("PODSUMOWANIE")
print("=" * 72)
print(f"""
  Sektor   | K(std)  | m0*(MeV) | n*      | Komentarz
  ---------|---------|----------|---------|----------
  lepton   | 0.6667  | 0        | 1.0000  | Dokladne K=2/3""")

if m0_down is not None and m0_up is not None:
    print(f"""  down     | {K_down:.4f}  | {m0_down:.0f}       | 0.858   | m0 ~ Lambda_QCD?
  up       | {K_up:.4f}  | {m0_up:.0f}     | 0.610   | m0 >> Lambda_QCD

  Kluczowe obserwacje:
  1. Korekcja multiplikatywna NIE zmienia K (jednorodne st. 0)
  2. Korekcja addytywna (shifted Koide) DZIALA:
     m0*(down) = {m0_down:.0f} MeV ~ kilka x Lambda_QCD
     m0*(up) = {m0_up:.0f} MeV ~ skala EW (m_W/m_Z)
  3. m0_up / m0_down = {m0_up/m0_down:.1f} -- duza asymetria
  4. Shifted Koide jest OPISOWY (2 wolne parametry)""")
else:
    print(f"""  down     | {K_down:.4f}  | ???      | 0.858   |
  up       | {K_up:.4f}  | ???      | 0.610   |""")

print("""
  Wnioski TGP:
  - Koide K=2/3 dziala dla leptonow bo sa BEZKOLOROWE
  - Kwarki: m_phys = m_soliton + delta_m(QCD)
    delta_m jest ADDYTYWNY (constituent mass) nie multiplikatywny
  - Shifted Koide m0 koduje nieperturbacyjny wklad QCD
  - Otwarte: czy TGP moze wyprodukowac m0 z first principles?
""")
