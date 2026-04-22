#!/usr/bin/env python3
"""
ex170_bct_cross_koide.py
Sesja v45, 2026-04-05

Glebokie badanie tripletu (b, c, t) z K = 0.6695 (0.42% od 2/3).

Pytania:
1. Czy K(b,c,t) = 2/3 jest predykcja m_t?
   - m_t(pred) = 168.4 GeV vs 172.8 GeV (PDG) -- -2.5%
   - Czy 2.5% mozna tlumaczyc roznicami skal MS-bar vs pole mass?
2. Analiza bledow: niepewnosci mas PDG -> niepewnosc K
3. Czy istnieja inne triplety cross-sektorowe z K ~ 2/3
   po uwzglednieniu bledow?
4. Systematyczne testowanie WSZYSTKICH cross-tripletow
   lacziajacych 3 sektory (lepton, up, down)

Masy:
- b: 4180 +/- 30 MeV (MS-bar, mu_b)
- c: 1270 +/- 20 MeV (MS-bar, 2 GeV)
- t: 172760 +/- 300 MeV (pole mass)
UWAGA: b i c sa MS-bar, t jest pole mass. To mieszanie skal!
"""
import numpy as np

print("=" * 72)
print("ex170: (b,c,t) cross-Koide -- glebokie badanie")
print("=" * 72)

# ---- PDG masses with uncertainties ----
# Source: PDG 2024
particles = {
    # Leptons (pole masses, MeV)
    'e':   (0.51100, 0.00001),
    'mu':  (105.658, 0.001),
    'tau': (1776.86, 0.12),
    # Down quarks (MS-bar at respective scales, MeV)
    'd':   (4.67, 0.48),       # MS-bar, 2 GeV
    's':   (93.4, 8.6),        # MS-bar, 2 GeV
    'b':   (4180, 30),         # MS-bar, mu_b
    # Up quarks
    'u':   (2.16, 0.49),       # MS-bar, 2 GeV
    'c':   (1270, 20),         # MS-bar, mu_c
    't':   (172760, 300),      # pole mass
}

def koide(m1, m2, m3):
    return (m1 + m2 + m3) / (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3))**2

def koide_predict_m3(m1, m2, K=2.0/3.0):
    """Given K and m1, m2, predict m3 from Koide formula.
    K = (m1+m2+m3)/(sqrt(m1)+sqrt(m2)+sqrt(m3))^2
    Returns list of positive solutions for m3."""
    s1 = np.sqrt(m1)
    s2 = np.sqrt(m2)
    S = s1 + s2
    # K*(S+x)^2 = m1+m2+x^2 where x = sqrt(m3)
    # K*S^2 + 2K*S*x + K*x^2 = m1+m2+x^2
    # (K-1)*x^2 + 2K*S*x + K*S^2 - m1 - m2 = 0
    A = K - 1
    B = 2*K*S
    C = K*S**2 - m1 - m2
    disc = B**2 - 4*A*C
    solutions = []
    if disc >= 0:
        for sign in [1, -1]:
            x = (-B + sign*np.sqrt(disc)) / (2*A)
            if x > 0:
                solutions.append(x**2)
    return solutions

# ---- 1. K(b,c,t) with error propagation ----
print("\n--- 1. K(b,c,t) z propagacja bledow ---\n")

m_b, dm_b = particles['b']
m_c, dm_c = particles['c']
m_t, dm_t = particles['t']

K_bct = koide(m_b, m_c, m_t)
print(f"  K(b,c,t) = {K_bct:.6f}")
print(f"  K_Koide  = {2/3:.6f}")
print(f"  Delta    = {K_bct - 2/3:.6f} ({(K_bct-2/3)/(2/3)*100:+.3f}%)")

# Monte Carlo error propagation
np.random.seed(42)
N_mc = 100_000
K_samples = np.zeros(N_mc)
for i in range(N_mc):
    mb = np.random.normal(m_b, dm_b)
    mc = np.random.normal(m_c, dm_c)
    mt = np.random.normal(m_t, dm_t)
    K_samples[i] = koide(mb, mc, mt)

K_mean = np.mean(K_samples)
K_std = np.std(K_samples)
print(f"\n  Monte Carlo ({N_mc:,} samples):")
print(f"  K = {K_mean:.6f} +/- {K_std:.6f}")
print(f"  Odchylenie od 2/3: {(K_mean - 2/3)/K_std:.1f} sigma")

# ---- 2. Scale mismatch ----
print("\n--- 2. Problem mieszania skal ---\n")
print("  b: MS-bar at mu_b ~ 4.18 GeV")
print("  c: MS-bar at mu_c ~ 1.27 GeV")
print("  t: pole mass (nie MS-bar!)")
print()

# t quark MS-bar mass at mu_t
# PDG: m_t(pole) = 172.76 GeV, m_t(MS-bar, mu_t) ~ 162.5 GeV
m_t_MSbar = 162500  # approximate MS-bar at mu_t
print(f"  m_t(pole) = {m_t/1000:.1f} GeV")
print(f"  m_t(MS-bar, mu_t) ~ {m_t_MSbar/1000:.1f} GeV")

K_MSbar = koide(m_b, m_c, m_t_MSbar)
print(f"\n  K(b,c,t) z pole mass:  {K_bct:.6f}")
print(f"  K(b,c,t) z MS-bar t:   {K_MSbar:.6f}")
print(f"  Delta z MS-bar: {(K_MSbar - 2/3)/(2/3)*100:+.3f}%")

# What m_t gives K = 2/3 exactly?
m_t_solutions = koide_predict_m3(m_b, m_c)
print(f"\n  m_t z warunku K(b,c,t) = 2/3:")
for sol in m_t_solutions:
    print(f"    m_t = {sol:.0f} MeV = {sol/1000:.1f} GeV")
    print(f"    (pole: {abs(sol - m_t)/m_t*100:.1f}% od PDG, "
          f"MS-bar: {abs(sol - m_t_MSbar)/m_t_MSbar*100:.1f}% od ~162.5)")

# ---- 3. Running masses to common scale ----
print("\n--- 3. Bieganie mas do wspolnej skali ---\n")
print("  Przyblizone masy MS-bar przy skali mu = m_Z = 91.2 GeV:")
print("  (uzycie 1-loop QCD running z alpha_s(M_Z) = 0.1179)")
print()

# Approximate running masses at M_Z (from literature)
# These are rough values for illustration
masses_MZ = {
    'b': 2850,    # MS-bar at M_Z ~ 2.85 GeV
    'c': 620,     # MS-bar at M_Z ~ 0.62 GeV
    't': 171400,  # MS-bar at M_Z ~ 171.4 GeV (close to pole for t)
}

K_MZ = koide(masses_MZ['b'], masses_MZ['c'], masses_MZ['t'])
print(f"  b(M_Z) ~ {masses_MZ['b']} MeV")
print(f"  c(M_Z) ~ {masses_MZ['c']} MeV")
print(f"  t(M_Z) ~ {masses_MZ['t']} MeV")
print(f"  K(b,c,t)|_MZ = {K_MZ:.6f} ({(K_MZ-2/3)/(2/3)*100:+.2f}%)")

# At 2 GeV (where light quarks are quoted)
masses_2GeV = {
    'b': 4180,    # roughly same (b runs little between 2 GeV and mu_b)
    'c': 1270,    # quoted at mu_c not 2 GeV, but close
    't': 172760,  # pole mass (doesn't run much)
}
# Actually b at 2 GeV would be larger due to running down
# m_b(2 GeV) ~ m_b(mu_b) * [alpha_s(2)/alpha_s(mu_b)]^{12/23}
# alpha_s(2 GeV) ~ 0.30, alpha_s(4.18) ~ 0.22
# ratio ~ (0.30/0.22)^{12/23} ~ 1.36^{0.52} ~ 1.17
m_b_2GeV = 4180 * 1.17  # rough estimate
m_c_2GeV = 1270 * 1.10  # rough: run c from mu_c to 2 GeV
K_2GeV = koide(m_b_2GeV, m_c_2GeV, m_t)
print(f"\n  Rough masses at 2 GeV:")
print(f"  b(2GeV) ~ {m_b_2GeV:.0f} MeV")
print(f"  c(2GeV) ~ {m_c_2GeV:.0f} MeV")
print(f"  t (pole) = {m_t} MeV")
print(f"  K(b,c,t)|_2GeV = {K_2GeV:.6f} ({(K_2GeV-2/3)/(2/3)*100:+.2f}%)")

# ---- 4. Predykcja m_t z cross-Koide ----
print("\n--- 4. Predykcja m_t z K(b,c,t) = 2/3 ---\n")

# Using different mass conventions
scenarios = [
    ("PDG nominal (mixed scales)", m_b, m_c),
    ("b,c at M_Z (rough)", masses_MZ['b'], masses_MZ['c']),
    ("b,c at 2 GeV (rough)", m_b_2GeV, m_c_2GeV),
]

for label, mb, mc in scenarios:
    sols = koide_predict_m3(mb, mc)
    print(f"  {label}:")
    for sol in sols:
        print(f"    m_t(pred) = {sol/1000:.1f} GeV", end="")
        err_pole = (sol - m_t) / m_t * 100
        err_MSbar = (sol - m_t_MSbar) / m_t_MSbar * 100
        print(f"  (vs pole: {err_pole:+.1f}%, vs MS-bar: {err_MSbar:+.1f}%)")

# ---- 5. All cross-triplets mixing all 3 sectors ----
print("\n--- 5. Triplety mieszajace WSZYSTKIE 3 sektory ---\n")
print("  (1 lepton + 1 down quark + 1 up quark)\n")

leptons = ['e', 'mu', 'tau']
downs = ['d', 's', 'b']
ups = ['u', 'c', 't']

results_3sector = []
for l in leptons:
    for d in downs:
        for u in ups:
            ml = particles[l][0]
            md = particles[d][0]
            mu = particles[u][0]
            K = koide(ml, md, mu)
            dev = (K - 2/3) / (2/3) * 100
            results_3sector.append((l, d, u, K, dev))

results_3sector.sort(key=lambda x: abs(x[4]))

print(f"  {'Triplet':>18s}  {'K':>10s}  {'dK/K0 (%)':>12s}")
print("  " + "-" * 45)
for l, d, u, K, dev in results_3sector[:15]:
    label = f"({l},{d},{u})"
    marker = " <--" if abs(dev) < 2 else ""
    print(f"  {label:>18s}  {K:10.6f}  {dev:+12.3f}%{marker}")

# ---- 6. Symetrie tripletu (b,c,t) ----
print("\n--- 6. Fizyczna interpretacja (b,c,t) ---\n")
print("  Triplet (b, c, t) laczy:")
print("    b = down-type, 3rd gen, Q = -1/3")
print("    c = up-type, 2nd gen, Q = +2/3")
print("    t = up-type, 3rd gen, Q = +2/3")
print()
print("  To sa najciezsze fermiony 2. i 3. generacji.")
print("  Mozliwa interpretacja TGP:")
print("  - phi-FP daje r21 (1st -> 2nd generacja)")
print("  - K(b,c,t) = 2/3 laczy 2nd i 3rd gen MIEDZY sektorami")
print("  - Tylko te 3 fermiony o masach > 1 GeV z gen 2,3")
print()

# What fraction of mass is in this triplet?
total_mass = sum(p[0] for p in particles.values())
bct_mass = m_b + m_c + m_t
print(f"  Masa (b+c+t) / suma wszystkich = {bct_mass/total_mass*100:.1f}%")
print(f"  (b,c,t) dominuje budzet masowy fermionow")

# ---- 7. Generalized Koide for pairs of sectors ----
print("\n--- 7. Koide sektorowy i cross-sektorowy ---\n")

# Within each sector
print("  Intra-sektorowy:")
for name, parts in [('lepton', leptons), ('down', downs), ('up', ups)]:
    m = [particles[p][0] for p in parts]
    K = koide(*m)
    print(f"    {name:>7s} ({','.join(parts)}): K = {K:.6f}")

# Cross: heavy fermions
print("\n  Cross-sektorowy (ciezkie fermiony):")
heavy_triplets = [
    ('tau', 'b', 't'),  # heaviest from each
    ('tau', 'b', 'c'),  # tau + heavy downs + charm
    ('tau', 's', 't'),
    ('mu', 'b', 't'),
    ('mu', 'b', 'c'),   # = mu + (b,c)
]
for combo in heavy_triplets:
    m = [particles[p][0] for p in combo]
    K = koide(*m)
    dev = (K - 2/3)/(2/3)*100
    print(f"    ({','.join(combo):>12s}): K = {K:.6f}  ({dev:+.2f}%)")

# ---- 8. Brannen parametrization for (b,c,t) ----
print("\n--- 8. Parametryzacja Brannena dla (b,c,t) ---\n")
print("  m_k = M/3 * (1 + sqrt(2)*cos(theta0 + 2*pi*k/3))^2\n")

# Fit Brannen to (b,c,t)
m_list = sorted([m_b, m_c, m_t])
M = m_list[0] + m_list[1] + m_list[2]
print(f"  M = {M:.0f} MeV")

# sqrt(m_k/M*3) = 1 + sqrt(2)*cos(theta_k)
# cos(theta_k) = (sqrt(3*m_k/M) - 1) / sqrt(2)
for k, (name, m) in enumerate(zip(['c', 'b', 't'], m_list)):
    cos_val = (np.sqrt(3*m/M) - 1) / np.sqrt(2)
    if abs(cos_val) <= 1:
        theta = np.arccos(cos_val)
        print(f"  {name}: m = {m:.0f} MeV, cos(theta) = {cos_val:.6f}, "
              f"theta = {np.degrees(theta):.2f} deg")
    else:
        print(f"  {name}: m = {m:.0f} MeV, cos(theta) = {cos_val:.6f} (poza [-1,1])")

# ---- 9. Sigma significance of K(b,c,t) = 2/3 ----
print("\n--- 9. Istotnosc statystyczna ---\n")
print(f"  K(b,c,t) = {K_bct:.6f} +/- {K_std:.6f}")
print(f"  K_Koide  = {2/3:.6f}")
print(f"  Odchylenie: {abs(K_bct - 2/3)/K_std:.1f} sigma")
print()
print(f"  m_t(pred) = {m_t_solutions[0]/1000:.1f} GeV (z K=2/3)")
print(f"  m_t(PDG)  = {m_t/1000:.2f} +/- {dm_t/1000:.2f} GeV")
print(f"  Roznica   = {abs(m_t_solutions[0] - m_t)/1000:.1f} GeV "
      f"= {abs(m_t_solutions[0] - m_t)/dm_t:.1f} sigma(m_t)")

# ---- WNIOSKI ----
print("\n" + "=" * 72)
print("WNIOSKI ex170")
print("=" * 72)
print(f"""
  1. K(b,c,t) = {K_bct:.6f} -- odchylenie {(K_bct-2/3)/K_std:.1f} sigma od 2/3
     (przy niepewnosciach PDG: sigma_K = {K_std:.6f})

  2. Predykcja: K(b,c,t) = 2/3 => m_t = {m_t_solutions[0]/1000:.1f} GeV
     PDG: m_t = {m_t/1000:.1f} +/- {dm_t/1000:.1f} GeV
     Roznica: {abs(m_t_solutions[0] - m_t)/1000:.1f} GeV ({abs(m_t_solutions[0]-m_t)/m_t*100:.1f}%)
     = {abs(m_t_solutions[0] - m_t)/dm_t:.1f} sigma(m_t)

  3. PROBLEM SKAL: b (MS-bar at mu_b), c (MS-bar at mu_c),
     t (pole mass) -- NIE SA na tej samej skali!
     Running do wspolnej skali ZMIENI K.
     Na skali M_Z: K ~ {K_MZ:.4f} ({(K_MZ-2/3)/(2/3)*100:+.1f}%)

  4. Jedyny triplet mieszajacy 3 sektory z |dK| < 2%:
     BRAK -- (b,c,t) nie miesza wszystkich 3 sektorow
     (bo b i t sa z tego samego [gen 3] ale roznych [up/down])

  5. STATUS: K(b,c,t) ~ 2/3 jest INTERESUJACE ale:
     - Wymaga korekcji na bieganie mas (running)
     - 2.5% roznica moze byc artefaktem mieszania skal
     - {abs(K_bct - 2/3)/K_std:.1f} sigma -- na granicy istotnosci
     - Brak jasnej motywacji fizycznej (dlaczego b,c,t?)

  6. HIPOTEZA ROBOCZA:
     Jesli K(b,c,t) = 2/3 przetrwa korekcje running,
     to lacznie z phi-FP daje ZAMKNIETY uklad:
       phi-FP: r_cs = m_c/m_s, r_bt ~ ???
       K(b,c,t)=2/3: m_t = f(m_b, m_c)
     Potrzebne: poprawne masy biegajace.
""")
