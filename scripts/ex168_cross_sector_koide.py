#!/usr/bin/env python3
"""
ex168_cross_sector_koide.py
Sesja v45, 2026-04-05

Systematyczny skan cross-sektorowych tripletow Koide.
Masy PDG (MS-bar, 2 GeV) + masy biegunowe leptonów.

Pytanie: czy istnieją triplety międzysektorowe z K ~ 2/3?
"""
import numpy as np
from itertools import combinations

# ── Masy PDG ──
masses = {
    # Leptony (pole masses, MeV)
    'e':   0.51100,
    'mu':  105.658,
    'tau': 1776.86,
    # Down-type quarks (MS-bar, 2 GeV, MeV)
    'd':   4.67,
    's':   93.4,
    'b':   4180.0,
    # Up-type quarks (MS-bar, 2 GeV, MeV)
    'u':   2.16,
    'c':   1270.0,
    't':   172760.0,
}

K_target = 2.0 / 3.0

def koide(m1, m2, m3):
    """K = (m1+m2+m3) / (√m1 + √m2 + √m3)²"""
    s = m1 + m2 + m3
    sq = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s / sq**2

print("=" * 72)
print("ex168: Systematyczny skan cross-sektorowych tripletow Koide")
print("=" * 72)

# ── 1. Standard Koide w kazdym sektorze ──
print("\n--- 1. Standardowe triplety (intra-sector) ---\n")
sectors = [
    ('e', 'mu', 'tau', 'leptony'),
    ('d', 's', 'b', 'down'),
    ('u', 'c', 't', 'up'),
]
for a, b, c, name in sectors:
    K = koide(masses[a], masses[b], masses[c])
    dev = (K - K_target) / K_target * 100
    ok = 'OK' if abs(dev) < 1 else 'NO'
    print(f"  ({a}, {b}, {c}) [{name:>7s}]: K = {K:.6f}  "
          f"dK/K0 = {dev:+.2f}%  {ok}")

# ── 2. WSZYSTKIE triplety z 9 cząstek ──
print("\n--- 2. Wszystkie triplety (9 choose 3 = 84) ---\n")

particles = list(masses.keys())
results = []

for combo in combinations(particles, 3):
    m1, m2, m3 = masses[combo[0]], masses[combo[1]], masses[combo[2]]
    K = koide(m1, m2, m3)
    dev = (K - K_target) / K_target * 100
    results.append((combo, K, dev))

# Sort by |deviation|
results.sort(key=lambda x: abs(x[2]))

print(f"  {'Triplet':>18s}  {'K':>10s}  {'dK/K0 (%)':>12s}  {'Status':>6s}")
print("  " + "-" * 55)
for combo, K, dev in results[:20]:
    label = f"({', '.join(combo)})"
    status = 'OK' if abs(dev) < 1 else ('~' if abs(dev) < 5 else 'NO')
    print(f"  {label:>18s}  {K:10.6f}  {dev:+12.3f}%  {status:>6s}")

print(f"\n  ... (84 tripletow lacznie, pokazano 20 najblizszych K=2/3)")

# ── 3. Wyrozone triplety cross-sektor ──
print("\n--- 3. Wyrozone triplety cross-sektorowe ---\n")
cross_special = [
    ('u', 's', 'tau'),   # H7 z dodatku X
    ('tau', 'b', 't'),   # H7 z dodatku X
    ('e', 'd', 'u'),     # lightest from each sector
    ('mu', 's', 'c'),    # middle from each
    ('tau', 'b', 't'),   # heaviest from each
    ('e', 'c', 'b'),     # mixed
    ('mu', 'c', 'tau'),  # mixed
    ('e', 's', 'b'),     # e + down quarks 2,3
]

seen = set()
for combo in cross_special:
    key = tuple(sorted(combo))
    if key in seen:
        continue
    seen.add(key)
    m1, m2, m3 = masses[combo[0]], masses[combo[1]], masses[combo[2]]
    K = koide(m1, m2, m3)
    dev = (K - K_target) / K_target * 100
    status = 'OK' if abs(dev) < 1 else ('~' if abs(dev) < 5 else 'NO')
    print(f"  ({', '.join(combo):>8s}): K = {K:.6f}  "
          f"dK/K0 = {dev:+.3f}%  {status}")

# ── 4. Szukanie tripletow z K ~ = 2/3 ± 1% ──
print("\n--- 4. Triplety z |dK/K0| < 1% ---\n")
close = [(c, K, d) for c, K, d in results if abs(d) < 1.0]
if close:
    for combo, K, dev in close:
        label = f"({', '.join(combo)})"
        # Classify: intra or cross?
        sector_map = {'e': 'L', 'mu': 'L', 'tau': 'L',
                      'd': 'D', 's': 'D', 'b': 'D',
                      'u': 'U', 'c': 'U', 't': 'U'}
        sectors_in = set(sector_map[p] for p in combo)
        cross = "CROSS" if len(sectors_in) > 1 else "intra"
        print(f"  {label:>18s}  K = {K:.6f}  dK = {dev:+.4f}%  [{cross}]")
else:
    print("  Brak tripletow z |dK/K0| < 1%.")

# ── 5. Triplety z |dK/K0| < 5% ──
print("\n--- 5. Triplety z |dK/K0| < 5% ---\n")
moderate = [(c, K, d) for c, K, d in results if abs(d) < 5.0]
if moderate:
    for combo, K, dev in moderate:
        label = f"({', '.join(combo)})"
        sector_map = {'e': 'L', 'mu': 'L', 'tau': 'L',
                      'd': 'D', 's': 'D', 'b': 'D',
                      'u': 'U', 'c': 'U', 't': 'U'}
        sectors_in = set(sector_map[p] for p in combo)
        cross = "CROSS" if len(sectors_in) > 1 else "intra"
        print(f"  {label:>18s}  K = {K:.6f}  dK = {dev:+.3f}%  [{cross}]")

# ── 6. Statystyka: ile tripletow na przedział ──
print("\n--- 6. Rozklad K ---\n")
Ks = [K for _, K, _ in results]
bins = [0.33, 0.50, 0.60, 0.65, 0.67, 0.70, 0.75, 0.80, 0.90, 1.00]
for i in range(len(bins) - 1):
    count = sum(1 for K in Ks if bins[i] <= K < bins[i + 1])
    bar = '#' * count
    print(f"  [{bins[i]:.2f}, {bins[i+1]:.2f}): {count:2d}  {bar}")

# ── 7. Test prawdopodobieństwa: losowe triplety ──
print("\n--- 7. Prawdopodobienstwo losowego K ~ 2/3 ---\n")
np.random.seed(42)
N_random = 100_000
count_close = 0
for _ in range(N_random):
    # Log-uniform masses in range [0.5 MeV, 200 GeV]
    log_m = np.random.uniform(np.log(0.5), np.log(200000), 3)
    m = np.exp(log_m)
    K = koide(m[0], m[1], m[2])
    if abs((K - K_target) / K_target) < 0.01:  # 1%
        count_close += 1
prob = count_close / N_random * 100
print(f"  P(|dK| < 1%) z losowych mas log-uniform: {prob:.1f}%")
print(f"  (z {N_random:,} losowan)")

# ── Wnioski ──
print("\n" + "=" * 72)
print("WNIOSKI ex168")
print("=" * 72)

n_close_1 = len(close)
n_moderate = len(moderate)

print(f"""
  1. Triplety z |dK/K0| < 1%: {n_close_1} (z 84)
  2. Triplety z |dK/K0| < 5%: {n_moderate} (z 84)
  3. Losowe P(|dK| < 1%) = {prob:.1f}%

  Interpretacja:
  - Jesli {n_close_1} cross-sektorowych tripletow ma K ~ 2/3,
    a losowa szansa to {prob:.1f}%, to:
    Oczekiwana liczba = 84 × {prob/100:.3f} = {84*prob/100:.1f}
  - Porownanie z obserwowana liczba {n_close_1} pozwala ocenic,
    czy cross-sektorowe K ~ 2/3 jest statystycznie istotne.
""")
