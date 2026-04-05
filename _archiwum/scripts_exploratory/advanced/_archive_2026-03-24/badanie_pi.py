"""
TGP — Badanie zwiazku Phi0_generacji z liczba pi
=================================================

Obserwacja: dla alpha=5.9, a_gam=0.05:
  Phi0_2 ~ 6.227  vs  2*pi = 6.283  (roznica ~0.9%)
  Phi0_3 ~ 26.53  vs  ???

Pytania:
  1. Czy Phi0_2 -> 2*pi w jakims limicie (a_gam->0? alpha->alfa_krit?)?
  2. Czy Phi0_3 ma rownie prosta postac?
  3. Gdzie w TGP pojawia sie pi algebraicznie?
  4. Czy Koide formula (sqrt-relacja mas) wynika z TGP?

Nota: Lambda = 4*pi ZAWSZE (z geometrii sferycznej).
      Phi0_2 ~ Lambda/2 = 2*pi  byloby bardzo eleganckie!
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

GAMMA = 1.0; PHI0_BASE = 1.0
LAM = 5.514e-7  # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)

# ======================================================================
# ENERGIA i SAMOSPOJNOSC
# ======================================================================
def energy(K, Phi0, alpha, a_gam):
    msp = np.sqrt(max(GAMMA * PHI0_BASE / Phi0, 1e-9))
    r_max = max(60.0 / msp, 20.0)
    r = np.linspace(a_gam, r_max, 3000)
    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1.0 + alpha / (Phi0 * phi)) * r**2, r)
    psi = phi / Phi0
    Ep  = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4 + LAM/6*(psi-1.0)**6) * r**2, r)
    return Ek + Ep

def g(M, Phi0, alpha, a_gam):
    if M < 1e-12: return 0.0
    return energy(M / (4*np.pi*Phi0), Phi0, alpha, a_gam) - M

def find_root(Phi0, alpha, a_gam, M_min=1e-5, M_max=50000, N=400):
    M_arr = np.concatenate([
        np.linspace(M_min, 2,     int(N*0.25)),
        np.linspace(2,     100,   int(N*0.25)),
        np.linspace(100,   2000,  int(N*0.25)),
        np.linspace(2000,  M_max, int(N*0.25))])
    M_arr = np.unique(M_arr)
    g_arr = np.array([g(M, Phi0, alpha, a_gam) for M in M_arr])
    for i in range(len(g_arr)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                return brentq(lambda M: g(M, Phi0, alpha, a_gam),
                              M_arr[i], M_arr[i+1], xtol=1e-8)
            except Exception:
                pass
    return None

def find_Phi0_for_ratio(target_ratio, M1, alpha, a_gam,
                         Phi0_lo=1.001, Phi0_hi=200.0):
    """Brentq po Phi0: szuka Phi0 takie ze M(Phi0)/M1 = target_ratio."""
    def h(p):
        Mn = find_root(p, alpha, a_gam)
        if Mn is None: return target_ratio  # brak korzenia
        return Mn/M1 - target_ratio
    try:
        return brentq(h, Phi0_lo, Phi0_hi, xtol=1e-5, maxiter=40)
    except Exception:
        return None


# ======================================================================
# 1. GDZIE PI POJAWIA SIE W TGP ALGEBRAICZNIE
# ======================================================================
print("="*65)
print("MIEJSCA WYSTAPIENIA PI W TGP (algebraicznie):")
print("="*65)
print()
print("  Lambda = 4*pi         (stala samospojnosci, z calki sferycznej)")
print("  K = M / (4*pi*Phi0)   (sprzezenie do pola)")
print("  Phi0_2 ~ Lambda/2     => Phi0_2 ~ 2*pi  [hipoteza]")
print()
print(f"  Lambda = 4*pi = {4*np.pi:.6f}")
print(f"  Lambda/2 = 2*pi = {2*np.pi:.6f}")
print()

# ======================================================================
# 2. SKAN: jak Phi0_2 zalazy od (alpha, a_gam)?
# ======================================================================
print("="*65)
print("SKAN: Phi0_2 (dla M2/M1=207) vs alpha i a_gam")
print(f"  2*pi = {2*np.pi:.5f},   Lambda = {4*np.pi:.5f}")
print("="*65)
print(f"{'alpha':>7} {'a_gam':>7} {'Phi0_2':>10} {'Phi0_2/2pi':>12} {'diff %':>8} {'Phi0_3':>10} {'Phi0_3/Phi0_2':>14}")
print("-"*75)

results_pi = []
for alpha in [4.0, 5.0, 5.9, 6.0, 7.0, 8.0, 10.0]:
    for a_gam in [0.01, 0.02, 0.03, 0.05, 0.08, 0.10]:
        M1 = find_root(1.0, alpha, a_gam)
        if M1 is None:
            continue
        P2 = find_Phi0_for_ratio(207.0,  M1, alpha, a_gam)
        P3 = find_Phi0_for_ratio(3477.0, M1, alpha, a_gam)
        if P2 is None:
            continue
        ratio_2pi = P2 / (2*np.pi)
        diff_pct  = 100*(P2 - 2*np.pi)/(2*np.pi)
        ratio_P3P2 = P3/P2 if P3 else None
        s_P3 = f"{P3:.4f}" if P3 else "-"
        s_r  = f"{ratio_P3P2:.4f}" if ratio_P3P2 else "-"
        print(f"{alpha:>7.1f} {a_gam:>7.3f} {P2:>10.5f} {ratio_2pi:>12.6f} {diff_pct:>8.3f}% {s_P3:>10} {s_r:>14}")
        results_pi.append({'alpha': alpha, 'a_gam': a_gam,
                           'P2': P2, 'P3': P3, 'M1': M1})

# ======================================================================
# 3. LIMIT a_gam -> 0: czy Phi0_2 -> 2*pi?
# ======================================================================
print()
print("="*65)
print("LIMIT a_gam -> 0 (alpha=5.9):")
print("="*65)
print(f"{'a_gam':>8} {'Phi0_2':>10} {'Phi0_2/2pi':>12} {'Phi0_3':>10}")
for a_gam in [0.10, 0.08, 0.05, 0.03, 0.02, 0.01]:
    M1 = find_root(1.0, 5.9, a_gam)
    if M1 is None: continue
    P2 = find_Phi0_for_ratio(207.0,  M1, 5.9, a_gam)
    P3 = find_Phi0_for_ratio(3477.0, M1, 5.9, a_gam)
    s_P3 = f"{P3:.5f}" if P3 else "-"
    if P2:
        print(f"{a_gam:>8.3f} {P2:>10.5f} {P2/(2*np.pi):>12.7f} {s_P3:>10}")

# ======================================================================
# 4. POSZUKIWANIE POSTACI Phi0_3
# ======================================================================
print()
print("="*65)
print("POSTAC Phi0_3 — sprawdzanie hipotez:")
print("="*65)
# Uzyj wynikow dla a_gam=0.05, alpha=5.9
r5 = [r for r in results_pi if r['alpha']==5.9 and r['a_gam']==0.05]
if r5:
    P2 = r5[0]['P2']; P3 = r5[0]['P3']
    if P3:
        print(f"\n  alpha=5.9, a_gam=0.05:")
        print(f"  Phi0_2 = {P2:.6f}")
        print(f"  Phi0_3 = {P3:.6f}")
        print()
        # Hipotezy dla Phi0_3:
        hypotheses = {
            "4*pi":          4*np.pi,
            "(2*pi)^(3/2)":  (2*np.pi)**1.5,
            "8*pi - 2":      8*np.pi - 2,
            "4*pi^2/pi":     4*np.pi,  # = 4pi
            "2*pi*sqrt(Phi0_2)": 2*np.pi*np.sqrt(P2),
            "Phi0_2^2/(2*pi)":   P2**2/(2*np.pi),
            "Phi0_2 * pi":       P2 * np.pi,
            "Phi0_2 * 4":        P2 * 4.0,
            "Phi0_2 * e":        P2 * np.e,
            "sqrt(207)*pi":      np.sqrt(207)*np.pi,
            "pi^3":              np.pi**3,
            "2*pi^2":            2*np.pi**2,
            "3*pi^2/sqrt(pi)":   3*np.pi**2/np.sqrt(np.pi),
        }
        for name, val in hypotheses.items():
            diff = 100*(P3 - val)/P3
            marker = "  <<< BLISKO!" if abs(diff) < 2.0 else ""
            print(f"  {name:30s} = {val:8.4f}  (diff {diff:+.2f}%){marker}")

# ======================================================================
# 5. KOIDE FORMULA — czy TGP moze ja wyprowadzic?
# ======================================================================
print()
print("="*65)
print("KOIDE FORMULA: (me+mu+mt)/(sqrt(me)+sqrt(mu)+sqrt(mt))^2 = 2/3")
print("="*65)
print()
# Znane wartosci eksperymentalne [MeV]
me, mmu, mta = 0.51099895, 105.6583755, 1776.86
K_koide_exp = (me+mmu+mta) / (np.sqrt(me)+np.sqrt(mmu)+np.sqrt(mta))**2
print(f"  Eksperyment: K = {K_koide_exp:.8f}  (= 2/3 = {2/3:.8f})")
print()
# Sprawdz: jezeli masy sa proporcjonalne 1:207:3477
m1, m2, m3 = 1.0, 207.0, 3477.0
K_tgp = (m1+m2+m3) / (np.sqrt(m1)+np.sqrt(m2)+np.sqrt(m3))**2
print(f"  Stosunki TGP (1:207:3477): K = {K_tgp:.8f}  (cel 2/3 = {2/3:.8f})")
diff_koide = 100*(K_tgp - 2/3)/(2/3)
print(f"  Roznica od 2/3: {diff_koide:+.4f}%")
print()

# Koide dla naszych numerycznych wartosci Phi0-based
if r5 and r5[0]['P3']:
    M1n = r5[0]['M1']
    M2n = M1n * 207.0
    M3n = M1n * 3477.0
    K_num = (M1n+M2n+M3n)/(np.sqrt(M1n)+np.sqrt(M2n)+np.sqrt(M3n))**2
    print(f"  Masy TGP (alpha=5.9): K = {K_num:.8f}")
    print()

# Pytanie: dla jakich stosunkow r21, r31 zachodzi Koide?
# (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2 = 2/3
# => 3*(1 + r21 + r31) = 2*(1 + sqrt(r21) + sqrt(r31))^2
# Sprawdzamy przy r21=207, jaki r31 spelnia Koide?
def koide_residual(r31, r21=207.0):
    num = 1 + r21 + r31
    den = (1 + np.sqrt(r21) + np.sqrt(r31))**2
    return num/den - 2/3

try:
    r31_koide = brentq(koide_residual, 100.0, 100000.0, xtol=1e-6)
    print(f"  Przy r21=207: Koide wymaga r31 = {r31_koide:.2f}  (eksperyment: 3477.0)")
    print(f"  Roznica od 3477: {100*(r31_koide-3477)/3477:+.2f}%")
except Exception:
    print("  Nie znaleziono rozwiazania Koide dla r21=207")

print()
print("="*65)
print("PODSUMOWANIE OBSERWACJI PI:")
print("="*65)
print()
print("  Lambda = 4*pi  (zawsze, z geometrii)")
print(f"  Phi0_2 / (2*pi) jest bliskie 1 dla wielu (alpha, a_gam)")
print(f"  Zalezy od alpha i a_gam => to nie jest DOKLADNE 2*pi")
print(f"  ALE: w limicie a_gam->0 moze konwergowac do 2*pi?")
print()
print("  Pi pojawia sie NATURALNIE przez Lambda = 4*pi.")
print("  Phi0_2 ~ Lambda/2 = 2*pi byloby eleganckim wynikiem:")
print("    'tlo drugiej generacji = polowa stalej samospojnosci'")
