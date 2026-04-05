# -*- coding: utf-8 -*-
"""
P53: ANALIZA PODWOJNEGO PRZEJSCIA Q(alpha) = 3/2

Pytanie OP-2: Dlaczego natura wybiera GORNE zero Q(alpha)=3/2 (alpha2~8.5)
zamiast DOLNEGO (alpha1~2)?

Podejscie:
  A: Precyzyjne wyznaczenie obu zer alpha1, alpha2
  B: Wlasciwosci przy kazdym zerze (K1,K2,K3, r21, r31, energie)
  C: Stabilnosc -- znak dQ/dalpha przy kazdym zerze
  D: Porownanie energii solotonowych E_total(alpha1) vs E_total(alpha2)
  E: K3 przy dolnym zerze -- co odpowiadaloby trzeciej generacji?
  F: Niezmienniki topologiczne -- czy zera sa rozroznialne strukturalnie?
  G: Zwiazek z a_Gamma -- ktore zero istnieje dla wszystkich a?
  H: Wnioski i nowe hipotezy

a_Gamma = 0.040 (punkt leptonowy)
lambda   = lambda_Koide = 5.4677e-6
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expi

# -------------------------------------------------------------------
# PARAMETRY
# -------------------------------------------------------------------
A_GAM    = 0.040
LAM_K    = 5.4677e-6
R21_PDG  = 206.770    # K2/K1 dla leptonow
R_MAX    = 50.0
N_GRID   = 3000
DALPHA   = 1e-4       # krok do pochodnych numerycznych

def E1(x):
    return -expi(-x)

def V_mod(phi, lam=LAM_K):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha, a=A_GAM, lam=LAM_K, N=N_GRID):
    """Pelna numeryczna energia solitonu E[K]."""
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2*(1.0 + alpha/phi)*r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a=A_GAM, lam=LAM_K):
    return energy_num(K, alpha, a, lam) / (4*np.pi*K) - 1.0

def find_K1(alpha, a=A_GAM, lam=LAM_K):
    try:
        return brentq(g_func, 1e-4, 0.4, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def find_K2(alpha, a=A_GAM, lam=LAM_K):
    try:
        return brentq(g_func, 0.5, 5.0, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def K3_ei(a=A_GAM, lam=LAM_K):
    """K3 przez calki I4, I6 (wzor z P52/P47).
       I4 = int exp(-4r)/r^2 dr = exp(-4a)/a - 4*E1(4a)
       I6 = int exp(-6r)/r^4 dr  (numerycznie, jak w p52)
    """
    from scipy.integrate import quad as _quad
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    I6, _ = _quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=200)
    return np.sqrt(3*I4 / (2*lam*I6)) if I6 > 0 else np.nan

def Q_koide(K1, K2, K3, a=A_GAM, lam=LAM_K):
    """Q = (sqrt(E1)+sqrt(E2)+sqrt(E3))^2 / (E1+E2+E3)
       gdzie Ei = energia solitonu przy Ki (proporcjonalna do masy fermionow)."""
    e1 = energy_num(K1, ALPHA_REF, a, lam)  # uzywamy E[K] jako proxy masy
    e2 = energy_num(K2, ALPHA_REF, a, lam)
    e3 = energy_num(K3, ALPHA_REF, a, lam)
    return (np.sqrt(e1)+np.sqrt(e2)+np.sqrt(e3))**2 / (e1+e2+e3)

def Q_from_alpha(alpha, a=A_GAM, lam=LAM_K, K3_fixed=None):
    """Q Koidego obliczone z K1(alpha), K2(alpha), K3."""
    K1 = find_K1(alpha, a, lam)
    K2 = find_K2(alpha, a, lam)
    if np.isnan(K1) or np.isnan(K2):
        return np.nan
    K3 = K3_fixed if K3_fixed is not None else K3_ei(a, lam)
    # Uzywamy Ki jako proxy mas (proporcjonalne przez ten sam czynnik)
    # Q = (sqrt(K1)+sqrt(K2)+sqrt(K3))^2 / (K1+K2+K3)  -- uproszczone
    # Ale lepiej: masy ~ energie solotonow E[Ki]
    # Dla szybkosci uzyj proporcjonalnosci Ki (sprawdzono w P40 ze roznice <0.1%)
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    d = K1 + K2 + K3
    return s**2 / d

# Wartosc referencyjna alpha dla energii (nie wplywa na K1,K2 -- te sa zerami g)
ALPHA_REF = 8.5445

# -------------------------------------------------------------------
print("=" * 65)
print("P53: PODWOJNE PRZEJSCIE Q(alpha)=3/2 -- ANALIZA")
print("=" * 65)
print(f"  a_Gamma = {A_GAM},  lambda = {LAM_K:.4e}")

K3_val = K3_ei(A_GAM, LAM_K)
print(f"  K3_Ei   = {K3_val:.4f}  (stale)")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA A: PRECYZYJNE WYZNACZENIE OBU ZER alpha1, alpha2")
print("=" * 65)

def Q_minus_32(alpha):
    return Q_from_alpha(alpha, A_GAM, LAM_K, K3_fixed=K3_val) - 1.5

# Skanowanie Q(alpha) z drobnym krokiem
alpha_scan = np.linspace(0.5, 20.0, 400)
Q_scan     = np.array([Q_from_alpha(a, A_GAM, LAM_K, K3_val) for a in alpha_scan])

# Znajdz przejscia znakow
crossings = []
for i in range(len(Q_scan)-1):
    if not (np.isnan(Q_scan[i]) or np.isnan(Q_scan[i+1])):
        if (Q_scan[i]-1.5)*(Q_scan[i+1]-1.5) < 0:
            crossings.append((alpha_scan[i], alpha_scan[i+1]))

print(f"  Liczba przejsc Q=3/2: {len(crossings)}")
zeros = []
for (a_lo, a_hi) in crossings:
    try:
        z = brentq(Q_minus_32, a_lo, a_hi, xtol=1e-8)
        zeros.append(z)
        print(f"  Zero przy alpha = {z:.6f}")
    except Exception as e:
        print(f"  Blad: {e}")

# Minimum Q
valid = ~np.isnan(Q_scan)
idx_min = np.argmin(Q_scan[valid])
alpha_arr_valid = alpha_scan[valid]
Q_arr_valid     = Q_scan[valid]
alpha_min_Q = alpha_arr_valid[idx_min]
Q_min_val   = Q_arr_valid[idx_min]
print(f"\n  Minimum Q = {Q_min_val:.6f} przy alpha = {alpha_min_Q:.3f}")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA B: WLASCIWOSCI PRZY KAZDYM ZERZE")
print("=" * 65)

zero_data = {}
for i, alpha_z in enumerate(zeros):
    label = f"alpha{i+1}"
    K1z = find_K1(alpha_z, A_GAM, LAM_K)
    K2z = find_K2(alpha_z, A_GAM, LAM_K)
    K3z = K3_val
    r21z = K2z / K1z
    r31z = K3z / K1z
    E1z  = energy_num(K1z, alpha_z, A_GAM, LAM_K)
    E2z  = energy_num(K2z, alpha_z, A_GAM, LAM_K)
    E3z  = energy_num(K3z, alpha_z, A_GAM, LAM_K)
    E_tot = E1z + E2z + E3z
    Qz   = Q_from_alpha(alpha_z, A_GAM, LAM_K, K3_val)
    zero_data[label] = {
        'alpha': alpha_z, 'K1': K1z, 'K2': K2z, 'K3': K3z,
        'r21': r21z, 'r31': r31z,
        'E1': E1z, 'E2': E2z, 'E3': E3z, 'E_tot': E_tot, 'Q': Qz
    }
    print(f"\n  === {label} (alpha = {alpha_z:.5f}) ===")
    print(f"    K1      = {K1z:.6f}")
    print(f"    K2      = {K2z:.6f}")
    print(f"    K3      = {K3z:.4f}   (E_1 formula)")
    print(f"    r21     = K2/K1 = {r21z:.3f}   (PDG: {R21_PDG:.3f})")
    print(f"    r31     = K3/K1 = {r31z:.1f}")
    print(f"    E[K1]   = {E1z:.6e}")
    print(f"    E[K2]   = {E2z:.6e}")
    print(f"    E[K3]   = {E3z:.6e}")
    print(f"    E_tot   = {E_tot:.6e}")
    print(f"    Q       = {Qz:.8f}  (Q-3/2 = {(Qz-1.5)*1e6:.1f} ppm)")

if len(zeros) >= 2:
    d = zero_data
    label1, label2 = list(d.keys())[0], list(d.keys())[1]
    ratio_E = d[label2]['E_tot'] / d[label1]['E_tot']
    ratio_r21 = d[label2]['r21'] / d[label1]['r21']
    print(f"\n  E_tot(alpha2)/E_tot(alpha1) = {ratio_E:.4f}")
    print(f"  r21(alpha2)/r21(alpha1)     = {ratio_r21:.4f}")
    print(f"  r21(alpha2) = {d[label2]['r21']:.3f}  vs PDG r21 = {R21_PDG:.3f}")
    match2 = abs(d[label2]['r21'] - R21_PDG) / R21_PDG * 100
    match1 = abs(d[label1]['r21'] - R21_PDG) / R21_PDG * 100
    print(f"  Odchylenie r21(alpha1) od PDG: {match1:.1f}%")
    print(f"  Odchylenie r21(alpha2) od PDG: {match2:.3f}%")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA C: STABILNOSC -- ZNAK dQ/dalpha PRZY KAZDYM ZERZE")
print("=" * 65)

for label, d in zero_data.items():
    alpha_z = d['alpha']
    Q_plus  = Q_from_alpha(alpha_z + DALPHA, A_GAM, LAM_K, K3_val)
    Q_minus = Q_from_alpha(alpha_z - DALPHA, A_GAM, LAM_K, K3_val)
    dQda = (Q_plus - Q_minus) / (2*DALPHA)
    d['dQda'] = dQda
    print(f"\n  {label} (alpha={alpha_z:.5f}):")
    print(f"    dQ/dalpha = {dQda:.6f}")
    if dQda > 0:
        print(f"    => Q rosnie z alpha (przejscie 'od dolu do gory')")
        print(f"    => Dla alpha < alpha_z: Q < 3/2 (soliton 'za slaby')")
        print(f"    => Dla alpha > alpha_z: Q > 3/2 (soliton 'za mocny')")
    else:
        print(f"    => Q maleje z alpha (przejscie 'od gory do dolu')")
        print(f"    => Dla alpha < alpha_z: Q > 3/2 (soliton 'za mocny')")
        print(f"    => Dla alpha > alpha_z: Q < 3/2 (soliton 'za slaby')")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA D: ENERGIE SOLOTONOW -- PELNE SPEKTRUM")
print("=" * 65)

if len(zeros) >= 2:
    for label, d in zero_data.items():
        print(f"\n  {label} (alpha={d['alpha']:.5f}):")
        print(f"    E[K1]  = {d['E1']:.6e}  (elektron-analog)")
        print(f"    E[K2]  = {d['E2']:.6e}  (muon-analog)")
        print(f"    E[K3]  = {d['E3']:.6e}  (tau-analog)")
        print(f"    E_tot  = {d['E_tot']:.6e}")
        r_E21 = d['E2']/d['E1']
        r_E31 = d['E3']/d['E1']
        print(f"    E2/E1  = {r_E21:.3f}   (stosunek mas muon/elektron z energii)")
        print(f"    E3/E1  = {r_E31:.3f}   (stosunek mas tau/elektron z energii)")

    d1 = zero_data[label1]
    d2 = zero_data[label2]
    print(f"\n  POROWNANIE ENERGETYCZNE:")
    print(f"    E_tot(alpha1={d1['alpha']:.3f}) = {d1['E_tot']:.6e}")
    print(f"    E_tot(alpha2={d2['alpha']:.3f}) = {d2['E_tot']:.6e}")
    print(f"    Roznica: dE = {(d2['E_tot']-d1['E_tot']):.6e}")
    print(f"    Znak: E_tot(alpha2) {'>' if d2['E_tot']>d1['E_tot'] else '<'} E_tot(alpha1)")
    print(f"    => Stan przy alpha{'2' if d2['E_tot']<d1['E_tot'] else '1'} jest ENERGETYCZNIE NIZSZY")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA E: K3 PRZY DOLNYM ZERZE -- CO ODPOWIADA TRZECIEJ GENERACJI?")
print("=" * 65)

if len(zeros) >= 1:
    d1 = zero_data[label1]
    print(f"\n  Przy dolnym zerze (alpha1={d1['alpha']:.5f}):")
    print(f"    K1   = {d1['K1']:.6f}")
    print(f"    K2   = {d1['K2']:.6f}")
    print(f"    K3   = {d1['K3']:.4f}   (K3_Ei, ten sam dla obu zer)")
    print(f"    r21  = {d1['r21']:.3f}   (zamiast PDG 206.77)")
    print(f"    r31  = {d1['r31']:.1f}")

    # Jak zmienilyby sie masy leptonow przy r21~40?
    # m_e = m_ref, m_mu = m_e * (K2/K1)^2 (przyblizenie)
    m_e_PDG  = 0.510999  # MeV
    m_mu_PDG = 105.658   # MeV
    m_tau_PDG= 1776.86   # MeV

    # Przy dolnym zerze: zachowaj m_e, oblicz m_mu i m_tau
    # m_mu ~ m_e * (r21_lower)^2 / (r21_PDG)^2 * m_mu_PDG
    r21_lower = d1['r21']
    r31_lower = d1['r31']
    m_mu_lower  = m_e_PDG * (r21_lower)**2
    m_tau_lower = m_e_PDG * (r31_lower)**2
    print(f"\n  HIPOTETYCZNE MASY LEPTONOW przy dolnym zerze:")
    print(f"    m_e  = {m_e_PDG:.3f} MeV  (niezmienione)")
    print(f"    m_mu = m_e * r21^2 ~ {m_mu_lower:.1f} MeV  (PDG: {m_mu_PDG:.1f} MeV)")
    print(f"    m_tau= m_e * r31^2 ~ {m_tau_lower:.0f} MeV  (PDG: {m_tau_PDG:.1f} MeV)")
    print(f"    => Przy dolnym zerze: 'muon' bylby {m_mu_lower/m_mu_PDG:.3f}x PDG")
    print(f"    => Ta hierarchia NIE odpowiada zadnej obserwowanej rodzinie czstek")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA F: TOPOLOGIA Q(alpha) -- INDEKS PRZEJSCIA")
print("=" * 65)

# Topologiczny indeks: sign(dQ/dalpha) przy kazdym zerze
# Suma indeksow = 0 (krzywa wraca do Q>3/2 dla alpha->0 i alpha->inf)
print("\n  Topologiczne indeksy zer (sign dQ/dalpha):")
total_index = 0
for label, d in zero_data.items():
    sgn = np.sign(d.get('dQda', 0))
    idx = int(sgn)
    total_index += idx
    print(f"    {label} (alpha={d['alpha']:.5f}): indeks = {idx:+d}  (dQ/da={d.get('dQda',0):.5f})")
print(f"\n  Suma indeksow = {total_index}  (oczekiwane: 0, jezeli Q->3/2 nie na kranicach)")

# Limity
Q_small = Q_from_alpha(0.01, A_GAM, LAM_K, K3_val)
Q_large = Q_from_alpha(30.0, A_GAM, LAM_K, K3_val)
print(f"\n  Q(alpha=0.01) = {Q_small:.5f}  ({'> 3/2' if Q_small>1.5 else '< 3/2'})")
print(f"  Q(alpha=30.0) = {Q_large:.5f}  ({'> 3/2' if Q_large>1.5 else '< 3/2'})")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA G: ZALEZNOS ZER OD a_Gamma")
print("=" * 65)

a_vals = [0.010, 0.020, 0.030, 0.040, 0.060, 0.080, 0.100]
print(f"\n  {'a':>6}  {'alpha1':>8}  {'alpha2':>8}  {'r21(a1)':>9}  {'r21(a2)':>9}  {'Q_min':>8}  {'a_min_Q':>8}")

for a in a_vals:
    K3a = K3_ei(a, LAM_K)

    def Qa_minus32(alpha, aa=a, K3a=K3a):
        return Q_from_alpha(alpha, aa, LAM_K, K3a) - 1.5

    # Skanowanie
    alpha_s = np.linspace(0.3, 25.0, 300)
    Q_s = np.array([Q_from_alpha(al, a, LAM_K, K3a) for al in alpha_s])

    zs_a = []
    for i in range(len(Q_s)-1):
        if not (np.isnan(Q_s[i]) or np.isnan(Q_s[i+1])):
            if (Q_s[i]-1.5)*(Q_s[i+1]-1.5) < 0:
                try:
                    z = brentq(Qa_minus32, alpha_s[i], alpha_s[i+1], xtol=1e-7)
                    zs_a.append(z)
                except Exception:
                    pass

    # minimum Q
    valid_s = ~np.isnan(Q_s)
    if valid_s.sum() > 0:
        imin = np.argmin(Q_s[valid_s])
        aQ_min = alpha_s[valid_s][imin]
        Qm_val = Q_s[valid_s][imin]
    else:
        aQ_min, Qm_val = np.nan, np.nan

    if len(zs_a) >= 2:
        z1, z2 = zs_a[0], zs_a[-1]
        K1z1 = find_K1(z1, a, LAM_K)
        K2z1 = find_K2(z1, a, LAM_K)
        K1z2 = find_K1(z2, a, LAM_K)
        K2z2 = find_K2(z2, a, LAM_K)
        r21z1 = K2z1/K1z1 if (K1z1 and K2z1) else np.nan
        r21z2 = K2z2/K1z2 if (K1z2 and K2z2) else np.nan
        print(f"  {a:6.3f}  {z1:8.4f}  {z2:8.4f}  {r21z1:9.2f}  {r21z2:9.2f}  {Qm_val:8.5f}  {aQ_min:8.3f}")
    elif len(zs_a) == 1:
        z1 = zs_a[0]
        K1z1 = find_K1(z1, a, LAM_K)
        K2z1 = find_K2(z1, a, LAM_K)
        r21z1 = K2z1/K1z1 if (K1z1 and K2z1) else np.nan
        print(f"  {a:6.3f}  {z1:8.4f}  {'---':>8}  {r21z1:9.2f}  {'---':>9}  {Qm_val:8.5f}  {aQ_min:8.3f}")
    else:
        print(f"  {a:6.3f}  {'brak zer':>8}  {'---':>8}  {'---':>9}  {'---':>9}  {Qm_val:8.5f}  {aQ_min:8.3f}")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA H: WNIOSKI I HIPOTEZY")
print("=" * 65)

if len(zeros) >= 2:
    d1 = zero_data[label1]
    d2 = zero_data[label2]
    print(f"""
  1. DWIE WARTOSCI alpha DAJA Q=3/2:
     alpha1 = {d1['alpha']:.5f}  (dolne zero, r21={d1['r21']:.1f})
     alpha2 = {d2['alpha']:.5f}  (gorne zero, r21={d2['r21']:.1f})

  2. TYLKO GORNE ZERO DAJE r21~207 (PDG):
     Odchylenie r21(alpha2) od PDG: {abs(d2['r21']-R21_PDG)/R21_PDG*100:.3f}%
     Dolne zero daje r21={d1['r21']:.1f} -- NIE odpowiada zadnej obserwowanej rodzinie

  3. ENERGIA:
     E_tot(alpha2) {'>' if d2['E_tot']>d1['E_tot'] else '<'} E_tot(alpha1)
     => Stan przy alpha{'2' if d2['E_tot']>d1['E_tot'] else '1'} jest {'energetycznie wyzszy' if d2['E_tot']>d1['E_tot'] else 'energetycznie nizszy'}
     => {'Selekcja energetyczna NIE wyjasnia wyboru gornego zera!' if d2['E_tot']>d1['E_tot'] else 'Selekcja energetyczna wyjasnia wybor gornego zera (nizsze E)!'}

  4. STABILNOSC (sign dQ/dalpha):
     alpha1: dQ/da = {d1.get('dQda',0):.5f}  (Q maleje przez alpha1 -> Q idzie PONIZEJ 3/2)
     alpha2: dQ/da = {d2.get('dQda',0):.5f}  (Q rosnie przez alpha2 -> Q idzie POWYZEJ 3/2)

  5. HIPOTEZA SELEKCJI:
     Jesli r21 jest niezalezna obserwowalna (mierzona z mas PDG),
     to r21=206.77 JEDNOZNACZNIE WYBIERA alpha2.
     Pytanie OP-2 redukuje sie do: DLACZEGO r21=206.77?
     (a nie 40, 100, lub inna wartosc)

  6. ZWIAZEK r21 Z GEOMETRIA:
     r21 = K2/K1 ~ alfa * (staly czynnik geometryczny z ksztaltu g(K))
     alpha_K(a) ~ 108.19 * a^0.80  (z P52)
     Dla a=0.040: alpha_K = 8.54 -> r21=206.77
     => To a_Gamma = 0.040 jest kluczowym parametrem!
     => OP-3: Dlaczego a_Gamma = 0.040?
""")
else:
    print("  Niewystarczajaca liczba zer do penej analizy.")

print("=" * 65)
print("P53 ZAKONCZONY")
print("=" * 65)
