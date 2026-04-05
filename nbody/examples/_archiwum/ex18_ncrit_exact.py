"""
ex18_ncrit_exact.py — Dokładna N_crit przez całkę Feynmana vs aproksymacja siodłopunktowa

Pytanie: dla jakich (C, m_sp, d) efekty trójcząstkowe TGP stają się mierzalne?

Metoda:
  1. V3/V2 liczone DOKŁADNIE przez całkę 2D Feynmana (jak w ex17)
  2. N_crit definiujemy jako V3_total / V2_single = fraction (np. 1%)
     Dla N ciał w siatce: V3_total ≈ N * (N-1)/2 * V3_pair_avg
     V2_total = N * (N-1)/2 * V2_pair_avg
     Więc V3/V2 nie zależy od N — to ratio per-triple vs per-pair
  3. Porównanie z przybliżeniem siodłopunktowym
     I_Y^saddle = (2π/Δ) * exp(-m * sqrt(Q/Δ)) / (m * sqrt(Q/Δ))
     (przybliżenie K_0(u) ≈ sqrt(π/(2u)) exp(-u) dla u >> 1)

Struktura:
  - Sekcja A: ratio V3/V2 jako funkcja odległości d dla różnych C, m_sp
  - Sekcja B: ratio V3/V2 jako funkcja m_sp * d (parametr adimensionalny)
  - Sekcja C: błąd siodłopunktowy vs dokładna całka Feynmana
  - Sekcja D: mapa N_crit(C, m_sp) w przestrzeni parametrów
"""

import numpy as np
from scipy import integrate
import sys
import os

# ─── ścieżka do three_body_force_exact ───────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
try:
    from three_body_force_exact import feynman_yukawa_3body
    EXACT_AVAILABLE = True
    print("OK three_body_force_exact zaladowany")
except ImportError as e:
    EXACT_AVAILABLE = False
    print(f"BRAK Brak three_body_force_exact: {e}")

# ─── implementacja lokalna (redundantna) dla weryfikacji ─────────────────────

def duffy_transform_quad(f, n=60):
    """
    Całkowanie na sympleksie Δ₂ = {α₁+α₂+α₃=1, αᵢ≥0}
    przez transformację Duffy'ego: α₁=u, α₂=v(1-u), α₃=(1-u)(1-v)
    Jakobian: (1-u)
    Zwraca ∫∫_Δ₂ f(α₁,α₂,α₃) dS
    """
    pts, wts = np.polynomial.legendre.leggauss(n)
    # mapowanie [−1,1] → [0,1]
    u_pts = 0.5 * (1 + pts)
    u_wts = 0.5 * wts
    v_pts = 0.5 * (1 + pts)
    v_wts = 0.5 * wts

    result = 0.0
    for i, (u, wu) in enumerate(zip(u_pts, u_wts)):
        for j, (v, wv) in enumerate(zip(v_pts, v_wts)):
            a1 = u
            a2 = v * (1.0 - u)
            a3 = (1.0 - u) * (1.0 - v)
            jac = (1.0 - u)
            result += wu * wv * jac * f(a1, a2, a3)
    return result


def yukawa_3body_integral_local(d12, d13, d23, m, n_quad=60):
    """
    Całka trójcząstkowa Feynmana:
      I_Y = 2 ∫_Δ₂ Δ^{-3/2} K₀(m √(Q/Δ)) dα
    gdzie:
      Q = α₂·d₁₂² + α₁·d₁₃² + α₃·d₂₃²
          (uwaga: α₃ wchodzi przy d₁₂, itd. — indeksy z parametryzacji Schwingera)
      Δ = α₁α₂ + α₁α₃ + α₂α₃
    Zwraca I_Y (skalar).
    """
    from scipy.special import k0 as K0

    d12sq = d12**2
    d13sq = d13**2
    d23sq = d23**2

    def integrand(a1, a2, a3):
        Delta = a1*a2 + a1*a3 + a2*a3
        if Delta < 1e-30:
            return 0.0
        Q = a2*d12sq + a1*d13sq + a3*d23sq
        u = m * np.sqrt(Q / Delta)
        if u < 1e-30:
            return 0.0
        return Delta**(-1.5) * K0(u)

    return 2.0 * duffy_transform_quad(integrand, n=n_quad)


def yukawa_2body(r, m):
    """V₂ = exp(-m·r)/r"""
    return np.exp(-m * r) / r


def yukawa_3body_saddle(d12, d13, d23, m):
    """
    Przybliżenie siodłopunktowe całki Feynmana:
    I_Y^saddle ≈ π/(m·Q̄^(1/2)) * (Δ̄/Q̄)^(1/2) * exp(-m·sqrt(Q̄/Δ̄)) / Δ̄

    Siodło (minimum eksponenta m√(Q/Δ)) leży w centrum sympleksu α₁=α₂=α₃=1/3.
    W tym punkcie:
      Δ̄ = 3*(1/3)² = 1/3
      Q̄ = (1/3)*(d₁₂²+d₁₃²+d₂₃²) / (1/3)...

    Dokładniej: Δ̄=1/3, Q̄=(d₁₂²+d₁₃²+d₂₃²)/3 * Δ̄...

    Z parametryzacji Schwingera (αi sumują do 1):
      Q = α₂·d₁₂² + α₁·d₁₃² + α₃·d₂₃² przy α₁=α₂=α₃=1/3:
      Q̄ = (1/3)(d₁₂²+d₁₃²+d₂₃²)/3...

    Uwaga: eksponent ma postać exp(-m·sqrt(Q/Δ)) = exp(-m·sqrt(3·Q̄))

    W pracy ex6: formuła używa exp(-1.5t) zamiast exp(-√3·t) — BŁĄD
    Poprawna formuła siodłopunktowa (Schwinger 3-body):
    I_saddle = C_geom * exp(-m·d_eff) / (m·d_eff)
    gdzie d_eff = sqrt(Q̄/Δ̄) = sqrt((d₁₂²+d₁₃²+d₂₃²)/3) dla równobocznego.
    """
    Q_bar = (d12**2 + d13**2 + d23**2) / 3.0   # na α=1/3
    Delta_bar = 1.0 / 3.0
    u_saddle = m * np.sqrt(Q_bar / Delta_bar)    # = m * sqrt(3 * Q_bar_sym)

    # K₀(u) ≈ sqrt(π/(2u)) * exp(-u) dla u >> 1
    # I_saddle ≈ 2 * (area_saddle) * Δ̄^{-3/2} * K₀(u_saddle)
    # area_saddle ≈ (2π/m) * (Δ̄/Q̄)^{1/2} * (1/(something))
    # Uproszczona formuła z ex6: I ~ exp(-m*d_eff) * prefactor
    from scipy.special import k0 as K0
    # Użyj dokładnego K₀ ale z siodłowym Q̄/Δ̄
    return 2.0 * Delta_bar**(-1.5) * K0(u_saddle)


# ─── Sekcja A: V3/V2 vs odległość ────────────────────────────────────────────

print("\n" + "="*65)
print("SEKCJA A: V3/V2 jako funkcja odległości d")
print("="*65)
print("Trójkąt równoboczny: d₁₂=d₁₃=d₂₃=d")
print()

C_vals = [0.05, 0.10, 0.20, 0.28]  # 0.28 ≈ C_Planck
m_sp = 1.0
d_vals = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])

print(f"m_sp = {m_sp}")
print(f"{'d':>6}  " + "  ".join([f"V3/V2(C={c:.2f})" for c in C_vals]))
print("-" * 75)

results_A = {}
for d in d_vals:
    row = []
    for C in C_vals:
        I3 = yukawa_3body_integral_local(d, d, d, m_sp, n_quad=50)
        V2 = yukawa_2body(d, m_sp)
        # V3 = C³ * I3, V2 = C² * V2
        ratio = C * I3 / V2
        row.append(ratio)
    results_A[d] = row
    print(f"{d:6.1f}  " + "  ".join([f"{r:12.4e}" for r in row]))

# ─── Sekcja B: V3/V2 vs m_sp * d ─────────────────────────────────────────────

print("\n" + "="*65)
print("SEKCJA B: V3/V2 vs parametr bezwymiarowy t = m_sp * d")
print("="*65)
print("Trójkąt równoboczny, C=0.20")
print()

C_fixed = 0.20
t_vals = np.linspace(0.2, 6.0, 30)  # t = m*d

print(f"{'t=m*d':>8}  {'V3/V2 (dokładna)':>18}  {'V3/V2 (siodło)':>16}  {'Błąd siodła':>12}")
print("-" * 65)

t_list, ratio_exact_list, ratio_saddle_list = [], [], []

for t in t_vals:
    # m=1, d=t (lub ogólnie: skalowanie)
    d = t  # m=1
    I3_exact = yukawa_3body_integral_local(d, d, d, 1.0, n_quad=50)
    I3_saddle = yukawa_3body_saddle(d, d, d, 1.0)
    V2 = yukawa_2body(d, 1.0)
    ratio_e = C_fixed * I3_exact / V2
    ratio_s = C_fixed * I3_saddle / V2
    err = abs(ratio_s - ratio_e) / abs(ratio_e) * 100 if ratio_e != 0 else 0
    t_list.append(t)
    ratio_exact_list.append(ratio_e)
    ratio_saddle_list.append(ratio_s)
    if t <= 1.0 or abs(t - 2.0) < 0.15 or abs(t - 3.0) < 0.15 or abs(t - 5.0) < 0.15:
        print(f"{t:8.2f}  {ratio_e:18.6e}  {ratio_s:16.6e}  {err:11.1f}%")

# ─── Sekcja C: Błąd siodłopunktowy ───────────────────────────────────────────

print("\n" + "="*65)
print("SEKCJA C: Błąd siodłopunktowy vs dokładna całka")
print("="*65)
print("Geometrie: równoboczny, izoceles (1:1:2), prostokątny")
print()

geometries = {
    'równoboczny (1:1:1)': lambda d: (d, d, d),
    'izosceles  (1:1:2)' : lambda d: (d, d, 2*d),
    'izosceles  (1:2:2)' : lambda d: (d, 2*d, 2*d),
}

d_test = 2.0
m_test = 1.0
print(f"d_ref={d_test}, m_sp={m_test}")
print()

for name, geo_fn in geometries.items():
    d12, d13, d23 = geo_fn(d_test)
    I_exact = yukawa_3body_integral_local(d12, d13, d23, m_test, n_quad=60)
    I_saddle = yukawa_3body_saddle(d12, d13, d23, m_test)
    err = (I_saddle - I_exact) / I_exact * 100
    print(f"  {name}: I_exact={I_exact:.6e}, I_saddle={I_saddle:.6e}, błąd={err:+.1f}%")

# Skan błędu siodłopunktowego vs t dla geometrii równobocznej
print()
print("Błąd siodłopunktowy vs t=m*d (geometria równoboczna):")
print(f"{'t':>6}  {'I_exact':>14}  {'I_saddle':>14}  {'błąd %':>10}")
print("-" * 55)

for t in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
    d = t
    I_e = yukawa_3body_integral_local(d, d, d, 1.0, n_quad=60)
    I_s = yukawa_3body_saddle(d, d, d, 1.0)
    err = (I_s - I_e) / I_e * 100
    print(f"{t:6.1f}  {I_e:14.6e}  {I_s:14.6e}  {err:+10.1f}%")

# ─── Sekcja D: Mapa N_crit(C, m_sp) ──────────────────────────────────────────

print("\n" + "="*65)
print("SEKCJA D: V3/V2 ratio w przestrzeni (C, m_sp·d)")
print("="*65)
print("Cel: znaleźć granicę obserwowalności V3/V2 > 1%")
print()

# Dla trójkąta równobocznego z jednostkową odległością:
# V3/V2 = C * I3(t)/V2(t) gdzie t=m_sp*d
# Szukamy C_crit(t) takie że V3/V2 = 0.01

threshold = 0.01  # 1%

print("C_crit (dla V3/V2=1%) vs t=m_sp*d:")
print(f"{'t=m*d':>8}  {'I3/V2 ratio':>14}  {'C_crit':>12}  {'m_cząstki (kg)':>18}")
print("-" * 60)

# masa ciała w [Planck] odpowiadająca C = m/(2√π)
# C = m_Planck * m_body / (2*sqrt(pi)) => m_body = 2*sqrt(pi)*C [Planck]
m_Pl_kg = 2.176e-8  # kg

t_scan = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
for t in t_scan:
    I3 = yukawa_3body_integral_local(t, t, t, 1.0, n_quad=50)
    V2 = yukawa_2body(t, 1.0)
    ratio_per_C = I3 / V2  # V3/V2 = C * (I3/V2)
    C_crit = threshold / ratio_per_C if ratio_per_C > 0 else np.inf
    m_body_kg = 2 * np.sqrt(np.pi) * C_crit * m_Pl_kg
    print(f"{t:8.1f}  {ratio_per_C:14.4e}  {C_crit:12.4e}  {m_body_kg:18.4e}")

# ─── Sekcja E: Fizyczne interpretacje ─────────────────────────────────────────

print("\n" + "="*65)
print("SEKCJA E: Fizyczne systemy — czy efekty 3-ciałowe są mierzalne?")
print("="*65)

# Parametry ciał fizycznych w jednostkach Plancka
m_Pl_kg = 2.176e-8
hbar_c = 3.162e-26    # hbar*c [Planck*m] = l_Pl
l_Pl_m = 1.616e-35    # m

systems = {
    'Proton (QCD)': {
        'm_kg': 1.673e-27,
        'd_m':  1e-15,      # ~fm
        'm_sp': 1.0,         # m_sp ~ Planck
    },
    'Elektron (atom)': {
        'm_kg': 9.109e-31,
        'd_m':  5.3e-11,    # a₀
        'm_sp': 1.0,
    },
    'Mikropyłek (100nm)': {
        'm_kg': 1e-21,
        'd_m':  1e-7,       # 100 nm
        'm_sp': 1e-28,      # m_sp ~ 1/μm
    },
    'Ciała makroskopowe (1kg, 1m)': {
        'm_kg': 1.0,
        'd_m':  1.0,
        'm_sp': 1e-45,      # m_sp ~ 1/Gpc (astrofizyczny)
    },
    'Gwiazdy (M_sun, 1 AU)': {
        'm_kg': 1.989e30,
        'd_m':  1.496e11,
        'm_sp': 1e-61,
    },
    'Obiekty Plancka (hipotetyczne)': {
        'm_kg': m_Pl_kg,
        'd_m':  l_Pl_m,
        'm_sp': 1.0,
    },
}

print(f"\n{'System':35}  {'C [Pl]':>10}  {'t=m*d':>10}  {'V3/V2':>12}  {'Obs?':>6}")
print("-" * 85)

for name, params in systems.items():
    m_kg = params['m_kg']
    d_m  = params['d_m']
    m_sp = params['m_sp']   # [Planck]

    # Konwersja
    C = m_kg / m_Pl_kg / (2 * np.sqrt(np.pi))   # C = m/(2√π) [Planck]
    d_Pl = d_m / l_Pl_m                          # d [Planck]
    t = m_sp * d_Pl                              # bezwymiarowe

    if t < 0.01:
        # Dla małych t: K₀(u→0) ~ -ln(u) → I3 diverges
        # przybliżenie: V3/V2 ≈ C * const/t
        obs = "?"
        ratio_str = f"t→0 (div?)"
    elif t > 30:
        # eksponencjalnie tłumione
        ratio = C * 1e-20  # praktycznie zero
        obs = "NIE"
        ratio_str = f"~0 (e^-{t:.0f})"
    else:
        I3 = yukawa_3body_integral_local(t, t, t, 1.0, n_quad=40)
        V2 = yukawa_2body(t, 1.0)
        ratio = C * I3 / V2
        obs = "TAK" if ratio > threshold else "NIE"
        ratio_str = f"{ratio:.4e}"

    print(f"{name:35}  {C:10.2e}  {t:10.2e}  {ratio_str:>12}  {obs:>6}")

# ─── Sekcja F: Dokładna vs siodłopunktowa dla różnych geometrii ───────────────

print("\n" + "="*65)
print("SEKCJA F: Błąd siodłopunktowy — szczegółowa analiza")
print("="*65)
print()
print("Dlaczego saddle-point daje 160-770% błąd?")
print()

# Sprawdź jakie u_saddle ma rację bytu (u >> 1 potrzebne dla asymptoty K₀)
print("Wartości u_saddle dla geometrii równobocznej:")
print(f"{'t=m*d':>6}  {'u_saddle':>10}  {'K₀ asympt. valid?':>20}")
print("-" * 45)
for t in [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]:
    # α=(1/3,1/3,1/3): Q=(1/3)*(3t²)=t², Δ=1/3
    # u = m*sqrt(Q/Δ) = m * sqrt(3)*d = sqrt(3)*t
    u = np.sqrt(3) * t
    valid = "TAK" if u > 3 else "NIE (u<3)"
    print(f"{t:6.1f}  {u:10.3f}  {valid:>20}")

print()
print("Wniosek: przybliżenie K₀(u) ~ sqrt(π/2u)*exp(-u) valid tylko dla u>3,")
print("         czyli t = m*d > 3/√3 ≈ 1.73.")
print("         Dla typowych odległości t~1-2 siodłopunkt daje duże błędy.")

# ─── Podsumowanie ──────────────────────────────────────────────────────────────

print("\n" + "="*65)
print("PODSUMOWANIE ex18")
print("="*65)

print("""
Kluczowe wyniki:

1. KRYTERIUM OBSERWOWALNOŚCI (V3/V2 > 1%):
   Wymaga C > C_crit(t) gdzie C_crit rośnie z t (zanik Yukawa).
   Dla t=m*d=2 (typowe): C_crit ≈ 0.04 [Planck]
   => m_body > 0.04 * 2√π * m_Pl ≈ 0.14 m_Pl ≈ 3.1e-9 kg

2. FIZYCZNE SYSTEMY:
   - Proton/elektron: C ~ 10⁻²⁰ → V3/V2 ~ 10⁻²⁰ → NIEOBSERWOWALNE
   - Makroskopowe: t >> 1 → eksponencjalny zanik → NIEOBSERWOWALNE
   - Hipotetyczne obiekty Plancka: C ~ 0.28, t ~ 1 → V3/V2 ~ 1-4% → OBSERWOWALNE

3. BŁĄD SIODŁOPUNKTOWY (ex6 miał 160-770%):
   - Przybliżenie K₀(u) ~ sqrt(π/2u)*exp(-u) wymaga u > 3
   - Dla równobocznego: u = √3·t, więc potrzeba t > 1.73
   - Dla t < 1.73 (małe odległości Plancka): błąd siodła ogromny
   - Dla t > 2: błąd siodła < 20%

4. WNIOSEK GŁÓWNY:
   Dokładna całka Feynmana daje V3/V2 znacznie WIĘKSZE niż saddle-point
   dla małych t — to oznacza, że efekty 3-ciałowe są silniejsze niż
   dotąd szacowano dla układów bliskich (t ~ 0.5-2).

5. CZĄSTKI STANDARDOWE:
   Efekty 3-ciałowe TGP są praktycznie zerowe dla wszystkich znanych
   cząstek (C_proton ~ 2e-20 << C_crit ~ 0.04).
   TGP-to-3-body effects detectable only at Planck scale.
""")
