#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
"""
p43_gmo_koide.py
=================
Analiza P43: Relacja Gell-Mann–Okubo a warunek Koide Q=3/2

Pytanie P43 (z P42): Czy relacja GMO implikuje Q≈3/2 lub odwrotnie?
Czy perturbacja mas baryonów w duchu GMO przesuwa Q → 3/2?

Zakres:
  A. Weryfikacja GMO dla oktetu i dekupletu baryonów (PDG 2024)
  B. Jaki Q dają progresje: arytmetyczna, geometryczna, harmoniczna?
     → Przy jakiej progresji Q=3/2?
  C. Przecięcie GMO ∩ Koide — wyprowadzenie warunku analitycznego
  D. Perturbacja mas baryonów wzdłuż kierunku GMO → jak zmienia się Q?
  E. Poszukiwanie tripletów spełniających jednocześnie GMO i Q≈3/2
  F. TGP interpretacja: czy „rodzina Koide" jest szczególnym przypadkiem GMO?
"""

import numpy as np
from scipy.optimize import brentq, minimize

# =========================================================================
# PDG masy baryonów [MeV] 2024
# =========================================================================
PDG = {
    # Oktet JP=1/2+
    'p':      938.272,
    'n':      939.565,
    'Lambda': 1115.683,
    'Sigma+': 1189.37,
    'Sigma0': 1192.642,
    'Sigma-': 1197.449,
    'Xi0':    1314.86,
    'Xi-':    1321.71,
    # Dekuplet JP=3/2+
    'Delta':       1232.0,
    'Sigma_s+':    1382.80,
    'Sigma_s0':    1383.7,
    'Sigma_s-':    1387.2,
    'Xi_s0':       1531.80,
    'Xi_s-':       1535.0,
    'Omega':       1672.45,
    # Leptony (referencja)
    'e':    0.51099895,
    'mu':   105.6583755,
    'tau':  1776.86,
}

def Q(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s * s / (m1 + m2 + m3)

def dQ_pct(m1, m2, m3):
    return (Q(m1, m2, m3) - 1.5) / 1.5 * 100.0

def star(q):
    d = abs(q - 1.5)
    if d < 0.005: return " ★★★ KOIDE"
    if d < 0.02:  return " ★★  (~1%)"
    if d < 0.05:  return " ★   (~3%)"
    return ""

# =========================================================================
# SEKCJA A: Weryfikacja relacji Gell-Mann–Okubo
# =========================================================================
def section_A():
    print("\n" + "=" * 72)
    print("SEKCJA A: WERYFIKACJA RELACJI GELL-MANN–OKUBO")
    print("=" * 72)

    p = PDG
    MN   = (p['p'] + p['n']) / 2.0      # izospin avg
    ML   = p['Lambda']
    MS   = (p['Sigma+'] + p['Sigma0'] + p['Sigma-']) / 3.0
    MX   = (p['Xi0'] + p['Xi-']) / 2.0

    print(f"\n  A1: Oktet JP=1/2+  (izospinowe uśrednienie)")
    print(f"  M_N  = {MN:.3f} MeV  (p,n avg)")
    print(f"  M_Λ  = {ML:.3f} MeV")
    print(f"  M_Σ  = {MS:.3f} MeV  (Σ+,Σ⁰,Σ- avg)")
    print(f"  M_Ξ  = {MX:.3f} MeV  (Ξ⁰,Ξ- avg)")

    # Relacja GMO: 2(M_N + M_Ξ) = 3M_Λ + M_Σ
    lhs = 2.0 * (MN + MX)
    rhs = 3.0 * ML + MS
    print(f"\n  Relacja GMO: 2(M_N+M_Ξ) = 3M_Λ+M_Σ")
    print(f"    LHS = 2·({MN:.1f}+{MX:.1f}) = {lhs:.2f}")
    print(f"    RHS = 3·{ML:.1f}+{MS:.1f} = {rhs:.2f}")
    print(f"    Odchylenie = {(lhs-rhs)/rhs*100:+.3f}%")

    # Relacja GMO dla dekupletu: równorozmieszczone masy
    MD  = (p['Delta'])
    MSs = (p['Sigma_s+'] + p['Sigma_s0'] + p['Sigma_s-']) / 3.0
    MXs = (p['Xi_s0'] + p['Xi_s-']) / 2.0
    MO  = p['Omega']

    print(f"\n  A2: Dekuplet JP=3/2+  (równorozmieszczenie mas)")
    print(f"  M_Δ  = {MD:.2f} MeV")
    print(f"  M_Σ* = {MSs:.2f} MeV")
    print(f"  M_Ξ* = {MXs:.2f} MeV")
    print(f"  M_Ω  = {MO:.2f} MeV")

    # Arytmetyczna progresja: M(Δ), M(Σ*), M(Ξ*), M(Ω)
    d1 = MSs - MD
    d2 = MXs - MSs
    d3 = MO  - MXs
    print(f"\n  Różnice kolejnych mas: Δd₁={d1:.2f}, Δd₂={d2:.2f}, Δd₃={d3:.2f} MeV")
    print(f"  Progresja arytmetyczna? → d̄={np.mean([d1,d2,d3]):.2f}, std={np.std([d1,d2,d3]):.2f} MeV")
    print(f"  → {'TAK ≈ arytmetyczna' if np.std([d1,d2,d3])/np.mean([d1,d2,d3]) < 0.10 else 'NIE (std/mean > 10%)'}")

    # Q dla różnych tripletów z dekupletu
    print(f"\n  Q dla tripletów dekupletu:")
    for nm, a, b, c in [("Δ/Σ*/Ξ*", MD, MSs, MXs),
                         ("Σ*/Ξ*/Ω",  MSs, MXs, MO),
                         ("Δ/Σ*/Ω",   MD,  MSs, MO),
                         ("Δ/Ξ*/Ω",   MD,  MXs, MO)]:
        qq = Q(a, b, c)
        print(f"    {nm:<16} Q={qq:.5f}  δQ={dQ_pct(a,b,c):+.2f}%{star(qq)}")

    # Fit GMO do oktetu
    print(f"\n  A3: Parametry GMO dla oktetu (fit)")
    # M(Y, I) = M0 + a*Y + b*(I(I+1) - Y²/4)
    # N: Y=1, I=1/2  → M0 + a + b/2
    # Λ: Y=0, I=0    → M0
    # Σ: Y=0, I=1    → M0 + 2b
    # Ξ: Y=-1, I=1/2 → M0 - a + b/2
    M0   = ML
    a_gmo = (MN - MX) / 2.0
    b_gmo = (MN + MX) / 2.0 - ML
    print(f"    M₀ = {M0:.2f} MeV,  a_GMO = {a_gmo:.2f} MeV,  b_GMO = {b_gmo:.2f} MeV")
    MN_fit = M0 + a_gmo + b_gmo/2
    ML_fit = M0
    MS_fit = M0 + 2*b_gmo
    MX_fit = M0 - a_gmo + b_gmo/2
    print(f"    Fit: N={MN_fit:.1f}  Λ={ML_fit:.1f}  Σ={MS_fit:.1f}  Ξ={MX_fit:.1f}")
    print(f"    PDG: N={MN:.1f}    Λ={ML:.1f}  Σ={MS:.1f}  Ξ={MX:.1f}")


# =========================================================================
# SEKCJA B: Jaki Q dają różne typy progresji?
# =========================================================================
def section_B():
    print("\n" + "=" * 72)
    print("SEKCJA B: Q DLA RÓŻNYCH TYPÓW PROGRESJI MAS")
    print("=" * 72)
    print("""
  Progresja arytmetyczna: m_k = m₀ + k·d  (k=0,1,2)
  Progresja geometryczna: m_k = m₀·r^k    (k=0,1,2)
  Progresja harmoniczna:  1/m_k = 1/m₀ + k·δ
  Ogólna: m_k = m₀ + a·k + b·k²  (kwadratowa)
    """)

    # --- Progresja arytmetyczna ---
    print("  B1: Progresja arytmetyczna  m = (m₀, m₀+d, m₀+2d)")
    print(f"  {'d/m₀':>8}  {'Q':>10}  {'δQ':>8}")
    print(f"  {'-'*35}")
    for ratio in [0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
        m0, d = 1.0, ratio
        qq = Q(m0, m0+d, m0+2*d)
        print(f"  {ratio:>8.3f}  Q={qq:.6f}  δQ={dQ_pct(m0, m0+d, m0+2*d):+.2f}%{star(qq)}")

    # Analityczne minimum Q dla progresji arytmetycznej przy t→∞:
    # Q(1, 1+t, 1+2t) → (√t + √(2t))²/(3t) = (1+√2)²/3 = (3+2√2)/3
    Q_min_arith = (3.0 + 2.0*np.sqrt(2.0)) / 3.0
    print(f"\n  ★ TWIERDZENIE (analityczne):")
    print(f"    Dla progresji arytmetycznej m=(M, M+d, M+2d):")
    print(f"    Q ∈ ({Q_min_arith:.6f}, 3.000000]  dla d/M ∈ [0, ∞)")
    print(f"    Q_min = (3+2√2)/3 ≈ {Q_min_arith:.6f}  przy d/M → ∞")
    print(f"    Q=3/2 = 1.500000 < Q_min ≈ {Q_min_arith:.4f}")
    print(f"    → Q=3/2 jest NIEMOŻLIWE dla ŻADNEJ progresji arytmetycznej!")
    print(f"    → GMO (liniowa relacja masowa) jest FUNDAMENTALNIE NIEZGODNA z Koide Q=3/2.")

    # --- Progresja geometryczna ---
    print(f"\n  B2: Progresja geometryczna  m = (m₀, m₀·r, m₀·r²)")
    print(f"  {'r':>8}  {'Q':>10}  {'δQ':>8}")
    print(f"  {'-'*35}")
    for r in [1.1, 2.0, 5.0, 10.0, 14.4, 20.0, 50.0, 100.0, 206.8]:
        qq = Q(1.0, r, r*r)
        print(f"  {r:>8.2f}  Q={qq:.6f}  δQ={dQ_pct(1, r, r*r):+.2f}%{star(qq)}")

    # Q(1, r, r²) = (1+√r+r)²/(1+r+r²) — analitycznie dla Q=3/2
    def Q_geom_m1(r_val):
        return Q(1.0, r_val, r_val**2) - 1.5
    try:
        r_koide_g = brentq(Q_geom_m1, 1.001, 1000.0)
        print(f"\n  ★ Progresja geometryczna Q(1,r,r²)=3/2: r = {r_koide_g:.6f}")
        print(f"    (r=√(m_tau/m_e)? → {np.sqrt(PDG['tau']/PDG['e']):.4f})")
        print(f"    Sprawdzenie: Q(1,{r_koide_g:.4f},{r_koide_g**2:.4f}) = {Q(1,r_koide_g,r_koide_g**2):.6f}")
    except Exception as e:
        print(f"  Brak zera: {e}")

    # --- Leptony — jaki typ progresji? ---
    print(f"\n  B3: Leptony — diagnoza typu progresji")
    me, mm, mt = PDG['e'], PDG['mu'], PDG['tau']
    print(f"  Masy: {me:.5f}, {mm:.4f}, {mt:.2f} MeV")
    print(f"  Różnice: {mm-me:.4f}, {mt-mm:.2f} (arytmetyczna? NIE)")
    print(f"  Ilorazy: {mm/me:.2f}, {mt/mm:.2f} (geometryczna? NIE)")
    # Sprawdzamy r₁₂ = √(m₂/m₁), r₂₃ = √(m₃/m₂)
    r12 = np.sqrt(mm/me)
    r23 = np.sqrt(mt/mm)
    print(f"  √(m₂/m₁) = {r12:.4f}, √(m₃/m₂) = {r23:.4f} (√-geometryczna? {'blisko' if abs(r12-r23)/r12 < 0.05 else 'NIE'})")
    # Sprawdzamy cbrt
    r12c = (mm/me)**(1/3)
    r23c = (mt/mm)**(1/3)
    print(f"  ∛(m₂/m₁) = {r12c:.4f}, ∛(m₃/m₂) = {r23c:.4f} (∛-geometryczna? {'blisko' if abs(r12c-r23c)/r12c < 0.05 else 'NIE'})")


# =========================================================================
# SEKCJA C: Przecięcie GMO ∩ Koide — warunek analityczny
# =========================================================================
def section_C():
    print("\n" + "=" * 72)
    print("SEKCJA C: PRZECIĘCIE GMO ∩ KOIDE — WARUNEK ANALITYCZNY")
    print("=" * 72)
    print("""
  GMO (liniowy):   m_k = M₀ + k·Δ   (k=0,1,2, Δ=stała odstępu)
  Koide (Q=3/2):   (√m₀+√m₁+√m₂)² = (3/2)(m₀+m₁+m₂)

  Podstawiamy: m_k = M + k·d (M = M₀, d = Δ)
  Szukamy: przy jakim stosunku x = d/M progresja arytmetyczna ma Q=3/2?
    """)

    # Q(M, M+d, M+2d) = 3/2
    # (√M + √(M+d) + √(M+2d))² = (3/2)(3M + 3d) = (9/2)(M+d)
    # Podstawiamy t = d/M (względny odstęp)
    # (1 + √(1+t) + √(1+2t))² = (3/2)(3 + 3t) = (9/2)(1+t)
    # Szukamy t numerycznie

    def Q_arith_centered(t):
        """Q dla m=(1, 1+t, 1+2t)."""
        return Q(1.0, 1.0+t, 1.0+2.0*t) - 1.5

    # Analiza: Q maleje z t (od 3 przy t=0 do 1 przy t→∞)
    t_vals = np.logspace(-3, 3, 1000)
    Q_vals = np.array([Q(1, 1+t, 1+2*t) for t in t_vals])

    # Analityczne: Q_min dla progresji arytmetycznej = (3+2√2)/3 ≈ 1.943
    Q_min_arith = (3.0 + 2.0*np.sqrt(2.0)) / 3.0

    # Sprawdzenie wartości Q dla dużych t
    print(f"  Wartości Q dla dużych t = d/M₀:")
    for t_test in [10, 100, 1000, 10000, 1e6]:
        qq = Q(1.0, 1.0+t_test, 1.0+2.0*t_test)
        print(f"    t={t_test:.0e}  Q={qq:.8f}  (granica: {Q_min_arith:.8f})")

    print(f"""
  ★★★ TWIERDZENIE GŁÓWNE P43:

    Q_min(progresja arytmetyczna) = (3+2√2)/3 ≈ {Q_min_arith:.6f} > 3/2

    Dowód:
      m = (M, M+d, M+2d); t = d/M; Q(t) = (1+√(1+t)+√(1+2t))² / (3+3t)
      lim_{{t→∞}} Q(t) = lim (√t+√(2t))²/(3t) = (1+√2)²/3 = (3+2√2)/3 ≈ 1.943
      Q(t) jest ściśle malejące: Q ∈ ((3+2√2)/3, 3)
      3/2 = 1.5 < 1.943 ≤ Q(t)  dla wszystkich t ≥ 0

    WNIOSEK: Q=3/2 (Koide) jest NIEMOŻLIWE dla JAKIEJKOLWIEK progresji arytmetycznej.
    GMO i Koide to FUNDAMENTALNIE NIEKOMPATYBILNE relacje masowe.
    """)

    return Q_min_arith  # zwracamy Q_min zamiast t_koide


# =========================================================================
# SEKCJA D: Perturbacja mas baryonów wzdłuż kierunku GMO → Q(t)
# =========================================================================
def section_D():
    print("\n" + "=" * 72)
    print("SEKCJA D: PERTURBACJA MAS BARYONÓW → ZMIANA Q(t)")
    print("=" * 72)
    print("""
  Startujemy od oktetu (p, Λ, Ξ⁰). Definiujemy „kierunek GMO":
  δM_k = k·Δ (perturbacja liniowa w indeksie strangeness).
  Pytamy: czy istnieje Δ ≠ 0 takie, że Q(p+0·Δ, Λ+1·Δ, Ξ+2·Δ) = 3/2?
    """)

    Mp, ML, MX = PDG['p'], PDG['Lambda'], PDG['Xi0']

    print(f"  Punkt startowy: p={Mp:.1f}, Λ={ML:.1f}, Ξ⁰={MX:.1f} MeV")
    print(f"  Q₀ = {Q(Mp, ML, MX):.5f}  δQ₀ = {dQ_pct(Mp, ML, MX):+.2f}%")

    # Q(Δ) = Q(Mp, ML+Δ, MX+2Δ)
    def Q_perturb(Delta):
        return Q(Mp, ML + Delta, MX + 2.0*Delta)

    Delta_vals = np.linspace(-500, 500, 1001)
    Q_vals = np.array([Q_perturb(d) for d in Delta_vals])

    print(f"\n  Q(Δ) dla perturbacji GMO (p stały, Λ→Λ+Δ, Ξ→Ξ+2Δ):")
    print(f"  {'Δ [MeV]':>10}  {'Q':>10}  {'δQ':>8}")
    print(f"  {'-'*35}")
    for dv in [-400, -300, -200, -100, 0, 100, 200, 300, 400]:
        qq = Q_perturb(dv)
        print(f"  {dv:>10.0f}  Q={qq:.5f}  δQ={dQ_pct(Mp, ML+dv, MX+2*dv):+.2f}%{star(qq)}")

    # Szukamy Δ dającego Q=3/2
    # Q maleje z Δ (dla Δ>0 masy rosną → Q maleje; dla Δ<0 masy maleją → Q rośnie)
    min_Q_Delta = np.min(Q_vals)
    max_Q_Delta = np.max(Q_vals)
    print(f"\n  Q_min = {min_Q_Delta:.4f}, Q_max = {max_Q_Delta:.4f}")

    if min_Q_Delta < 1.5 < max_Q_Delta:
        try:
            Delta_koide = brentq(lambda d: Q_perturb(d) - 1.5,
                                 Delta_vals[0], Delta_vals[-1], xtol=1e-8)
            print(f"\n  ★ Q=3/2 osiągane przy Δ = {Delta_koide:.2f} MeV")
            Ml_new = ML + Delta_koide
            Mx_new = MX + 2 * Delta_koide
            print(f"    Masy zmodyfikowane: p={Mp:.1f}, Λ={Ml_new:.1f}, Ξ={Mx_new:.1f} MeV")
            print(f"    ΔM_Λ = {Delta_koide:+.1f} MeV  ({Delta_koide/ML*100:+.1f}%)")
            print(f"    ΔM_Ξ = {2*Delta_koide:+.1f} MeV  ({2*Delta_koide/MX*100:+.1f}%)")
            print(f"    Weryfikacja: Q = {Q(Mp, Ml_new, Mx_new):.6f}")
        except Exception as e:
            print(f"  Nie znaleziono: {e}")
    else:
        print(f"\n  Q=3/2 NIE jest osiągane dla perturbacji GMO w zakresie Δ∈[-500, 500]")
        print(f"  Q zawsze > {min_Q_Delta:.3f} i < {max_Q_Delta:.3f}")

    # Sprawdzamy też oktet N, Σ, Ξ (bez Λ)
    print(f"\n  D2: Perturbacja GMO dla N/Σ/Ξ (z izospinowym uśrednieniem):")
    MN = (PDG['p'] + PDG['n']) / 2
    MS = (PDG['Sigma+'] + PDG['Sigma0'] + PDG['Sigma-']) / 3
    MX2 = (PDG['Xi0'] + PDG['Xi-']) / 2
    print(f"  N={MN:.1f}, Σ={MS:.1f}, Ξ={MX2:.1f} MeV → Q={Q(MN,MS,MX2):.5f}  δQ={dQ_pct(MN,MS,MX2):+.2f}%")

    # D3: Triplet pionowy przez pokolenia kwarkowe
    print(f"\n  D3: Triplety baryonowe 'pionowe' (przez pokolenia):")
    for nm, m1, m2, m3 in [
        ("p/Λc/Λb", PDG['p'], 2286.46, 5619.60),
        ("n/Λc/Λb", PDG['n'], 2286.46, 5619.60),
        ("Λ/Λc/Λb", PDG['Lambda'], 2286.46, 5619.60),
        ("Ξ⁰/Λc/Λb", PDG['Xi0'], 2286.46, 5619.60),
    ]:
        qq = Q(m1, m2, m3)
        print(f"    {nm:<18} Q={qq:.5f}  δQ={dQ_pct(m1,m2,m3):+.2f}%{star(qq)}")


# =========================================================================
# SEKCJA E: Szukanie tripletów GMO z Q≈3/2 — przestrzeń (M₀, a, b)
# =========================================================================
def section_E():
    print("\n" + "=" * 72)
    print("SEKCJA E: PARAMETRY GMO DAJĄCE Q=3/2 — PRZESTRZEŃ (M₀, Δ)")
    print("=" * 72)
    print("""
  Parameteryzacja: m(k) = M₀ + k·Δ  (k=0,1,2)
  Warunek Q=3/2: (√M₀ + √(M₀+Δ) + √(M₀+2Δ))² = (3/2)(3M₀+3Δ)

  Twierdzenie (P43, sekcja C): Q_min(arytm.) = (3+2√2)/3 ≈ 1.943 > 3/2
  → progresja arytmetyczna (GMO) NIGDY nie daje Q=3/2 dla żadnych M₀, Δ!
    """)

    # Z twierdzenia P43: Q_min(arytm) = (3+2√2)/3 ≈ 1.943 > 3/2
    # → progresja arytmetyczna NIE może dać Q=3/2
    Q_min_arith = (3.0 + 2.0*np.sqrt(2.0)) / 3.0

    print(f"  ★ Twierdzenie P43: Q_min(arytm.) = (3+2√2)/3 = {Q_min_arith:.6f} > 3/2")
    print(f"  → Progresja arytmetyczna (GMO) NIE może dać Q=3/2 DLA ŻADNEGO M₀, Δ!")
    print(f"\n  Zamiast tego: geometryczna Q=3/2 przy r = 22.956")
    r_geom = 22.956439
    print(f"  Dla geometrycznej: masy ∝ (1, r, r²) = (1, {r_geom:.2f}, {r_geom**2:.1f})")

    # Pokaż co daje GMO dla faktycznych mas baryonowych vs. Q=3/2
    print(f"\n  Porównanie GMO dla baryonów vs. q=3/2 (geometryczna):")
    print(f"  {'Typ mas':<35} {'Q':>10}")
    print(f"  {'-'*48}")
    for M0 in [938.0, 1115.0, 1232.0]:
        for t_ratio in [0.2, 0.4, 0.693]:
            Delta = t_ratio * M0
            m0, m1, m2 = M0, M0+Delta, M0+2*Delta
            lbl = f"arytm. M₀={M0:.0f}, d/M₀={t_ratio:.2f}"
            print(f"  {lbl:<35} Q={Q(m0,m1,m2):.5f}")

    # Faktyczny przykład: czy jakakolwiek seria baryonów bliska tej krzywej?
    print(f"\n  Odniesienie do faktycznych baryonów PDG:")
    real_series = [
        ("N, Σ, Ξ (oktet)",      (PDG['p']+PDG['n'])/2,
                                  (PDG['Sigma+']+PDG['Sigma0']+PDG['Sigma-'])/3,
                                  (PDG['Xi0']+PDG['Xi-'])/2),
        ("Δ, Σ*, Ξ* (dekuplet)", PDG['Delta'],
                                  (PDG['Sigma_s+']+PDG['Sigma_s0']+PDG['Sigma_s-'])/3,
                                  (PDG['Xi_s0']+PDG['Xi_s-'])/2),
    ]
    for nm, m0, m1, m2 in real_series:
        Delta_real = (m2 - m0) / 2.0
        t_real = Delta_real / m0
        qq = Q(m0, m1, m2)
        Q_min_str = f"Q_min(arytm)={(3+2*np.sqrt(2))/3:.3f}"
        print(f"    {nm:<30} Δ/M₀={t_real:.3f}  ({Q_min_str})")
        print(f"    {'':30} Q={qq:.5f}  δQ={dQ_pct(m0,m1,m2):+.2f}%")


# =========================================================================
# SEKCJA F: TGP interpretacja — co wyróżnia leptony?
# =========================================================================
def section_F():
    print("\n" + "=" * 72)
    print("SEKCJA F: INTERPRETACJA TGP — CO WYRÓŻNIA LEPTONY?")
    print("=" * 72)

    # Charakterystyczny iloraz mas letonowych
    me, mm, mt = PDG['e'], PDG['mu'], PDG['tau']
    r21_lep = mm / me
    r31_lep = mt / me

    # Leptony: r₂₁≈207, r₃₁≈3477
    # Kwarki MS-bar: r₂₁(d/s/b) ≈ 20, 895; r₂₁(u/c/t) ≈ 588, 79949
    print(f"""
  Leptony (e/μ/τ):   r₂₁ = {r21_lep:.1f},  r₃₁ = {r31_lep:.1f}
  Kwarki u/c/t:      r₂₁ ≈ 588,    r₃₁ ≈ 79949
  Kwarki d/s/b:      r₂₁ ≈ 20,     r₃₁ ≈ 895

  Zatem leptony NIE leżą na krzywej GMO (progresja arytmetyczna r₂₁≈1+t, t<<1),
  ani na krzywej geometrycznej prostej (r₁₂≠r₂₃).

  Leptony leżą na SPECYFICZNEJ krzywej Q=3/2 o parametrach:
    r₃₁ = r₃₁_Koide(r₂₁) = [2(1+√r₂₁) + √(3(1+4√r₂₁+r₂₁))]²

  TGP WYJAŚNIA dlaczego ta krzywa:
    K₃ = C·a/√λ jest universalne (niezależne od α_f)
    K₁ = K₁(α_f, a) jest specyficzne dla rodziny
    r₃₁ = K₃/K₁ → przy λ=λ_Koide, para (r₂₁, r₃₁) leży na krzywej Q=3/2

  GMO i Koide to NIEZALEŻNE relacje:
    GMO: masa ~ liniowa funkcja ładunku (strangeness, charm, beauty)
    Koide: √masa ~ liniowa relacja (fenomenologiczna, wywodzona z TGP)

  Połączenie możliwe w P44:
    Jeśli baryony spełniają GMO, a jednocześnie masa efektywna kwarka
    (konstytuent) jest solitonem TGP z Q=3/2 → baryon może się „rozpaść"
    na trzy solitony Koide → byłoby to geometryczne wyjaśnienie kompozytowości.
    """)

    # Numeryczna weryfikacja: porównaj Q leptonów z Q baryonów
    print(f"  Porównanie Q:")
    print(f"    Q(e/μ/τ)        = {Q(me, mm, mt):.6f}  ← Q=3/2 ✓ (TGP soliton)")
    MN = (PDG['p'] + PDG['n']) / 2
    ML = PDG['Lambda']
    MX = (PDG['Xi0'] + PDG['Xi-']) / 2
    MS = (PDG['Sigma+'] + PDG['Sigma0'] + PDG['Sigma-']) / 3
    print(f"    Q(N/Λ/Ξ) oktet = {Q(MN, ML, MX):.6f}  ← Q>>3/2 (bariiony)")
    print(f"    Q(N/Σ/Ξ) oktet = {Q(MN, MS, MX):.6f}  ← Q>>3/2 (bariiony)")

    # Twierdzenie P43: Q_min(arytm.) = (3+2√2)/3 ≈ 1.943 > 3/2
    # → progresja arytmetyczna NIGDY nie daje Q=3/2
    Q_min_arith = (3.0 + 2.0*np.sqrt(2.0)) / 3.0
    Delta_N_Xi = (MX - MN) / 2.0
    t_baryon = Delta_N_Xi / MN
    print(f"\n  Q_min(arytm. progresja) = (3+2√2)/3 = {Q_min_arith:.5f}  > 3/2 = 1.50000")
    print(f"  t_baryon (N/Λ/Ξ GMO)   = Δ/M_N = {t_baryon:.5f}  → Q={Q(1, 1+t_baryon, 1+2*t_baryon):.4f}")
    print(f"\n  WNIOSEK: Progresja arytmetyczna ma Q ∈ ({Q_min_arith:.3f}, 3.000) → Q=3/2 niemożliwe.")
    print(f"  Relacja Q=3/2 i GMO są FUNDAMENTALNIE NIEKOMPATYBILNE (twierdzenie analityczne).")


# =========================================================================
# SEKCJA G: Podsumowanie P43
# =========================================================================
def section_G():
    print("\n" + "=" * 72)
    print("SEKCJA G: PODSUMOWANIE P43")
    print("=" * 72)
    print("""
  PYTANIE P43: Czy relacja GMO implikuje Q≈3/2?
              Czy perturbacja GMO przesuwa Q → 3/2 dla baryonów?

  ODPOWIEDŹ:

  1. GMO i Koide Q=3/2 są NIEZALEŻNYMI relacjami:
     • GMO: masa ~ liniowa funkcja quantum number (strangeness)
     • Koide: wymaga SPECYFICZNYCH ilorazów mas (r₂₁≈207, r₃₁≈3477 dla leptonów)
     • ★ TWIERDZENIE: Q_min(progresja arytmetyczna) = (3+2√2)/3 ≈ 1.943 > 3/2
       → Q=3/2 jest NIEMOŻLIWE dla ŻADNEJ progresji arytmetycznej (dowód analityczny)
       → GMO (liniowe masy) i Koide (Q=3/2) są fundamentalnie niekompatybilne

  2. Perturbacja GMO dla baryonów (p/Λ/Ξ) nie osiąga Q=3/2:
     Q ∈ (2.87, 3.00) dla perturbacji Δ∈(-500,+500) MeV — nigdy blisko 3/2.

  3. Bariiony ze spełnioną GMO mają Q ≈ 2.97–3.00 (blisko górnej granicy),
     bo ich masy są prawie równe w skali bezwzględnej.

  WNIOSEK GŁÓWNY P43:
     GMO → Q>>3/2  (bariiony)
     Koide Q=3/2  → wymagane r₂₁, r₃₁ charakterystyczne dla leptonów
     ŻADNA liniowa relacja masowa (GMO-type) nie daje Q=3/2 dla typowych hadronów.

  INTERPRETACJA TGP:
     Koide Q=3/2 to właściwość SOLITONU TGP, nie układu złożonego.
     Leptony: fundamentalne solitony (wolne cząstki) → Q=3/2 z definicji λ_Koide
     Hadrony: złożone z kwarków (sfrustrowanych solitonów) → GMO, Q>>3/2

  ZAMKNIĘCIE CYKLU P38–P43:
     P38: Q=3/2 tylko leptony (wstępna diagnoza)
     P39: d/s/b osiągalne w TGP przy α_f=0.22
     P40: K₃ universalne; Q_TGP≈1.50–1.53 dla WSZYSTKICH rodzin
     P41: hadrony NIE mają Q≈3/2 (brute-force 15k kombinacji)
     P42: masy konstituentów NIE dają Q≈3/2; diagram (r₂₁,r₃₁) ustalony
     P43: GMO i Koide NIEZALEŻNE; liniowa progresja masowa NIE prowadzi do Q=3/2

  NOWY KIERUNEK P44:
     Analityczne wyprowadzenie DLACZEGO Q=3/2 jest charakterystyczne wyłącznie
     dla solitonu TGP z λ=λ_Koide — zakończenie cyklu analitycznego.
     Ewentualnie: rozszerzenie TGP o kolor (wielociałowe solitony).
    """)


# =========================================================================
# MAIN
# =========================================================================
if __name__ == '__main__':
    print("=" * 72)
    print("P43: RELACJA GELL-MANN–OKUBO A WARUNEK KOIDE Q=3/2")
    print("Czy GMO implikuje Q≈3/2? Czy baryony mogą mieć Q=3/2?")
    print("=" * 72)

    section_A()
    section_B()
    t_k = section_C()
    section_D()
    section_E()
    section_F()
    section_G()

    print("\n" + "=" * 72)
    print("KONIEC ANALIZY P43")
    print("=" * 72)
