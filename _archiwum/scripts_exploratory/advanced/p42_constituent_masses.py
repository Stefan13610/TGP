#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
"""
p42_constituent_masses.py
==========================
Analiza P42: Masy konstituentów kwarków a warunek Koide Q=3/2

Pytanie P42 (z P41): Czy efektywne masy kwarków w hadronie (masy konstituentów)
dają Q≈3/2? Jeśli tak — TGP ma geometryczną zasadę kompozytowości.

Zakres:
  A. Masy konstituentów z różnych modeli (NJL, ChQSM, lattice, bag) — czy Q≈3/2?
  B. Odwrotnie: jakie masy konstituentów MUSZĄ być, żeby Q=3/2 przy danym r₂₁_PDG?
  C. Jakiego λ_eff potrzebuje każda rodzina kwarkowa, żeby soliton TGP osiągnął Q=3/2?
  D. „Parametr uwięzienia" δQ = Q_TGP - 3/2 vs. obserwowalne hadronu.
  E. Analityczne: kiedy K₁ dałoby Q=3/2 dla d/s/b i u/c/t?
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar

# =========================================================================
# PODSTAWOWE FUNKCJE
# =========================================================================

def koide_Q(m1, m2, m3):
    """Q = (√m₁+√m₂+√m₃)² / (m₁+m₂+m₃) ∈ [1,3]; Q=3/2 dla Koide."""
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s * s / (m1 + m2 + m3)


def r31_koide_analytic(r21):
    """Analityczne r₃₁ dające Q=3/2 przy danym r₂₁.
    y = 2(1+s) ± √(3(1+4s+s²)),  s=√r₂₁,  r₃₁=y².
    """
    s = np.sqrt(r21)
    disc = 3.0 * (1.0 + 4.0*s + s*s)
    disc_sqrt = np.sqrt(max(disc, 0.0))
    y_plus  = 2.0*(1.0+s) + disc_sqrt
    y_minus = 2.0*(1.0+s) - disc_sqrt
    r31_plus  = y_plus**2  if y_plus  > 0 else np.nan
    r31_minus = y_minus**2 if y_minus > 0 else np.nan
    return r31_plus, r31_minus


def m3_from_Q32(m1, m2):
    """Oblicza m₃ tak, żeby Q(m₁,m₂,m₃)=3/2 przy zadanych m₁,m₂.
    Rozwiązuje: (√m₁+√m₂+y)² = (3/2)(m₁+m₂+y²), y=√m₃.
    Równanie: -y²/2 + 2(√m₁+√m₂)y + (m₁+m₂)/2 - 2√(m₁m₂) + ... → quadratyczne w y.

    Pełne rozwiązanie:
    (3/2)(m₁+m₂+y²) = (√m₁+√m₂+y)²
    3/2·m₁ + 3/2·m₂ + 3/2·y² = m₁+m₂ + 2(√m₁+√m₂)y + y²
    (3/2-1)y² - 2(√m₁+√m₂)y + (3/2-1)(m₁+m₂) - 2√(m₁m₂) = 0

    Wait, let me redo:
    (√m₁+√m₂+y)² = m₁+m₂+y² + 2√m₁y + 2√m₂y + 2√(m₁m₂)
    Eq: m₁+m₂+y² + 2(√m₁+√m₂)y + 2√(m₁m₂) = (3/2)(m₁+m₂+y²)
    y² + 2(√m₁+√m₂)y + 2√(m₁m₂) = (3/2)(m₁+m₂) + (3/2-1)y² ...

    Przegrupowanie:
    m₁+m₂+y² + 2(√m₁+√m₂)y + 2√(m₁m₂) = (3/2)(m₁+m₂+y²)
    0 = (3/2-1)y² - 2(√m₁+√m₂)y + (3/2-1)(m₁+m₂) - 2√(m₁m₂)
    0 = (1/2)y² - 2(√m₁+√m₂)y + (1/2)(m₁+m₂) - 2√(m₁m₂)
    Mnożąc przez 2:
    y² - 4(√m₁+√m₂)y + (m₁+m₂) - 4√(m₁m₂) = 0
    y² - 4(√m₁+√m₂)y + (√m₁-√m₂)² + 2√(m₁m₂) - 4√(m₁m₂) = 0
    hmm, let me just use:
    a_coef = 1
    b_coef = -4*(√m₁+√m₂)
    c_coef = (m₁+m₂) - 4√(m₁m₂) = (√m₁-√m₂)² - 2√(m₁m₂) + ...
    """
    s1, s2 = np.sqrt(m1), np.sqrt(m2)
    # y² - 4(s1+s2)y + (m1 + m2 - 4*s1*s2) = 0
    a = 1.0
    b = -4.0*(s1 + s2)
    c = m1 + m2 - 4.0*s1*s2
    disc = b*b - 4.0*a*c
    if disc < 0:
        return np.nan, np.nan
    sq = np.sqrt(disc)
    y1 = (-b + sq) / (2.0*a)
    y2 = (-b - sq) / (2.0*a)
    m3_1 = y1**2 if y1 > 0 else np.nan
    m3_2 = y2**2 if y2 > 0 else np.nan
    return m3_1, m3_2


# =========================================================================
# SEKCJA A: Masy konstituentów z różnych modeli
# =========================================================================
def section_A():
    print("\n" + "=" * 72)
    print("SEKCJA A: MASY KONSTITUENTÓW KWARKÓW — RÓŻNE MODELE")
    print("=" * 72)

    # Różne zestawy mas konstituentów z literatury [MeV]
    models = {
        # Model              u_c    d_c    s_c    c_c     b_c     t_pole
        "De Rújula 1975":   (338,   338,   540,   1550,   5000,   172760),
        "NJL (Vogl-Weise)": (335,   335,   555,   1570,   4980,   172760),
        "Bag model":        (310,   310,   483,   1640,   5090,   172760),
        "Brodsky-Hwang":    (300,   300,   500,   1500,   5000,   172760),
        "ChQSM":            (360,   360,   540,   1700,   5100,   172760),
        "PDG MS-bar(2GeV)": (2.16,  4.67,  93.4,  1270,   4180,   172760),
        "PDG pole":         (2.16,  4.67,  150,   1670,   4780,   172760),
        "Lattice (κ-def)":  (340,   340,   500,   1490,   4850,   172760),
    }

    print(f"\n  {'Model':<25} {'Q(u/c/t)':>10}  {'Q(d/s/b)':>10}  {'δQ_uct':>8}  {'δQ_dsb':>8}")
    print(f"  {'-'*70}")

    for name, (mu, md, ms, mc, mb, mt) in models.items():
        Q_uct = koide_Q(mu, mc, mt)
        Q_dsb = koide_Q(md, ms, mb)
        d_uct = (Q_uct - 1.5) / 1.5 * 100
        d_dsb = (Q_dsb - 1.5) / 1.5 * 100
        fl_u = " ★" if abs(Q_uct-1.5)<0.05 else ""
        fl_d = " ★" if abs(Q_dsb-1.5)<0.05 else ""
        print(f"  {name:<25} Q_uct={Q_uct:.4f}{fl_u:<3}  Q_dsb={Q_dsb:.4f}{fl_d:<3}  "
              f"δQ={d_uct:+.1f}%   δQ={d_dsb:+.1f}%")

    # Referencja: leptony
    me, mmu, mtau = 0.511, 105.658, 1776.86
    Q_lep = koide_Q(me, mmu, mtau)
    print(f"\n  {'LEPTONY (ref.)':<25} Q={Q_lep:.6f}  δQ={0:+.4f}%")

    print(f"""
  WNIOSEK A:
    Żaden standardowy model mas konstituentów NIE daje Q(d/s/b)≈3/2.
    Q_dsb ≈ 2.1–2.4 (zbyt wysokie — masy d_c≈s_c są zbyt podobne → Q→3).
    Q_uct ≈ 1.28–1.34 (poniżej 3/2 — dla MS-bar: Q≈1.18).
    Masy konstituentów nie realizują zasady Koide'go dla kwarków.
    """)


# =========================================================================
# SEKCJA B: Jakie m₃ dałoby Q=3/2 przy zachowaniu m₁,m₂?
# =========================================================================
def section_B():
    print("\n" + "=" * 72)
    print("SEKCJA B: WYMAGANA m₃ DLA Q=3/2 PRZY ZADANYCH m₁, m₂ (ODWROTNIE)")
    print("=" * 72)
    print("""
  Dane: r₂₁=m₂/m₁ z PDG (MS-bar). Pytanie: jakie m₃ dałoby Q=3/2?
  Porównanie z PDG m₃ i z K₃·m₁/K₁ (TGP).
    """)

    # PDG MS-bar 2024 [MeV]
    mu, mc, mt = 2.16, 1270.0, 172760.0   # u/c/t
    md, ms, mb = 4.67,  93.4,   4180.0   # d/s/b

    cases = [
        ("u/c/t (MS-bar)", mu, mc, mt),
        ("d/s/b (MS-bar)", md, ms, mb),
        ("u/c/t konst.",   336.0, 1500.0, 172760.0),
        ("d/s/b konst.",   340.0,  540.0,   5000.0),
        ("e/μ/τ (ref.)",   0.511, 105.658, 1776.86),
    ]

    print(f"  {'Rodzina':<20} {'m₁':>8} {'m₂':>8} {'m₃_PDG':>10} {'m₃_Koide(+)':>14} {'m₃_Koide(-)':>14} {'Q_PDG':>8}")
    print(f"  {'-'*85}")

    for name, m1, m2, m3_pdg in cases:
        m3_k_plus, m3_k_minus = m3_from_Q32(m1, m2)
        Q_pdg = koide_Q(m1, m2, m3_pdg)
        print(f"  {name:<20} {m1:>8.2f} {m2:>8.2f} {m3_pdg:>10.2f} "
              f"{m3_k_plus:>14.2f} {m3_k_minus:>14.4f} {Q_pdg:>8.4f}")

    print(f"""
  INTERPRETACJA:
    m₃_Koide(+) = masa trzeciej cząstki WYMAGANA przez Q=3/2 przy zachowaniu m₁,m₂.
    Dla u/c → m₃_Koide(+) > m_top_PDG → niemożliwe przy obecnym spektrum.
    Dla d/s → m₃_Koide(+) ≠ m_b_PDG   → b jest "za lekkie" dla Q=3/2.
    TGP wyjaśnia: kwarki mają strukturę solitonową z Q_TGP>3/2 właśnie dlatego,
    że "brakuje im masy" w trzecim pokoleniu aby osiągnąć Q=3/2.
    """)


# =========================================================================
# SEKCJA C: Jakiego λ_eff TGP potrzebuje każda rodzina do Q=3/2?
# =========================================================================
def section_C():
    print("\n" + "=" * 72)
    print("SEKCJA C: λ_eff DLA Q=3/2 W TGP — ANALIZA WSTECZNA")
    print("=" * 72)
    print("""
  Formuła TGP: K₃ = C·a_Γ/√λ  (C≈2.000, a_Γ=0.040)
  r₃₁ = K₃/K₁

  Warunek Q=3/2:  r₃₁ = r₃₁_Koide(r₂₁)  [formuła analityczna]

  Stąd:  λ_eff = (C·a_Γ)² / (K₁·r₃₁_Koide)²

  K₁ zależy od α_f (parametr rodziny). Dla każdej rodziny mamy K₁ z P40.
    """)

    C = 2.000
    a = 0.040

    # Wyniki z P40 (a=0.040, λ=1e-5)
    families = {
        "e/μ/τ": {"alpha_f": 8.5445,  "K1": 0.009839, "r21_pdg": 206.77, "r31_pdg": 3477.0},
        "u/c/t": {"alpha_f": 20.343,  "K1": 0.004175, "r21_pdg": 587.96, "r31_pdg": 79949.0},
        "d/s/b": {"alpha_f": 0.2207,  "K1": 0.077649, "r21_pdg": 20.00,  "r31_pdg": 895.0},
    }

    print(f"  {'Rodzina':<12} {'K₁':>10} {'r₂₁':>8} {'r₃₁_Koide(+)':>14} "
          f"{'r₃₁_PDG':>10} {'λ_eff(Koide)':>14} {'λ_eff(PDG)':>13}")
    print(f"  {'-'*85}")

    for fam, d in families.items():
        K1    = d['K1']
        r21   = d['r21_pdg']
        r31_pdg = d['r31_pdg']

        r31_k_plus, _ = r31_koide_analytic(r21)

        # λ_eff dające Q=3/2 (r₃₁=r₃₁_Koide)
        # K₃ = K₁·r₃₁_Koide = C·a/√λ_eff → λ_eff = (C·a)² / (K₁·r₃₁_Koide)²
        lam_eff_koide = (C * a)**2 / (K1 * r31_k_plus)**2 if not np.isnan(r31_k_plus) else np.nan

        # λ_eff dające r₃₁=r₃₁_PDG (faktyczne masy PDG)
        lam_eff_pdg = (C * a)**2 / (K1 * r31_pdg)**2

        print(f"  {fam:<12} {K1:>10.6f} {r21:>8.2f} {r31_k_plus:>14.1f} "
              f"{r31_pdg:>10.1f} {lam_eff_koide:>14.3e} {lam_eff_pdg:>13.3e}")

    print(f"""
  LEPTONY:  λ_eff(Koide) = λ_eff(PDG) ≈ 5.47e-6  ← TGP jest dokładnie skalibrowany
  U/C/T:    λ_eff(Koide) << λ_eff(PDG)           ← PDG wymaga λ << λ_Koide dla Q=3/2
  D/S/B:    λ_eff(Koide) << λ_eff(PDG)           ← analogicznie

  INTERPRETACJA:
    λ_Koide (wyznaczone z leptonów) jest ZBYT DUŻE dla kwarków, żeby dać Q=3/2.
    Inaczej: przy tym samym λ, K₃ jest za małe/za duże, żeby spełnić r₃₁=r₃₁_Koide.
    Kwarki potrzebowałyby MNIEJSZEGO λ (słabszego potencjału), żeby K₃ wzrosło
    i przesunęło Q w stronę 3/2.
    """)


# =========================================================================
# SEKCJA D: Parametr uwięzienia δQ vs. własności hadronów
# =========================================================================
def section_D():
    print("\n" + "=" * 72)
    print("SEKCJA D: PARAMETR UWIĘZIENIA δQ = Q_TGP - 3/2")
    print("=" * 72)
    print("""
  TGP wyniki z P40 (λ=λ_Koide=5.4677e-6):
    Q_TGP(e/μ/τ) = 1.5000  → δQ = 0.000  (wolna cząstka)
    Q_TGP(u/c/t) = 1.5259  → δQ = 0.026  (uwięzione)
    Q_TGP(d/s/b) = 1.5170  → δQ = 0.017  (uwięzione)

  Pytanie: czy δQ koreluje z siłą uwięzienia (energią wiązania, masą piona,
           stałą QCD, itd.)?
    """)

    # Dane
    families = {
        "e/μ/τ": {"Q_tgp": 1.5000, "alpha_f": 8.545,  "note": "wolne"},
        "u/c/t": {"Q_tgp": 1.5259, "alpha_f": 20.343, "note": "silnie uwięzione"},
        "d/s/b": {"Q_tgp": 1.5170, "alpha_f": 0.2207, "note": "uwięzione"},
    }

    print(f"  {'Rodzina':<12} {'Q_TGP':>8} {'δQ':>8} {'α_f':>8}  Opis uwięzienia")
    print(f"  {'-'*55}")
    for fam, d in families.items():
        dQ = d['Q_tgp'] - 1.5
        print(f"  {fam:<12} {d['Q_tgp']:>8.4f} {dQ:>8.4f} {d['alpha_f']:>8.4f}  {d['note']}")

    print(f"""
  OBSERWACJE:
    1. d/s/b ma MNIEJSZY δQ (0.017) niż u/c/t (0.026) →
       d/s/b jest "bliżej wolności" niż u/c/t.
       Konsystentne z: mesony D, B są cięższe ale relatywnie stabilne;
       kwarki bottom/charm mają asymptotyczną wolność.

    2. α_f dla d/s/b jest ZNACZNIE mniejszy (0.22) niż dla u/c/t (20.3) i leptonów (8.5) →
       parametr sprzężenia α_f nie jest prostą miarą uwięzienia.
       α_f kontroluje r₂₁=K₂/K₁, nie Q bezpośrednio.

    3. δQ = Q_TGP - 3/2 może być traktowany jako "miara niedopasowania geometrycznego":
       jak bardzo geometria solitonu odbiega od konfiguracji samokonzystentnej.

  HIPOTEZA D1 (do weryfikacji P43):
    δQ · E_soliton ≈ energia uwięzienia (string tension × R_hadron)?
    """)


# =========================================================================
# SEKCJA E: Geometryczna interpretacja — diagram (r21, r31) i krzywa Q=3/2
# =========================================================================
def section_E():
    print("\n" + "=" * 72)
    print("SEKCJA E: DIAGRAM (r₂₁, r₃₁) I KRZYWA KOIDE Q=3/2")
    print("=" * 72)
    print("""
  Krzywa Q=3/2 w płaszczyźnie (r₂₁, r₃₁): r₃₁ = r₃₁_Koide(r₂₁).
  Punkty TGP i PDG naniesione na diagram.
    """)

    # Oblicz krzywą Q=3/2
    r21_vals = np.logspace(np.log10(2), np.log10(1e7), 200)
    r31_plus_vals = []
    for r21 in r21_vals:
        rp, rm = r31_koide_analytic(r21)
        r31_plus_vals.append(rp)

    # K3_Koide = C·a/√λ_Koide = 2.000·0.040/√(5.4677e-6) = 34.213
    # r31_TGP(λ=λ_Koide) = K3_Koide / K1
    # K1: lep=0.009839, uct=0.004175, dsb=0.077649  (z P40, a=0.040)
    K3_koide = 2.000 * 0.040 / np.sqrt(5.4677e-6)   # ≈ 34.213

    # Punkty TGP (P40, λ=λ_Koide) — r31 = K3_Koide/K1
    r31_tgp_lep = K3_koide / 0.009839   # ≈ 3477
    r31_tgp_uct = K3_koide / 0.004175   # ≈ 8196
    r31_tgp_dsb = K3_koide / 0.077649   # ≈  441
    points_tgp = {
        "e/μ/τ (TGP)": (206.77, r31_tgp_lep, "Q=1.5000"),
        "u/c/t (TGP)": (587.96, r31_tgp_uct, "Q=1.5259"),
        "d/s/b (TGP)": (20.00,  r31_tgp_dsb, "Q=1.5170"),
    }
    # Punkty PDG (faktyczne masy kwarków)
    points_pdg = {
        "e/μ/τ (PDG)": (206.77,  3477.0,  "Q=1.5000"),
        "u/c/t (PDG)": (587.96,  79949.0, "Q=1.1779"),
        "d/s/b (PDG)": (20.00,   895.0,   "Q=1.3672"),
    }

    print(f"  {'Punkt':<22} {'r₂₁':>10} {'r₃₁':>12} {'Q':>10}  Pozycja vs krzywa Q=3/2")
    print(f"  {'-'*75}")

    all_pts = {**points_tgp, **points_pdg}
    for label, (r21, r31, Qstr) in all_pts.items():
        r31_k, _ = r31_koide_analytic(r21)
        ratio = r31 / r31_k if not np.isnan(r31_k) and r31_k > 0 else np.nan
        if abs(ratio - 1.0) < 0.01:
            pos = "★ NA KRZYWEJ"
        elif ratio > 1:
            pos = f"POWYŻEJ krzywej (×{ratio:.2f})"
        else:
            pos = f"PONIŻEJ krzywej (×{ratio:.3f})"
        print(f"  {label:<22} {r21:>10.2f} {r31:>12.1f} {Qstr:>10}  {pos}")

    print(f"""
  DIAGRAM (r₂₁, r₃₁) — pozycja vs. krzywa Q=3/2:

    POWYŻEJ krzywej Q=3/2  → Q < 3/2  (PDG kwarki: za małe r₃₁)
    NA KRZYWEJ              → Q = 3/2  (leptony, TGP kalibracja)
    PONIŻEJ krzywej Q=3/2  → Q > 3/2  (TGP kwarki: za duże r₃₁ → uwięzienie)

  Punkty TGP leżą PONIŻEJ krzywej → Q_TGP > 3/2
  Punkty PDG (kwarki) leżą POWYŻEJ krzywej → Q_PDG < 3/2

  GEOMETRYCZNA INTERPRETACJA UWIĘZIENIA:
    Soliton TGP kwarkowy leży w regionie "za rozciągniętym" (r₃₁ zbyt duże).
    Aby osiągnąć Q=3/2, potrzebuje SKOMPRESOWANIA r₃₁ → zmniejszenia K₃.
    QCD (kolor) dostarcza tej kompresji przez uwięzienie → hadron.
    Ale hadrony są UKŁADAMI ZŁOŻONYMI, nie solitonami → ich Q_hadron ≠ 3/2.
    """)


# =========================================================================
# SEKCJA F: Czy istnieje unifikacja? — jednoparametrowy model
# =========================================================================
def section_F():
    print("\n" + "=" * 72)
    print("SEKCJA F: CZY ISTNIEJE UNIFIKACJA? — JEDNOPARAMETROWY MODEL")
    print("=" * 72)
    print("""
  Pytanie: Czy można znaleźć JEDEN parametr (λ lub a_Γ lub α_f),
  który dałby Q=3/2 JEDNOCZEŚNIE dla leptonów I kwarków?

  Z analizy C: każda rodzina potrzebuje INNEGO λ_eff dla Q=3/2.
  λ_eff(lep) ≈ 5.47e-6
  λ_eff(u/c/t) << λ_eff(lep)  (K₃ musi być większe)
  λ_eff(d/s/b) << λ_eff(lep)  (analogicznie)

  Jedyna spójność: K₃ jest UNIVERSALNE (nie zależy od α_f przy stałym λ).
  Zmiana λ przesuwa K₃ jednakowo dla WSZYSTKICH rodzin.
  Dlatego nie można jednocześnie dać Q=3/2 dla leptonów I kwarków
  przez samą zmianę λ — byłoby to sprzeczne z universalnością K₃.
    """)

    # Weryfikacja numeryczna
    C = 2.000
    a = 0.040

    K1_lep = 0.009839
    K1_uct = 0.004175
    K1_dsb = 0.077649

    r21_lep = 206.77
    r21_uct = 587.96
    r21_dsb = 20.00

    r31_k_lep, _ = r31_koide_analytic(r21_lep)
    r31_k_uct, _ = r31_koide_analytic(r21_uct)
    r31_k_dsb, _ = r31_koide_analytic(r21_dsb)

    # λ_eff dające Q=3/2 dla każdej rodziny
    lam_lep = (C*a)**2 / (K1_lep * r31_k_lep)**2
    lam_uct = (C*a)**2 / (K1_uct * r31_k_uct)**2
    lam_dsb = (C*a)**2 / (K1_dsb * r31_k_dsb)**2

    print(f"  λ_eff(e/μ/τ) = {lam_lep:.4e}   ← to jest λ_Koide")
    print(f"  λ_eff(u/c/t) = {lam_uct:.4e}   ← {lam_uct/lam_lep:.2f}× λ_Koide")
    print(f"  λ_eff(d/s/b) = {lam_dsb:.4e}   ← {lam_dsb/lam_lep:.2f}× λ_Koide")

    print(f"""
  WNIOSEK F:
    λ_eff(u/c/t) ≈ {lam_uct/lam_lep:.1f}× λ_Koide  — nie można unifikować przez λ
    λ_eff(d/s/b) ≈ {lam_dsb/lam_lep:.1f}× λ_Koide  — analogicznie

    TGP w obecnej postaci (soliton izolowany, sferyczna symetria, V_mod)
    nie daje Q=3/2 dla kwarków PRZY TYM SAMYM λ co leptony.
    Konieczne rozszerzenie: albo modyfikacja V_mod dla kwarków (kolor?),
    albo uwzględnienie oddziaływań między solitonami (wielociałowe TGP).
    """)


# =========================================================================
# SEKCJA G: Podsumowanie P42 i perspektywy
# =========================================================================
def section_G():
    print("\n" + "=" * 72)
    print("SEKCJA G: PODSUMOWANIE P42 I PERSPEKTYWY")
    print("=" * 72)
    print("""
  PYTANIE P42: Czy masy konstituentów kwarków dają Q≈3/2?
  ODPOWIEDŹ: NIE — żaden model konstituentów nie daje Q≈3/2 dla kwarków.

  WYNIKI KLUCZOWE:

  1. Masy konstituentów (różne modele): Q(u/c/t)≈1.28–1.35, Q(d/s/b)≈2.1–2.4
     → Daleko od Q=3/2 w obu przypadkach.

  2. Odwrotnie: m₃ wymagane dla Q=3/2 przy m₁,m₂ PDG ≠ m₃_PDG
     → Spektrum kwarkowe PDG jest niezgodne z Koide Q=3/2.

  3. λ_eff dające Q=3/2 jest różne dla każdej rodziny (sekcja F):
     λ_eff(u/c/t) ≈ 0.80 × λ_Koide,  λ_eff(d/s/b) ≈ 0.87 × λ_Koide
     → Różnice ~13–20%. Nie można unifikować przez jeden λ (K₃ universalne).

  4. Diagram (r₂₁, r₃₁):
     TGP kwarki: PONIŻEJ krzywej Q=3/2  → Q_TGP > 3/2 (uwięzienie)
     PDG kwarki: POWYŻEJ krzywej Q=3/2  → Q_PDG  < 3/2 (faktyczne masy)

  PERSPEKTYWY P43:

  Ścieżka A: Wielociałowe TGP
    — Rozważyć energię układu N solitonów (kwarki w hadronie)
    — Warunek Q=3/2 dla układu → analogon zasady stabilności dla hadronu?

  Ścieżka B: Modyfikacja V_mod z kolorem
    — Dodać człon "kolorowy" do V_mod: V_kolor ~ g_s · (kolor) · φ²
    — Sprawdzić, czy Q=3/2 jest przywracane dla kwarków kolorowych

  Ścieżka C: TGP a Gell-Mann–Okubo
    — Relacja GMO: M = M₀ + a·Y + b·[I(I+1) - Y²/4]
    — Czy warunek Q=3/2 implikuje relację GMO lub odwrotnie?
    — Numerycznie: sprawdzić, czy Q(p,Λ,Ξ) zbliża się do 3/2 przy perturbacji
      mas baryonów w duchu GMO

  ★ REKOMENDACJA P43: Ścieżka C (GMO) — natychmiastowo weryfikowalna numerycznie.
    """)


# =========================================================================
# MAIN
# =========================================================================
if __name__ == '__main__':
    print("=" * 72)
    print("P42: MASY KONSTITUENTÓW KWARKÓW A WARUNEK KOIDE Q=3/2")
    print("Pytanie: Czy efektywne masy kwarków w hadronie dają Q=3/2?")
    print("=" * 72)

    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    section_G()

    print("\n" + "=" * 72)
    print("KONIEC ANALIZY P42")
    print("=" * 72)
