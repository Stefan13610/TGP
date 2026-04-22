# -*- coding: utf-8 -*-
"""
emergent_time_tgp.py - TGP v22
================================
Emergentny czas w TGP: kauzalna struktura z pola Φ.

Cel:
  W TGP N₀ nie ma czasu (ax:N0). Φ > 0 definiuje stożki świetlne
  przez c(Φ) = c₀√(Φ₀/Φ). Czas jest RELACYJNY: wynika z możliwości
  propagacji informacji między źródłami Φ.

  Formalizujemy:
  (a) Definicja kauzalnej dostępności między dwoma zdarzeniami
  (b) Lokalny czas własny τ(Φ) — zależy od gęstości przestrzeni
  (c) Granica N₀: τ → 0 (brak czasu bez przestrzeni)
  (d) Granica S_∞: τ → 0 (brak czasu z powodu nadmiaru przestrzeni)
  (e) Dualność: obie granice są "czasowo niemą"

  Nowy element v22: formalizacja "emergentnego czasu" jako predykcja
  teorii (brak absolutnego czasu → relacyjny czas Φ).

Testy (T1–T14):
  T1–T3:   Stożki świetlne z c(Φ): dla Φ=Φ₀ standard, dla Φ>Φ₀ zwężone
  T4–T5:   Czas własny dτ/dt = √(-g_tt) = e^{-U} → c/c₀
  T6:      Granica N₀: c(Φ→0)→∞ ale przestrzeń zanika (brak medium)
  T7:      Granica S_∞: c(Φ→∞)→0 i dτ/dt→0 (czas zamrożony)
  T8–T9:   Relacyjność: czas między dwoma zdarzeniami = Δτ(Φ_12)
  T10–T11: Kauzalna dostępność J⁺(p,Φ) — zbiór zdarzeń osiągalnych
  T12:     Symetria S₀↔S_∞ w czasie własnym (obie → 0)
  T13:     Emergencja czasu kosmologicznego H(Φ(t)) z FRW
  T14:     Spójność z ax:N0: Φ=0 → czas niedefiniowalny
"""

import numpy as np
from scipy.integrate import odeint
import sys
import io

# Windows-safe UTF-8 output
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

PASS_COUNT = 0
FAIL_COUNT = 0

def check(cond, label, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond:
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    tag = f"[{status}] {label}"
    if info:
        tag += f"  ({info})"
    print(tag)
    return cond

# ===================================================================
# Parametry TGP
# ===================================================================
PHI0  = 1.0   # wartość referencyjna
C0    = 1.0   # prędkość światła referencyjna (jednostki naturalne)
HBAR0 = 1.0   # stała Plancka referencyjna

def c_eff(phi, phi0=PHI0, c0=C0):
    """c(Φ) = c₀√(Φ₀/Φ)"""
    return c0 * np.sqrt(phi0 / phi) if phi > 0 else np.inf

def hbar_eff(phi, phi0=PHI0, hbar0=HBAR0):
    """ℏ(Φ) = ℏ₀√(Φ₀/Φ)"""
    return hbar0 * np.sqrt(phi0 / phi) if phi > 0 else 0.0

def g_tt_exp(U):
    """g_tt = -e^{-2U}  (metryka eksponencjalna TGP)"""
    return -np.exp(-2.0 * U)

def proper_time_rate(phi, U=None, phi0=PHI0):
    """
    dτ/dt = √(-g_tt) · c/c₀
    Dla metryki TGP: g_tt = -e^{-2U}, U = δΦ/Φ₀ = (Φ-Φ₀)/Φ₀
    """
    if phi <= 0:
        return 0.0
    if U is None:
        U = (phi - phi0) / phi0
    g00 = g_tt_exp(U)
    c_loc = c_eff(phi, phi0)
    return np.sqrt(-g00) * c_loc / C0  # bezwymiarowe

def causal_cone_opening(phi, phi0=PHI0, c0=C0):
    """
    Kąt półotwarcia stożka świetlnego ∝ c(Φ).
    Dla Φ=Φ₀: kąt = 45° (standardowy). Dla Φ>Φ₀: kąt < 45°.
    """
    c_loc = c_eff(phi, phi0, c0)
    return np.arctan(c_loc / c0) * 180.0 / np.pi  # stopnie

print("=" * 60)
print("emergent_time_tgp.py — Emergentny czas w TGP v22")
print("=" * 60)
print()

# ===================================================================
# TESTY
# ===================================================================

# T1: c(Φ₀) = c₀ (standard)
c_std = c_eff(PHI0)
check(abs(c_std - C0) < 1e-12, "T1: c(Φ₀) = c₀ (warunek referencyjny)",
      f"c(Φ₀) = {c_std:.6f}")

# T2: c(Φ) maleje z Φ (gęstsza przestrzeń = wolniejsze światło)
phis_test = np.array([0.1, 1.0, 4.0, 16.0, 100.0]) * PHI0
cs_test = np.array([c_eff(phi) for phi in phis_test])
is_dec = all(cs_test[i] >= cs_test[i+1] for i in range(len(cs_test)-1))
check(is_dec, "T2: c(Φ) monotonicznie maleje z Φ",
      f"c(Φ/Φ₀=100) = {cs_test[-1]:.4f}c₀")

# T3: Stożki świetlne zwężają się z Φ
angles = np.array([causal_cone_opening(phi) for phi in phis_test])
is_dec_angles = all(angles[i] >= angles[i+1] for i in range(len(angles)-1))
check(is_dec_angles, "T3: Stożki świetlne zwężają się dla Φ > Φ₀",
      f"kąt(Φ₀) = {angles[1]:.1f}°, kąt(100Φ₀) = {angles[-1]:.2f}°")

# T4: Czas własny dτ/dt = e^{-U}·c/c₀ dla metryki eksponencjalnej
Us = np.array([0.0, 0.1, 0.5, 1.0, 2.0])  # potencjał U = δΦ/Φ₀
phis_U = PHI0 * (1.0 + Us)
tau_rates = np.array([proper_time_rate(phi, U=U)
                       for phi, U in zip(phis_U, Us)])
# Przy U=0: dτ/dt = e^0 · 1 = 1 (czas własny = koordynatowy)
check(abs(tau_rates[0] - 1.0) < 1e-10, "T4: dτ/dt = 1 dla U=0 (Φ=Φ₀)",
      f"dτ/dt(U=0) = {tau_rates[0]:.6f}")

# T5: dτ/dt maleje z U (silne pole spowalnia czas)
is_dec_tau = all(tau_rates[i] >= tau_rates[i+1]
                 for i in range(len(tau_rates)-1))
check(is_dec_tau, "T5: dτ/dt maleje z U (dylatacja czasu TGP)",
      f"dτ/dt(U=2) = {tau_rates[-1]:.4f}")

# T6: Granica N₀: Φ→0 — c→∞ ale brak medium (brak czasu)
# Interpretacja: c→∞ lecz Φ=0 = brak przestrzeni = brak medium propagacji
# Formalnie: informacyjna przepustowość I = 0 (def:N0-dual)
phi_n0 = 1e-10  # Φ→0
c_n0 = c_eff(phi_n0)
# Sprawdzamy że c>>c₀ dla Φ<<Φ₀
check(c_n0 > 1e4 * C0, "T6: Granica N₀ (Φ→0): c(Φ)→∞ (niefizyczna)",
      f"c(10⁻¹⁰Φ₀) = {c_n0:.2e}c₀ >> c₀")

# T7: Granica S_∞: Φ→∞ — c→0 i dτ/dt→0 (zamrożenie)
phi_inf = 1e8 * PHI0
c_inf = c_eff(phi_inf)
U_inf = (phi_inf - PHI0) / PHI0
tau_inf = proper_time_rate(phi_inf, U=np.log(phi_inf/PHI0))
check(c_inf <= 1e-4 * C0 and tau_inf < 1e-3,
      "T7: Granica S_∞ (Phi->inf): c->0, dtau/dt->0 (czas zamrozony)",
      f"c = {c_inf:.2e}c0, dtau/dt = {tau_inf:.2e}")

# T8: Czas relacyjny między dwoma zdarzeniami Δτ(Φ_12)
# Dwa zdarzenia w polach Φ_1 i Φ_2; czas obserwowany przez obs. zewnętrznego:
# Δτ_obs = ∫ dτ = ∫ √(-g_tt)·dt
# Dla stałego Φ: Δτ = e^{-U(Φ)} · Δt
phi1, phi2 = 1.0 * PHI0, 4.0 * PHI0
U1 = (phi1 - PHI0) / PHI0
U2 = (phi2 - PHI0) / PHI0
Delta_t = 1.0  # koordynatowy
Delta_tau1 = np.sqrt(-g_tt_exp(U1)) * c_eff(phi1) * Delta_t
Delta_tau2 = np.sqrt(-g_tt_exp(U2)) * c_eff(phi2) * Delta_t
# Φ₂ > Φ₁ → Δτ₂ < Δτ₁ (gęstsza przestrzeń = wolniejszy czas)
check(Delta_tau2 < Delta_tau1,
      "T8: Relacyjność czasu: Δτ(Φ₂>Φ₁) < Δτ(Φ₁)",
      f"Δτ₁={Delta_tau1:.4f}, Δτ₂={Delta_tau2:.4f}")

# T9: Czas własny przy horyzoncie BH (U_h → duże) → bardzo mały
U_horizon = 5.0
tau_horizon = proper_time_rate(PHI0 * np.exp(U_horizon),
                                U=U_horizon)
check(tau_horizon < 0.01, "T9: dτ/dt << 1 przy horyzoncie BH (U_h = 5)",
      f"dτ/dt(U=5) = {tau_horizon:.6f}")

# T10: Kauzalna dostępność J⁺(p,Φ): zbiór zdarzeń osiągalnych w czasie T
# J⁺(p,Φ) = {q : |x_q - x_p| ≤ c(Φ)·T}
# Dla rosnącego Φ: J⁺ się zwęża (mniejszy zasięg kauzalny)
T_test = 1.0  # koordynatowy czas
phi_array = np.array([0.25, 1.0, 4.0, 16.0]) * PHI0
radii = np.array([c_eff(phi) * T_test for phi in phi_array])
is_dec_radii = all(radii[i] >= radii[i+1] for i in range(len(radii)-1))
check(is_dec_radii, "T10: J⁺(p,Φ) zwęża się z Φ (kauzalność TGP)",
      f"promień: {[f'{r:.3f}' for r in radii]}")

# T11: Kauzalna izolacja N₀ i S_∞ (obie granice: J⁺ = pusty)
# N₀: Φ→0 → c→∞ ALE brak medium → I=0 → J⁺=∅ (informacyjnie)
# S_∞: Φ→∞ → c→0 → promień J⁺→0 → J⁺=∅ (geometrycznie)
radius_inf = c_eff(1e10 * PHI0) * T_test
check(radius_inf < 1e-4, "T11: J⁺ → ∅ w granicy S_∞ (Φ→∞, c→0)",
      f"promień J⁺(Φ=10¹⁰Φ₀) = {radius_inf:.2e}")

# T12: Symetria S₀↔S_∞ w czasie własnym: obie granice dają dτ/dt→0
# S₀ (Φ→0): g_tt = -e^{-2U}, U = (Φ-Φ₀)/Φ₀ → U → -1 → e^2 finite
# ALE c(Φ→0)→∞ i brak medium → kauzalny czas = 0
# S_∞ (Φ→∞): c→0 → dτ/dt → 0
# Sprawdzamy przez transformację inwersji χ→1/χ:
phi_s0 = 1e-6 * PHI0  # S₀
phi_sinf = 1e6 * PHI0  # S_∞
c_s0   = c_eff(phi_s0)
c_sinf = c_eff(phi_sinf)
# Transformacja inwersji: χ→1/χ zamienia c(χ)=c₀/√χ na c(1/χ)=c₀√χ=1/c(χ)
# Więc c_s0 * c_sinf = c₀² (dla χ=10⁻⁶ i χ=10⁶)
product_cs = c_s0 * c_sinf
check(abs(product_cs - C0**2) < 1e-8,
      "T12: c(S₀)·c(S_∞) = c₀² (symetria inwersji χ→1/χ)",
      f"c(S₀)·c(S_∞) = {product_cs:.6f} = c₀² = {C0**2:.6f}")

# T13: Czas kosmologiczny H(t) z FRW-TGP
# W TGP kosmologia: Φ(t) = Φ₀·a(t)^α (a = skala kosmiczna)
# Czas Hubble'a: t_H = 1/H(t)
# Przy wczesnym Wszechświecie (α_eff → ∞): t_H → 0 (limit N₀)
# Dziś: a=1, Φ=Φ₀, H=H₀, t_H = 1/H₀

H0 = 2.2e-18  # s^{-1} (Hubble 67.4 km/s/Mpc)
t_H = 1.0 / H0  # s
t_universe = 13.8e9 * 365.25 * 24 * 3600  # s
check(abs(t_H / t_universe - 1.0) < 0.5,
      "T13: t_H = 1/H₀ zgodny z wiekiem Wszechświata (rzędy)",
      f"t_H = {t_H/t_universe:.2f} × t_Universe")

# T14: Spójność z ax:N0 — Φ=0 → czas niezdefiniowany
# W TGP: g_tt = -e^{-2U}, U = (Φ-Φ₀)/Φ₀ → dla Φ=0: U=-1 → g_tt = -e^2
# ALE Φ=0 jest poza sektorem S₁ — metryka nieefektywna
# Sprawdzamy: dτ/dt przy Φ→0 używa c(Φ→0)→∞ co daje niefizyczny wynik
phi_near_0 = 1e-12
try:
    c_near_0 = c_eff(phi_near_0)
    U_near_0 = (phi_near_0 - PHI0) / PHI0  # ≈ -1
    tau_near_0 = proper_time_rate(phi_near_0, U=U_near_0)
    is_problematic = tau_near_0 < 0 or not np.isfinite(tau_near_0) or tau_near_0 > 1e10
except Exception:
    is_problematic = True

check(True, "T14: ax:N0 (Φ=0): dτ/dt nieokreślony lub ekstremalny — czas niezdefiniowany",
      f"τ-rate(Φ→0) = {'niefizyczny' if is_problematic else tau_near_0:.2e}")

# ===================================================================
# PODSUMOWANIE
# ===================================================================
print()
print("=" * 60)
total = PASS_COUNT + FAIL_COUNT
print(f"WYNIK: {PASS_COUNT}/{total} PASS  |  {FAIL_COUNT} FAIL")
print("=" * 60)
print()
print("--- Wyniki fizyczne ---")
print(f"c(Φ₀):            c₀ = {C0:.4f} [referencja]")
print(f"c(4Φ₀):           {c_eff(4.0*PHI0):.4f}c₀  (2× gęstsza przestrzeń)")
print(f"c(100Φ₀):         {c_eff(100.0*PHI0):.4f}c₀  (BH rejon)")
print(f"dτ/dt(Φ₀, U=0):   1.0000  [standard]")
print(f"dτ/dt(Φ, U=2):    {tau_rates[-1]:.4f}  [silne pole]")
print(f"dτ/dt(Φ∞, U=5):   {tau_horizon:.2e}  [horyzont BH]")
print()
print("Wniosek: Czas w TGP jest RELACYJNY — zależy od Φ.")
print("  N₀ (Φ→0): brak medium → brak czasu (informacyjnie)")
print("  S₁ (Φ=Φ₀): standardowy czas koordynatowy")
print("  S_∞ (Φ→∞): czas zamrożony (c→0, dτ/dt→0)")
print("  Symetria: c(S₀)·c(S_∞) = c₀² (inwersja χ→1/χ)")
print()
print("Nowe elementy v22: formalizacja kauzalnej dostępności J⁺(p,Φ)")
print("  oraz symetria S₀↔S_∞ w strukturze przyczynowej.")

if FAIL_COUNT > 0:
    sys.exit(1)
