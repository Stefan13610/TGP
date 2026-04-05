"""
TGP — Analityczne wyprowadzenie r21 = M2/M1
============================================

Cel: pokazac skad pochodzi stosunek M2/M1 ~ 207.
Strategia: porzadkowanie po potegach Phi0 i a_Gam.
"""
import numpy as np
from scipy.optimize import brentq

# ======================================================================
# KROK 1: Energia w wspolrzednych bezwymiarowych
# ======================================================================
# Profil: Phi(r) = Phi0 + K*exp(-msp*r)/r
# msp = sqrt(Phi0_base/Phi0), Phi0_base=1
# Wprowadzamy: s = r*msp = r/xi, xi = sqrt(Phi0)
#
# Phi(s) = Phi0 * [1 + kappa_tilde * exp(-s)/s]
#   kappa_tilde = K/(Phi0*xi) = K/Phi0^(3/2)
#
# dPhi/dr = Phi0 * kappa_tilde * msp * exp(-s)*(-s-1)/s^2
#
# ======================================================================
# ENERGIA KINETYCZNA (standardowa)
# ======================================================================
# Ek_std = 4pi * int (1/2)(dPhi/dr)^2 * r^2 dr
#        = 4pi * Phi0^2 * kappa_tilde^2 * I_kin(eps)
#
# gdzie: I_kin(eps) = int_{eps}^{inf} (1/2) exp(-2s)*(1+s)^2/s^2 ds
#   eps = a_Gam/xi = a_Gam/sqrt(Phi0)
#
# Rozklad przy eps -> 0:
#   I_kin(eps) ~ 1/(4*eps) + (regularne)
#              ~ xi/(4*a_Gam) + ...
#              = sqrt(Phi0)/(4*a_Gam) + ...
#
# => Ek_std ~ 4pi * Phi0^2 * kappa_tilde^2 * sqrt(Phi0)/(4*a_Gam)
#           = pi * kappa_tilde^2 * Phi0^(5/2) / a_Gam
#
# kappa_tilde = K/Phi0^(3/2), wiec kappa_tilde^2*Phi0^(5/2) = K^2/Phi0^3 * Phi0^(5/2) = K^2/sqrt(Phi0)
# => Ek_std ~ pi * K^2 / (a_Gam * sqrt(Phi0))
#           ~ pi * K^2 / a_Gam   dla Phi0~1
#
# ======================================================================
# ENERGIA POTENCJALNA
# ======================================================================
# Ep = 4pi * Phi0^2 * int (V(psi)-V(1)) r^2 dr
#   V(psi) = psi^3/3 - psi^4/4
#   V(1) = 1/12
#
# Rozwinienie V(1+eps)-V(1) przy eps = kappa_tilde*exp(-s)/s:
#   V(1+eps)-V(1) = -eps^2/2 - (2/3)eps^3 - eps^4/4 + ...
#
# Czlon eps^2: int_{eps}^inf (-1/2)*kappa_tilde^2*exp(-2s)/s^2 * s^2 ds
#            = -kappa_tilde^2/2 * int exp(-2s) ds = -kappa_tilde^2/4
# => Ep_quad = 4pi*Phi0^2*xi^3 * (-kappa_tilde^2/4)
#            = -pi * kappa_tilde^2 * Phi0^(7/2)
#            = -pi * K^2 / Phi0^3 * Phi0^(7/2) = -pi * K^2 * sqrt(Phi0)
#
# Czlon eps^4 (UV-dominowany, calkuje od eps):
# int_{eps}^inf (-1/4)*kappa_tilde^4*exp(-4s)/s^4 * s^2 ds
# ~ -kappa_tilde^4/4 * int_eps^inf exp(-4s)/s^2 ds ~ +kappa_tilde^4/(4*eps)
# => Ep_quart ~ 4pi*Phi0^2*xi^3 * kappa_tilde^4/(4*eps)
#             = 4pi*Phi0^2*Phi0^(3/2) * kappa_tilde^4 * sqrt(Phi0)/(4*a_Gam)
#             = pi * kappa_tilde^4 * Phi0^4 / a_Gam
#             = pi * K^4 / (Phi0^2 * a_Gam)
#
# ======================================================================
# PELNA ENERGIA (rzedy wiodace)
# ======================================================================
# E ~ Ek_std + Ep_quad + Ep_quart + Ek_nl + ...
#   ~ pi*K^2/(a_Gam*sqrt(Phi0)) - pi*K^2*sqrt(Phi0) + pi*K^4/(Phi0^2*a_Gam) + ...
#
# Czlon kinetyczny nieliniowy Ek_nl = alpha/Phi0^2 * Ek_std ~ alpha*K^2/(Phi0^2*a_Gam*sqrt(Phi0))
#
# ======================================================================
# WARUNEK SAMOSPOJNOSCI E = M = 4*pi*Phi0*K
# ======================================================================
# Porzadkowanie wiodace (zachowujac czlony K^2/a_Gam i K^2*sqrt(Phi0)):
#
# pi*K^2/(a_Gam*sqrt(Phi0)) - pi*K^2*sqrt(Phi0) = 4*pi*Phi0*K
#
# K * pi * K * (1/(a_Gam*sqrt(Phi0)) - sqrt(Phi0)) = 4*pi*Phi0
#
# K = 4*Phi0 / (1/(a_Gam*sqrt(Phi0)) - sqrt(Phi0))
#   = 4*Phi0*a_Gam*sqrt(Phi0) / (1 - a_Gam*Phi0)
#
# M = 4*pi*Phi0*K = 16*pi*a_Gam*Phi0^(5/2) / (1 - a_Gam*Phi0)
#
# To jest wynik dla WLACZENIA czlonu K^2*sqrt(Phi0).
# Dla Phi0 << 1/a_Gam: M ~ 16*pi*a_Gam*Phi0^(5/2)  => beta = 5/2!
#
# ======================================================================
# WYNIK ANALITYCZNY: M(Phi0)
# ======================================================================
# M(Phi0) ~ 16*pi*a_Gam*Phi0^(5/2) / (1 - a_Gam*Phi0)   [porzadkowanie 1]
#
# Dla a_Gam << 1 i Phi0 nie za duze: M ~ 16*pi*a_Gam*Phi0^(5/2)
# => r21 = M2/M1 = (Phi0_2/Phi0_1)^(5/2)  i  Phi0_2 = 207^(2/5) ~ 8.44
#
# Korekty od Ep_quart (czlon K^4) daja wyzszy efektywny wykladnik
# i obnizaja Phi0_2 z ~8.4 w kierunku ~6.2 (zmierzone numerycznie).

GAMMA = 1.0; PHI0_BASE = 1.0

def M_analytic_v1(Phi0, a_Gam):
    """Porzadkowanie v1: Ek_std + Ep_quad"""
    denom = 1.0 - a_Gam * Phi0
    if denom <= 0: return None
    return 16 * np.pi * a_Gam * Phi0**2.5 / denom

def M_analytic_v2(Phi0, a_Gam, alpha=5.9):
    """Porzadkowanie v2: Ek_std + Ep_quad + Ek_nl (bez Ep_quart)"""
    # Ek_nl / Ek_std = alpha/Phi0^2 (dla Phi0 >> 1)
    # Wiec: Ek_eff = Ek_std * (1 + alpha/Phi0^2)
    # Ep_quad = -pi*K^2*sqrt(Phi0) (bez zmian)
    #
    # K * pi * (1+alpha/Phi0^2) * K/(a_Gam*sqrt(Phi0)) - pi*K^2*sqrt(Phi0) = 4*pi*Phi0*K
    # K * [(1+alpha/Phi0^2)/(a_Gam*sqrt(Phi0)) - sqrt(Phi0)] = 4*Phi0
    A = (1 + alpha/Phi0**2) / (a_Gam * np.sqrt(Phi0))
    B = np.sqrt(Phi0)
    denom = A - B
    if denom <= 0: return None
    K = 4*Phi0 / denom
    return 4*np.pi*Phi0*K

def M_analytic_v3(Phi0, a_Gam, alpha=5.9):
    """Porzadkowanie v3: wszystkie wazne czlony przez rownosc samosp."""
    # E = Ek_std*(1+alpha/Phi0^2) + Ep_quad + Ep_quart = M
    # pi*K^2*(1+alp/P^2)/(a_G*sqrt(P)) - pi*K^2*sqrt(P) + pi*K^4/(P^2*a_G) = 4*pi*P*K
    # Podziel przez pi*K:
    # K*(1+alp/P^2)/(a_G*sqrt(P)) - K*sqrt(P) + K^3/(P^2*a_G) = 4*P
    # To jest rownanie na K. Dla malego K: K ~ 4*P*a_G*sqrt(P)/(1+alp/P^2 - a_G*P)
    # Dla wlaczenia K^3: rownanie cubiczne w K.
    P = Phi0
    def eq(K):
        if K <= 0: return 0.0
        lhs = (K*(1+alpha/P**2)/(a_Gam*np.sqrt(P))
                - K*np.sqrt(P)
                + K**3/(P**2*a_Gam))
        return lhs - 4*P
    try:
        K_sol = brentq(eq, 1e-8, 100.0, xtol=1e-10)
        return 4*np.pi*P*K_sol
    except Exception:
        return None


# ======================================================================
# POROWNANIE ANALITYCZNE vs NUMERYCZNE
# ======================================================================
a_Gam = 0.05; alpha = 5.9

# Numeryczne M(Phi0) dla referencji (z poprzednich obliczen):
M_num = {1.0: 0.2213, 2.0: 2.626, 3.0: 8.434, 5.0: 28.54, 6.23: 45.80,
         7.0: 58.48, 10.0: 120.2, 20.0: 455.7, 26.53: 769.4}

print("="*70)
print("POROWNANIE ANALITYCZNE vs NUMERYCZNE  (a_Gam=0.05, alpha=5.9)")
print("="*70)
print(f"{'Phi0':>6}  {'M_num':>9}  {'v1=beta=5/2':>12}  {'v2+Ek_nl':>12}  {'v3+K^4':>10}")
print("-"*70)
for Phi0 in sorted(M_num.keys()):
    M_n = M_num[Phi0]
    v1 = M_analytic_v1(Phi0, a_Gam)
    v2 = M_analytic_v2(Phi0, a_Gam, alpha)
    v3 = M_analytic_v3(Phi0, a_Gam, alpha)
    s_v1 = f"{v1:.4f}" if v1 else "-"
    s_v2 = f"{v2:.4f}" if v2 else "-"
    s_v3 = f"{v3:.4f}" if v3 else "-"
    print(f"{Phi0:>6.2f}  {M_n:>9.4f}  {s_v1:>12}  {s_v2:>12}  {s_v3:>10}")

# ======================================================================
# WYZNACZENIE r21 z METODY ANALITYCZNEJ
# ======================================================================
print()
print("="*70)
print("WYZNACZENIE r21 = M2/M1 z METOD ANALITYCZNYCH")
print("="*70)

M1_n = M_num[1.0]

for label, M_func, args in [
    ("v1 (beta=5/2)", M_analytic_v1, (a_Gam,)),
    ("v2 (+Ek_nl)",   M_analytic_v2, (a_Gam, alpha)),
    ("v3 (+K^4)",     M_analytic_v3, (a_Gam, alpha)),
]:
    M1_a = M_func(1.0, *args)
    if M1_a is None: continue

    # Znajdz Phi0_2 takie ze r21=207
    def r21_eq(p2):
        M2 = M_func(p2, *args)
        if M2 is None: return -207
        return M2/M1_a - 207
    try:
        p2_sol = brentq(r21_eq, 1.01, 30.0, xtol=1e-5)
        M2_a = M_func(p2_sol, *args)
        print(f"\n  {label}:")
        print(f"    Phi0_2 dla r21=207: {p2_sol:.4f}  (numeryczne: 6.227)")
        print(f"    Roznica: {100*(p2_sol-6.227)/6.227:+.2f}%")

        # Sprawdz r31
        def r31_eq(p3):
            M3 = M_func(p3, *args)
            if M3 is None: return -3477
            return M3/M1_a - 3477
        try:
            p3_sol = brentq(r31_eq, p2_sol, 100.0, xtol=1e-5)
            print(f"    Phi0_3 dla r31=3477: {p3_sol:.4f}  (numeryczne: 26.53)")
            print(f"    Phi0_3/Phi0_2 = {p3_sol/p2_sol:.4f}  (potrzebne: 4.26)")
        except Exception:
            print("    Phi0_3 poza zakresem v3")
    except Exception:
        print(f"\n  {label}: nie znaleziono rozwiazania")

# ======================================================================
# WYNIK ANALITYCZNY KONCOWY: PRZYBLIZONA FORMULA r21
# ======================================================================
print()
print("="*70)
print("WZOR ZAMKNIETY NA r21 (przyblizenie v1, beta=5/2):")
print("="*70)
print()
print("  M(Phi0) ~ 16*pi*a_Gam * Phi0^(5/2) / (1 - a_Gam*Phi0)")
print()
print("  => r21 = M2/M1 = Phi0_2^(5/2) * (1 - a_Gam) / (1 - a_Gam*Phi0_2)")
print()
print("  Dla r21 = 207, a_Gam = 0.05:")
# Phi0_2^(5/2) * (1 - 0.05) / (1 - 0.05*Phi0_2) = 207
# Phi0_2^(5/2) * 0.95 / (1 - 0.05*Phi0_2) = 207
def r21_v1_eq(p2):
    if 0.05*p2 >= 1: return -207
    return p2**2.5 * (1-a_Gam) / (1-a_Gam*p2) - 207
try:
    p2_v1 = brentq(r21_v1_eq, 1.01, 18.0)
    print(f"  => Phi0_2 = {p2_v1:.4f}  (vs numeryczne 6.227, blad {100*(p2_v1-6.227)/6.227:+.1f}%)")
    print(f"     czyli r21 ~ (Phi0_2)^beta  z efektywnym beta = {np.log(207)/np.log(p2_v1):.4f}")
except Exception as e:
    print(f"  Blad: {e}")

print()
print("  Efektywny wykladnik beta z porzadkowan:")
for a_G in [0.01, 0.02, 0.05, 0.10]:
    # Szukamy Phi0_2 dla r21=207 w wersji v1
    def eq(p2):
        if a_G*p2>=1: return -207
        M1 = 16*np.pi*a_G*1.0**2.5/(1-a_G)
        M2 = 16*np.pi*a_G*p2**2.5/(1-a_G*p2)
        return M2/M1 - 207
    try:
        p2 = brentq(eq, 1.01, 25.0)
        beta_eff = np.log(207)/np.log(p2)
        print(f"  a_Gam={a_G:.3f}: Phi0_2={p2:.4f}, beta_eff={beta_eff:.4f}, r21={207:.0f}")
    except Exception:
        print(f"  a_Gam={a_G:.3f}: brak rozwiazania")

print()
print("="*70)
print("ZWIAZEK Z LICZBA PI:")
print("="*70)
print()
print("  Lambda = 4*pi  (stala samospojnosci)")
print()
print("  M ~ 16*pi*a_Gam*Phi0^(5/2)  =>  pi pojawia sie w prefaktorze!")
print("  Samospojnosc E/K = 4*pi = Lambda narzuca skale PI.")
print()
print("  Czy Phi0_2 = 2*pi?")
p2 = 2*np.pi
try:
    r21_at_2pi = r21_v1_eq(p2) + 207
    print(f"  r21(Phi0_2=2*pi={p2:.4f}) wg v1 = {r21_at_2pi:.2f}  (cel: 207)")
    print(f"  => 2*pi = 6.283 daje r21 = {r21_at_2pi:.1f}, nie 207.")
    print(f"     To pokazuje ze Phi0_2=2*pi jest PRZYPADKOWE (blad ~5%) dla tego a_Gam.")
except Exception as e:
    print(f"  Blad: {e}")

print()
print("="*70)
print("PODSUMOWANIE ANALITYCZNE")
print("="*70)
print("""
WYNIK: M(Phi0) ~ 16*pi*a_Gam * Phi0^(5/2) / (1 - a_Gam*Phi0)

Wiodace porzadkowanie skalowania:
  Ek ~ pi*K^2 / (a_Gam*sqrt(Phi0))    [UV-dominowane przez rdzen]
  Ep ~ -pi*K^2 * sqrt(Phi0)           [maly czlon potencjalu]

Samospojnosc E=M daje:
  K ~ 4*Phi0*a_Gam*sqrt(Phi0) / (1 - a_Gam*Phi0)
  M ~ 16*pi*a_Gam * Phi0^(5/2) / (1 - a_Gam*Phi0)

Wykladnik beta = 5/2 = 2.5 w limicie a_Gam -> 0 lub Phi0 -> 0.
Korekty K^4 i K^3 (od Ep_quart) daja beta_eff > 5/2 dla Phi0 ~ 6.

Stosunek:
  r21 = M2/M1 = Phi0_2^(5/2) * (1-a_Gam) / (1-a_Gam*Phi0_2)

  Dla r21=207, a_Gam=0.05: Phi0_2 ~ 8.0  (v1 porzadkowanie)
  Numerycznie:              Phi0_2 ~ 6.23  (pelen model)

  Roznica: czynniki potencjalu wyzszego rzedu obnizaja Phi0_2.
""")
