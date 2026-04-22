"""
ex162_phi0_from_substrate.py
==============================
G5: Wyprowadzenie Phi_0 z parametrow substratu.

STRUKTURA PROBLEMU:
  Phi_0 ≈ 25 (z obserwacji: rho_DE = Phi_0 * H_0^2 / (96*pi*G_0))
  Pytanie: skad 25?

SCIEZKA WYPROWADZENIA:
  1. Substrat = siatka Z_2 (Ising 3D) z (J, z=6)
  2. T_c = 4.5115 * J/k_B (znane dokladnie)
  3. Wilson-Fisher fixed point: r* = -2.25, u* = 3.92 (znane)
  4. Ginzburg-Landau: F[phi] = integral [K(nabla phi)^2 + u2*phi^2 + u4*phi^4 + u6*phi^6]
  5. Identyfikacja: Phi = phi^2 / phi_ref^2 * Phi_0, psi = Phi/Phi_0
  6. V(psi) = beta/3 * psi^3 - gamma/4 * psi^4, beta = gamma (warunek prozni)

KLUCZOWA RELACJA:
  phi_0^2 = |r|/(2u) w fazie zlamanej (MF)
  Phi_0 = phi_0^2 / phi_ref^2 * Phi_0 ... (cirkularne?)

  Nie! Phi_0 jest NORMALIZACJA: Phi_0 = <Phi>_eq = lim_{T->0} <phi^2> / phi_ref^2 * Phi_0

  Prawidlowo: Phi_0 to SKONCZONY parametr wyznaczany przez:
    - Stale substratu (J, z, a_sub)
    - Procedura coarse-graining (l_cg)
    - Punkt staly WF (r*, u*)

PODEJSCIE NUMERYCZNE:
  Kluczowa estymacja z sek10: Phi_0 = r_0/lambda ~ 25
  gdzie r_0 i lambda sa parametrami GL po renormalizacji.

  Z WF: r_0 ~ |r*| * (T_c - T_sub)/T_c = |r*| (bo T_sub << T_c)
  lambda ~ u* / (4*pi^2 * nu^3) ... (zalezy od normalizacji)

  Probujemy kilka sciezek wyprowadzenia.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

# ===== Stale Isinga 3D =====
Z_COORD = 6              # siatka sc
T_C_ZJ = 4.5115          # T_c / J (Ising 3D, simple cubic)
NU = 0.6300              # wykladnik korelacji
ETA = 0.0363             # wykladnik anomalny
BETA_C = 0.3265          # wykladnik parametru porzadku
GAMMA_C = 1.2372         # podatnosc
DELTA_C = 4.789          # rownanie stanu

# Wilson-Fisher fixed point (w normalnych konwencjach)
# V = r*phi^2/2 + u*phi^4/24 (konwencja symetryczna)
R_STAR = -1.40   # (w konwencji field-theoretic: r* ~ -O(1))
U_STAR = 3.92    # (w konwencji field-theoretic, 4-epsilon)

# Alternatywna normalizacja (Wetterich FRG):
# lambda_2 = r_k (mass), lambda_4 = u_k (quartic)
# Na WF: u*/r*^2 = universal ratio
# Known: u*/(4*pi)^2 ~ 0.025 (w 3D)

print("=" * 72)
print("ex162: Phi_0 z parametrow substratu (G5)")
print("=" * 72)

# ===== Sciezka 1: Mean-field z poprawkami WF =====
print("\n--- Sciezka 1: Mean-field + Wilson-Fisher ---")
print("  W fazie zlamanej (T < T_c):")
print("  <phi^2> = phi_0^2 = |r|/(2u) = (T_c-T)/(T_c) * |r*|/(2*u*) * Lambda_UV^2")
print("  Phi_0 = <phi^2> / phi_ref^2 * Phi_0 ... (tautologia w tej formie)")
print()
print("  Prawidlowa interpretacja:")
print("  Phi_0 to STOSUNEK phi_0^2 / phi_ref^2, gdzie:")
print("    phi_0 = wartosc oczekiwana pola substratu w fazie porzadkowej")
print("    phi_ref = skala referencyjna z coarse-graining")

# phi_0^2 / phi_ref^2 w MF:
# phi_0 = sqrt(|r*|/(2*u*)) * Lambda_UV  (w jednostkach UV cutoff)
# phi_ref^2 ~ 1/(z*a^2)  (z gradientu)
# Phi_0 ~ z * a^2 * |r*|/(2*u*) * Lambda_UV^2

# W 3D: Lambda_UV ~ 1/a (inverse lattice spacing)
# Phi_0 ~ z * |r*|/(2*u*)

phi0_MF = Z_COORD * abs(R_STAR) / (2 * U_STAR)
print(f"\n  Phi_0(MF, sciezka 1) = z * |r*| / (2*u*)")
print(f"    = {Z_COORD} * {abs(R_STAR)} / (2 * {U_STAR})")
print(f"    = {phi0_MF:.4f}")
print(f"  (target: ~25)")

# ===== Sciezka 2: Relacja z wykladnikami krytycznymi =====
print("\n--- Sciezka 2: Wykladniki krytyczne ---")
# W poblizu T_c: <phi^2> ~ (T_c - T)^{2*beta_c}
# Phi ~ phi^2 ~ (T_c - T)^{2*beta_c}
# Phi_0 = lim_{T->0} <phi^2>/phi_ref^2

# Stosunek:
# Phi_0 / Phi(T_c-delta) = (T_c)^{2*beta_c} / delta^{2*beta_c}
# To daje: Phi_0 ~ (T_c/delta_T)^{2*beta_c} * prefactor

# Dla T_sub = 0: delta_T = T_c
# Phi_0 ~ 1^{2*beta_c} * prefactor = prefactor

# Prefactor z amplitudy: B w m ~ B * t^beta, t = (T_c-T)/T_c
# Dla Isinga 3D: B ≈ 1.656 (z MC, w konwencji phi w {-1,+1})
B_ISING = 1.656
print(f"  B (amplituda Isinga 3D) = {B_ISING}")
print(f"  <phi> = B * t^beta_c, t = (T_c-T)/T_c")
print(f"  Phi = <phi^2> ~ B^2 * t^{2*BETA_C}")
print(f"  Phi_0 = B^2 * 1^{2*BETA_C} * (normalizacja)")

phi0_B = B_ISING**2
print(f"  B^2 = {phi0_B:.4f}")
print(f"  (target: ~25 — potrzebna normalizacja!)")

# ===== Sciezka 3: Stosunek uniwersalny =====
print("\n--- Sciezka 3: Stosunki uniwersalne ---")
# W Isingu 3D istnieja uniwersalne BEZWYMIAROWE stosunki:
# R_chi = C+/C- (amplitudy podatnosci) = 4.77
# R_xi = f+/f- (amplitudy dlugosci korelacji) = 1.963
# Rc = alpha * C+ * (f+/a_0)^3 / B^2 ≈ ... (dlugosc korelacji / amplituda)

# Moze Phi_0 jest uniwersalnym stosunkiem?
# Phi_0 = f(R_chi, R_xi, ...) ?

R_CHI = 4.77    # C+/C-
R_XI = 1.963    # f+/f-
A_0 = 0.110     # amplituda ciepla wlasciwego (Ising 3D)

print(f"  Stosunki uniwersalne Isinga 3D:")
print(f"    R_chi = C+/C- = {R_CHI}")
print(f"    R_xi = f+/f- = {R_XI}")
print(f"    A_0 = {A_0}")
print(f"    B (amplituda m) = {B_ISING}")
print(f"    beta_c = {BETA_C}")

# Kombinacja: Phi_0 ~ z * R_chi * R_xi^2 / (4*pi)
phi0_univ = Z_COORD * R_CHI * R_XI**2 / (4*np.pi)
print(f"\n  Phi_0 ~ z * R_chi * R_xi^2 / (4*pi)")
print(f"    = {Z_COORD} * {R_CHI} * {R_XI}^2 / (4*pi)")
print(f"    = {phi0_univ:.4f}")

# Inne kombinacje
combos = {
    "z * B^2 * R_xi^2": Z_COORD * B_ISING**2 * R_XI**2,
    "z * B^2": Z_COORD * B_ISING**2,
    "B^2 * R_chi": B_ISING**2 * R_CHI,
    "z^2 * B^2 / (4*pi)": Z_COORD**2 * B_ISING**2 / (4*np.pi),
    "R_chi * R_xi^3": R_CHI * R_XI**3,
    "z * R_chi": Z_COORD * R_CHI,
    "z * B^2 * R_chi / (4*pi)": Z_COORD * B_ISING**2 * R_CHI / (4*np.pi),
    "(z*B*R_xi)^2 / (8*pi)": (Z_COORD * B_ISING * R_XI)**2 / (8*np.pi),
    "4*pi * B^2": 4*np.pi * B_ISING**2,
    "z * B * R_xi * delta_c / pi": Z_COORD * B_ISING * R_XI * DELTA_C / np.pi,
    "T_c/J * z * B^2 / (2*pi)": T_C_ZJ * Z_COORD * B_ISING**2 / (2*np.pi),
    "T_c/J * B^2 * R_chi / (2*pi)": T_C_ZJ * B_ISING**2 * R_CHI / (2*np.pi),
}

print(f"\n  Systematyczny test kombinacji:")
print(f"  {'Kombinacja':>45s}  {'Wartosc':>10s}  {'delta%':>8s}")
print("  " + "-" * 70)
for label, val in sorted(combos.items(), key=lambda x: abs(x[1]-25)):
    delta = abs(val - 25) / 25 * 100
    mark = " <<<" if delta < 15 else ""
    print(f"  {label:>45s}  {val:10.4f}  {delta:8.2f}%{mark}")

# ===== Sciezka 4: Hipoteza a_Gamma * Phi_0 = 1 =====
print(f"\n--- Sciezka 4: a_Gamma * Phi_0 = 1 ---")
print("  Jesli a_Gamma = a_sub = stala sieciowa substratu,")
print("  to Phi_0 = 1/a_sub (w jednostkach l_P)")
print()
print("  Z kosmologii TGP:")
print("  Omega_Lambda = 1/(36*a_Gamma)")
print("  a_Gamma = 1/Phi_0")
print("  => Omega_Lambda = Phi_0/36")

for Phi0 in [23.84, 24.66, 25.0, 26.0]:
    Omega_L = Phi0 / 36.0
    Omega_m = 1 - Omega_L
    print(f"  Phi_0 = {Phi0:6.2f}: Omega_L = {Omega_L:.4f}, Omega_m = {Omega_m:.4f}")

print(f"\n  PDG: Omega_L = 0.685, Omega_m = 0.315")
print(f"  => Phi_0 = 36 * Omega_L = 36 * 0.685 = {36*0.685:.2f}")
print(f"  DESI DR2: Omega_L = 0.6973 => Phi_0 = {36*0.6973:.2f}")

# ===== Sciezka 5: Renormalizacja grupy (ERG) =====
print(f"\n--- Sciezka 5: Wynik ERG (Wetterich) ---")
print("  W ERG, potencjal na WF fixed point:")
print("  V_k(phi) = lambda_2 * phi^2 + lambda_4 * phi^4/4! + ...")
print("  Na fixed point: lambda_2* ~ -0.186 (w 3D LPA)")
print("  lambda_4* ~ 9.20 (w 3D LPA)")
print("  Minimum: phi_min^2 = -6*lambda_2/lambda_4")

lambda2_star = -0.186   # LPA 3D
lambda4_star = 9.20     # LPA 3D
phi_min2 = -6 * lambda2_star / lambda4_star
print(f"  phi_min^2 = -6 * ({lambda2_star}) / {lambda4_star} = {phi_min2:.6f}")
print(f"  phi_min = {np.sqrt(phi_min2):.6f}")

# W TGP: Phi_0 = phi_min^2 * (Normalizacja z coarse-graining)
# Normalizacja: N ~ (Lambda_UV / k_IR)^{2-eta} ~ (1/a_sub * L_hub)^{2-eta}
# W logarytmach: ln(N) ~ (2-eta) * ln(L_hub/a_sub)
# L_hub/a_sub ~ 10^{60} (Hubble radius / Planck length)

# Ale Phi_0 musi byc SKONCZONY (nie zalezy od L_hub/a_sub)!
# Phi_0 jest wyrazony przez LOKALNE parametry substratu.

print(f"\n  WNIOSEK: phi_min^2 z ERG = {phi_min2:.4f}")
print(f"  To jest bezwymiarowe. Phi_0 ~ phi_min^2 * (czynnik geometryczny)")
print(f"  Potrzebny czynnik ~ {25/phi_min2:.1f} aby uzyskac Phi_0 ~ 25")

# ===== Sciezka 6: Podejscie czysto wymiarowe =====
print(f"\n--- Sciezka 6: Analiza wymiarowa ---")
print("  Phi = [pole substratu]^2 / [skala ref]^2 * Phi_0")
print("  [Phi] = bezwymiarowe w TGP (Phi_0 to skala)")
print()
print("  Z dzialania TGP:")
print("  S = integral Phi_0^4/kappa * psi^4 * [K(psi)(nabla psi)^2 - V(psi)] d^4x")
print("  Kappa = 3/(4*Phi_0) => S ~ Phi_0^5 * ...")
print()
print("  Gravitational coupling: 8*pi*G_N = kappa/Phi_0^2 = 3/(4*Phi_0^3)")
print("  G_N ~ 1/Phi_0^3 w jednostkach TGP")
print("  => Phi_0 ~ (M_P/m_sub)^{2/3}")
print()

# W jednostkach Plancka: G_N = l_P^2/m_P = 1 (naturalnie)
# 3/(4*Phi_0^3) = 8*pi => Phi_0 = (3/(32*pi))^{1/3}
phi0_grav = (3/(32*np.pi))**(1.0/3.0)
print(f"  Z 8*pi*G = 3/(4*Phi_0^3):")
print(f"  Phi_0 = (3/(32*pi))^{1/3} = {phi0_grav:.4f}")
print(f"  ALE: to zaklada G = 1 w jednostkach l_P, co daje Phi_0 ~ 0.2")
print(f"  Brak: Phi_0 ~ 25 wymaga G << 1 (czyli Phi_0 >> 1)")

# Interpretacja: Phi_0 = (M_P^2/(m_sub^2 * N_sites))^{1/3}
# gdzie N_sites to liczba wezlow substratu w objetosci Hubble'a
# Ale to jest cirkularne (N_sites zalezy od kosmologii)

# ===== Wnioski =====
print(f"\n{'='*72}")
print("WNIOSKI ex162")
print(f"{'='*72}")
print(f"""
  PROBLEM G5 (Phi_0 z pierwszych zasad) jest FUNDAMENTALNIE TRUDNY.

  Stan wiadomosci:
  1. Phi_0 ~ 25 z obserwacji (rho_Lambda) [NUM] - dobrze znane
  2. a_Gamma * Phi_0 = 1 [HIP] - wspierana DESI (1.03 sigma)
  3. Phi_0 = 36 * Omega_Lambda [ALGEBRAICZNE] - wynika z TGP
  4. Zadna kombinacja uniwersalnych stosunkow Isinga 3D nie daje 25

  GLOWNE PRZESZKODY:
  - Phi_0 koduje STOSUNEK skal: M_P / m_sub (masa Plancka / masa substratu)
  - Ten stosunek NIE jest uniwersalny — zalezy od MIKROFIZYKI substratu
  - W efektywnej teorii pola: Phi_0 jest parametrem Wilson'owskim
    (renormalizowany coupling na skali IR)
  - Obliczenie wymaga KOMPLETNEJ renormalizacji od UV do IR

  MOZLIWE SCIEZKI:
  A. Hipoteza a_Gamma * Phi_0 = 1 => Phi_0 = 1/a_Gamma (1 parametr)
     -> Wystarczy wyznaczyc a_Gamma z dynamiki kosmologicznej
     -> Ale a_Gamma = Omega_m^{{1/(1+q)}} zalezy od Omega_m (cirkularne)

  B. Phi_0 z punktu stalego FRG w 4D (nie 3D Ising)
     -> Punkt staly 4D kwantowej grawitacji? (program asympt. bezpieczenstwa)
     -> BARDZO spekulatywne

  C. Phi_0 z warunku samoorganizacji (SOC)
     -> Substrat samo-tunuje do stanu krytycznego
     -> Phi_0 wynika z atrakcyjnego fixed point
     -> Wymaga modelu dynamicznego substratu

  STATUS: G5 = OTWARTY. Hipoteza a_Gamma*Phi_0 = 1 redukuje problem
  do jednego parametru (Omega_m lub Omega_Lambda), co jest STANDARDOWYM
  problemem kosmologicznym, nie specyficznym dla TGP.
""")
