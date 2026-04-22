"""
as03_polarizability_probe.py — jeden ostatni test atomowy TGP

Po as01 (H trywialne) i as02 (A_orb SC NIE są atomowe), zostaje ostatnia
testowalna predykcja: **polaryzowalność atomowa α_pol**.

Dlaczego polaryzowalność?
-------------------------
W TGP atom = soliton w substracie. Zewnętrzne pole E perturbuje substrat →
przesuwa centrum solitonu. Indukowany moment dipolowy:
   p = α_pol · E

W standard QM (sum-over-states lub Dalgarno-Lewis):
   α_pol = 2·∑_n |⟨0|z|n⟩|² / (E_n - E_0)

Dla alkaliów α_pol jest BARDZO duże (miękka wiązana valencja):
   H   : 4.50 a.u.
   Li  : 164.1 a.u.
   Na  : 162.7 a.u.
   K   : 290.6 a.u.
   Rb  : 318.8 a.u.
   Cs  : 399.9 a.u.
   Fr  : 317.0 a.u.

Konkretny test TGP:
-------------------
Jeśli TGP A_tail elektronu walencyjnego ~ dominujący 'ogon' w dipole matrix
element, oczekujemy α_pol ∝ (Z_eff/n)^(-?), z jakąś potęgą.

Standardowy skalowanie: α_pol(alkali) ~ n^7 (Drude-like) — empirycznie:
   Li → Na: n 2→3, ratio 164→163 ≈ 1.00 (NIE skaluje!)
   Na → K:  n 3→4, ratio 163→291 ≈ 1.79
   K → Rb:  n 4→5, ratio 291→319 ≈ 1.10
   Rb → Cs: n 5→6, ratio 319→400 ≈ 1.25

Empirycznie: α_pol ~ const dla Li,Na i rośnie wolno potem. Coś o co TGP
powinien się zatrzymać.

Cel: czy (IE/Ry)^(const) daje α_pol(alkali)? Co ten const mówi o TGP?
"""

import math
import sys
import io
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# Stałe
Ry = 13.605693122994

# ---------------------------------------------------------------------------
# Dane: alkali IE₁ i polaryzowalność (a.u. = a_0³)
# ---------------------------------------------------------------------------
# Źródła: NIST (IE), Mitroy/Safronova 2010 (α_pol)
alkalis = [
    # name, Z, n, IE [eV], α_pol [a.u.]
    ("H",  1,  1, 13.59844,   4.5),
    ("Li", 3,  2,  5.39171, 164.1),
    ("Na", 11, 3,  5.13908, 162.7),
    ("K",  19, 4,  4.34066, 290.6),
    ("Rb", 37, 5,  4.17713, 318.8),
    ("Cs", 55, 6,  3.89390, 399.9),
    ("Fr", 87, 7,  4.07272, 317.0),  # relativistic contraction
]

print("=" * 70)
print("  as03 — Atomic polarizability vs TGP soliton picture")
print("=" * 70)

# ---------------------------------------------------------------------------
# Test 1: skalowanie α_pol z IE
# Standard QM: α_pol ~ 1/IE² (z Dalgarno-Lewis nieskończona suma).
# Dla wodorowego 1s: α_pol_H = 9/2·a_0³ = 4.5, IE = 13.6 eV.
# ---------------------------------------------------------------------------
print("\n[1] α_pol vs IE₁ — testowanie prawa skalowania:")
print(f"    {'atom':<4}{'n':>4}{'IE [eV]':>10}{'α_pol':>10}{'α·IE²':>12}{'α·IE³':>12}{'α^(1/4)/IE':>14}")
for name, Z, n, IE, alpha_pol in alkalis:
    p1 = alpha_pol * IE**2
    p2 = alpha_pol * IE**3
    ratio = (alpha_pol)**(0.25) / IE
    print(f"    {name:<4}{n:>4}{IE:>10.4f}{alpha_pol:>10.2f}{p1:>12.2f}{p2:>12.1f}{ratio:>14.5f}")

# Obs: α·IE² prawie stały dla Na,K,Rb,Cs?
print("\n    Produkty: czy jest niezmiennik?")
products = [alpha_pol * IE**2 for _,_,_,IE,alpha_pol in alkalis]
products_arr = np.array(products)
print(f"    α·IE² mean={products_arr.mean():.1f}, std={products_arr.std():.1f}, CV={100*products_arr.std()/products_arr.mean():.1f}%")

products3 = [alpha_pol * IE**3 for _,_,_,IE,alpha_pol in alkalis]
products3_arr = np.array(products3)
print(f"    α·IE³ mean={products3_arr.mean():.1f}, std={products3_arr.std():.1f}, CV={100*products3_arr.std()/products3_arr.mean():.1f}%")

# Tylko alkaliów (bez H)
alk_only = [(Z,n,IE,a) for name,Z,n,IE,a in alkalis if name != "H"]
prod2_alk = np.array([a*IE**2 for Z,n,IE,a in alk_only])
prod3_alk = np.array([a*IE**3 for Z,n,IE,a in alk_only])
prod25_alk = np.array([a*IE**(2.5) for Z,n,IE,a in alk_only])
print(f"\n    Tylko alkali (Li-Fr): α·IE² CV={100*prod2_alk.std()/prod2_alk.mean():.1f}%")
print(f"                          α·IE^2.5 CV={100*prod25_alk.std()/prod25_alk.mean():.1f}%")
print(f"                          α·IE³ CV={100*prod3_alk.std()/prod3_alk.mean():.1f}%")

# ---------------------------------------------------------------------------
# Test 2: fit log(α_pol) vs log(IE) — wyznaczenie wykładnika
# ---------------------------------------------------------------------------
print("\n[2] Fit log(α) = A + β·log(IE) dla alkaliów:")
ln_IE = np.array([math.log(IE/Ry) for _,_,IE,_ in alk_only])
ln_alpha = np.array([math.log(a) for _,_,_,a in alk_only])
# linear regression
beta, A_const = np.polyfit(ln_IE, ln_alpha, 1)
pred = A_const + beta*ln_IE
r = np.corrcoef(ln_IE, ln_alpha)[0,1]
print(f"    β = {beta:.4f}   (wykładnik)")
print(f"    A = {A_const:.4f}   (offset)")
print(f"    r  = {r:.4f}")
print(f"    → α_pol ≈ {math.exp(A_const):.2f}·(IE/Ry)^{beta:.2f}")

# Z TGP: jeśli masa m ∝ A⁴ i polaryzowalność ~ 1/(m·ω²) z ω ~ IE/ℏ:
# α ~ ℏ²/(m · IE²). m/m_e stała (e⁻ jest jedno). Więc α ~ 1/IE² -- naiwnie.
print("\n    Naiwne TGP/QM przewidywanie: α_pol ~ 1/IE² → β = -2")
print(f"    Obserwowany β = {beta:.2f} (różnica od -2: {beta+2:.2f} ≈ {(beta+2)/2*100:.0f}% relatywnie)")

# ---------------------------------------------------------------------------
# Test 3: hipoteza m_core screening — polaryzowalność ∝ a_eff⁴?
# W klasycznym modelu Drude dla alkaliów: α = a_eff³ (objętość "miękkiej orbity")
# gdzie a_eff = a_0·n*²/Z_ion. Z_ion = 1 (alkali).
# ---------------------------------------------------------------------------
print("\n[3] Drude-like: α_pol ≈ (n*)⁶·a_0³ dla alkali z Z_ion=1")
print(f"    {'atom':<4}{'n*':>8}{'(n*)⁶':>12}{'α_pred/α_obs':>18}")
for name, Z, n, IE, alpha_pol in alkalis:
    n_star = math.sqrt(Ry/IE)
    pred = n_star**6
    if name != "H":
        ratio = pred / alpha_pol
        print(f"    {name:<4}{n_star:>8.4f}{pred:>12.2f}{ratio:>18.4f}")
    else:
        print(f"    {name:<4}{n_star:>8.4f}{pred:>12.2f}  (ref: α_H=4.5)")

# ---------------------------------------------------------------------------
# Test 4: co TGP MÓGŁBY dodać do Drude?
# Jeśli ℏ = πχ·A_tail i A_tail ma skalowanie relatywistyczne (Fr),
# oczekujemy RELATYWISTYCZNEJ poprawki dla Cs/Fr.
#
# Cs/Fr paradoks: Fr ma MNIEJSZE α_pol (317) niż Cs (400), mimo że n=7>6.
# To jest znany efekt: Fr 7s contraction (relatywistyka).
# Czy TGP A_tail ma wbudowane to skalowanie?
# ---------------------------------------------------------------------------
print("\n[4] Cs→Fr relativistic contraction:")
print("    Cs (n=6): α = 399.9 a.u., IE = 3.894 eV")
print("    Fr (n=7): α = 317.0 a.u., IE = 4.073 eV  (↑ IE przy ↑ n — relativistic!)")
print()
print("    W TGP: A_tail(Z) ~ A_tail(Z=0)·(1 - β·Z²/c²)^(1/2)")
print("           (używane w ps25 dla f-orbital: η_f = -1.308)")
print()
# Fr → n*² scale: n*(Fr) = sqrt(Ry/4.073) = 1.830
# naive Drude: (1.830)^6 = 37.5. Obs = 317. Ratio ~ 8.5.
# Cs: n*(Cs) = sqrt(13.6/3.894) = 1.869. (1.869)^6 = 42.5. Obs = 400. Ratio ~ 9.4.
# Obie wymagają ~8-10× enhancement z "drowded state" (szeroki rozkład; okres kwantowania)
print("    Drude naive underestimates by ~9× (both Cs, Fr) — rozsądne ze względu na")
print("    uproszczenie n*⁶ vs pełna suma over n'p states.")

# ---------------------------------------------------------------------------
# Test 5: czy istnieje invariant TGP (niezależny od n)?
# ---------------------------------------------------------------------------
print("\n[5] Szukamy TGP invariant:")
print(f"    {'atom':<4}{'(α·IE⁴)^(1/4)':>18}{'(α·IE²)^(1/3)':>18}{'(α·IE)^(1/4)':>18}")
for name, Z, n, IE, alpha_pol in alkalis:
    v1 = (alpha_pol * IE**4)**(0.25)
    v2 = (alpha_pol * IE**2)**(1/3)
    v3 = (alpha_pol * IE)**(0.25)
    print(f"    {name:<4}{v1:>18.4f}{v2:>18.4f}{v3:>18.4f}")

# ---------------------------------------------------------------------------
# WERDYK
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("  AS03 VERDICT")
print("=" * 70)
print(f"""
  Znaleziony wykładnik: α_pol ∝ IE^{beta:.2f}
  Oczekiwane (naiwne):  α_pol ∝ IE^(-2)  (Drude/QM)
  Odchyłka:             {(beta+2)/2*100:.0f}% relatywnie

  Interpretacja:
  • Poprawka od struktury core (penetracja, exchange) powoduje odchyłki
    od prostego 1/IE².
  • TGP nie ma mechanizmu dostarczającego tej poprawki w chwili obecnej.
  • Standard HF/DFT radzi sobie z tym doskonale. TGP nie dodaje wartości.

  OGÓLNY WNIOSEK from as01-as03 (Atomic Shells Program):
  =====================================================

  1. TGP redukuje się do standardowej QM dla atomów (as01: H 1s trywialne).
  2. TGP NIE derywuje A_orb(s/sp/d/f) z atomów (as02: różnica 7× magnituda).
  3. TGP NIE dodaje wartości do opisu polaryzowalności (as03: potrzebuje HF).

  TGP w swojej obecnej formie opisuje:
  ✓ Masy cząstek (leptony dokładne do 0.013%)
  ✓ QM emergence (79/79 testów)
  ✓ Sektor partikularny (Koide dla leptonów, kwarki, neutrina)
  ✓ Galaxy dynamics (jako ν(y) phenomenology)
  ✓ Nadprzewodnictwo ograniczone (calibrated, 60% d-class closure)

  TGP nie opisuje obecnie:
  ✗ Atomy wieloelektronowe (brak aparatu)
  ✗ Chemia (brak aparatu)
  ✗ A_orb jako pierwsze zasady (empiryczne z SC/ρ)
  ✗ Lu superconductivity (4f descriptor brakuje)

  Następny front do uderzenia:
  • TGP + elektromagnetyzm: czy substrat produkuje EM? Jeśli tak, atom mógłby być
    naturalnym testem. Jeśli nie — TGP jest z definicji poza zakresem chemii.
  • Cohesive energy metali (nie atomów): może tutaj A_orb-jak-fazowy fermion-sea
    ma sens. Ale to cofa się do SC/ρ — krąg zamknięty.
  • Relatywistyczny ansatz dla 4f/5f (kontynuacja P7.5a η_f): lokalnie może zamknąć
    Lu. Ale to nie "chemia", to inny brak w SC.

  OSTATECZNA KONSTATACJA:
  Lit nie daje się TGP łatwiej niż cięższe elementy — BO TGP NIE MA APARATU
  atomowego. Sprawa nie jest w "masie pierwiastka", lecz w fundamentalnym
  zakresie teorii. TGP jest teorią GRAWITACJI + MAS CZĄSTEK, nie teorią atomów.
""")
