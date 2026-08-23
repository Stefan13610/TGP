"""
AUDYT ssec:dwa-zrodla (sek08_formalizm, lin. 2578-2916).
Nie ufam zacytowanym asymptotykom -- licze calki nakladania WPROST.

eq:Ebeta   E_beta  = (2 beta/Phi0)   * INT Phi1 Phi2            > 0   (odpychanie)
eq:Egamma  E_gamma = -(3 gamma/Phi0^2) * INT Phi1 Phi2 (Phi1+Phi2) < 0 (studnia)

Rdzen twierdzi:
  E_beta  ~ (1/d) ln(d/r0)     [z regularyzacja "przez nieliniowa strukture profilu"]
  E_gamma ~ -1/d^3             [BEZ zadnej wzmianki o regularyzacji]
Sprawdzam oba na regularyzowanym profilu Phi = A/(r + r0).
"""
import numpy as np

A, r0, L, N = 1.0, 1.0, 60.0, 260

def integrals(d):
    """INT Phi1 Phi2  oraz  INT Phi1 Phi2 (Phi1+Phi2)  na siatce 3D."""
    x = np.linspace(-L, L, N); h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r1 = np.sqrt((X + d/2)**2 + Y**2 + Z**2)
    r2 = np.sqrt((X - d/2)**2 + Y**2 + Z**2)
    P1 = A/(r1 + r0); P2 = A/(r2 + r0)
    I2 = np.sum(P1*P2) * h**3
    I3 = np.sum(P1*P2*(P1+P2)) * h**3
    return I2, I3

ds = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
print("="*78)
print("Calki nakladania (profil regularyzowany Phi = A/(r+r0), r0=1)")
print("="*78)
print(f"{'d':>8}{'I2=INT P1P2':>16}{'I3=INT P1P2(P1+P2)':>22}{'I3/I2':>12}")
I2s, I3s = [], []
for d in ds:
    I2, I3 = integrals(d)
    I2s.append(I2); I3s.append(I3)
    print(f"{d:>8.2f}{I2:>16.5f}{I3:>22.5f}{I3/I2:>12.5f}")
I2s, I3s = np.array(I2s), np.array(I3s)

print("\n" + "="*78)
print("Wykladniki potegowe (fit na duzych d, gdzie asymptotyka ma obowiazywac)")
print("="*78)
m = ds >= 4.0
p2 = np.polyfit(np.log(ds[m]), np.log(I2s[m]), 1)[0]
p3 = np.polyfit(np.log(ds[m]), np.log(I3s[m]), 1)[0]
print(f"  I2 ~ d^({p2:+.3f})   rdzen: ~ (1/d)*ln(d/r0)  => wykladnik ~ -1 (z log)")
print(f"  I3 ~ d^({p3:+.3f})   rdzen: E_gamma ~ -1/d^3   => oczekiwane -3")
print(f"\n  ==> I3 skaluje sie jak d^({p3:+.2f}), a NIE jak d^-3.")
print("      Zacytowane '-1/d^3' NIE wynika z tej calki.")

print("\n" + "="*78)
print("Zachowanie przy d -> 0: czy calki dywerguja?")
print("="*78)
print(f"{'d':>8}{'I2':>16}{'I3':>22}")
for d in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.0]:
    I2, I3 = integrals(d)
    print(f"{d:>8.4f}{I2:>16.5f}{I3:>22.5f}")
print("\n  ==> Z regularyzacja r0 obie calki sa SKONCZONE przy d=0.")
print("      Czyli E_gamma NIE leci do -inf, o ile profil jest regularyzowany.")
print("      Ale wtedy rowniez nie ma '-1/d^3' z prop:trzy-rezimy-beta-gamma.")

print("\n" + "="*78)
print("SPRZECZNOSC SKAL: thm:three-regimes vs prop:trzy-rezimy-beta-gamma")
print("="*78)
print("  thm:three-regimes  eq:scales :  r_rep ~ (beta/gamma)*(qM/Phi0)   r_well ~ (alpha/beta)*r0")
print("  prop:trzy-rezimy             :  d_rep = 2b - sqrt(4b^2-18bC)     d_well = 2b + sqrt(4b^2-18bC)")
print()
b = 1.0
print(f"  {'C':>8}{'d_rep':>12}{'d_well':>12}{'d_rep<d_well?':>16}{'r_well~(a/b)r0':>18}")
for C in [0.02, 0.05, 0.10, 0.20]:
    dr = 2*b - np.sqrt(4*b**2 - 18*b*C)
    dw = 2*b + np.sqrt(4*b**2 - 18*b*C)
    print(f"  {C:>8.3f}{dr:>12.5f}{dw:>12.5f}{str(dr<dw):>16}{'~r0 (stale)':>18}")
print("\n  Sprzecznosci:")
print("   1. prop: d_well > d_rep.  Ale tabela ssec:diagram i thm:three-regimes")
print("      wymagaja r_well < r_rep.  => porzadek skal ODWROTNY.")
print("   2. eq:scales: r_well ~ (alpha/beta)*r0 -- NIEZALEZNE od masy M.")
print("      prop:      d_well = 2b + sqrt(4b^2-18bC) -- MALEJE z masa (przez C).")
print("      => rozna zaleznosc od masy dla tej samej wielkosci.")
print("   3. eq:scales: r_well zalezy od ALPHA. Ale w eq:Eint-decomp czlon E_alpha(d)")
print("      NIGDY NIE JEST OBLICZONY (brak paragrafu (iv)), a prop:trzy-rezimy")
print("      uzywa tylko 3 czlonow bez alpha.  => skala zalezy od nieobliczonego czlonu.")

print("\n" + "="*78)
print("ZAKRES WAZNOSCI zadeklarowany przez sam rdzen (lin. 2660-2662)")
print("="*78)
print("  'Funkcjonal E_int z f(Phi) jest wiarygodny jedynie w rezimie Phi >~ 0.8 Phi0")
print("   (tj. |dPhi/Phi0| <~ 0.2).'")
print()
print("  A rezim III (lin. 2839): 'Dla malych d Phi_1,2 sa DUZE w regionie nakladania'")
print("  => rezim III z definicji wymaga |dPhi/Phi0| ~ 1, czyli 5x poza zakresem.")
print()
print("  Dodatkowo f(Phi)=1+4ln(Phi/Phi0) ZERUJE SIE przy Phi/Phi0 = e^(-1/4) = 0.779")
print(f"    e^(-1/4) = {np.exp(-0.25):.5f}  -- dokladnie na granicy deklarowanego zakresu 0.8.")
print("  Ponizej: f<0 => czlon kinetyczny zmienia znak (duch) w tym funkcjonale.")
print("\n  ==> Studnia (rezim III) jest wyprowadzona z funkcjonalu W OBSZARZE,")
print("      w ktorym rdzen SAM deklaruje go za niewiarygodny, i w ktorym jego")
print("      czlon kinetyczny jest ujemny.")
