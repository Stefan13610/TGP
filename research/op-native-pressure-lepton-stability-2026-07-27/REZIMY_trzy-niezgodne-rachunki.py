"""
two_source_potential.py raportuje 13/13 PASS i test T1a:
   "Trzy rezimy BEZ E_gamma (E_lin+E_beta): 2 zmiany znaku"
Czyli dostaje trzy rezimy z DWOCH czlonow, bez czlonu kwartycznego,
ktory wg rdzenia MA TWORZYC studnie.

Sprawdzam:
 [1] skad biora sie trzy rezimy przy profilu Yukawy z samych E_lin + E_beta,
 [2] czy "studnia" w tym rachunku to cokolwiek innego niz osobliwosc -1/d zrodla,
 [3] zestawienie TRZECH niezgodnych rachunkow tej samej wielkosci.
"""
import numpy as np
from scipy.integrate import quad

m = 1.0   # masa Yukawy (jak w skrypcie rdzenia: 1/m = 1)

def I2_yukawa(d):
    """INT d^3x [e^{-m r1}/r1][e^{-m r2}/r2] -- analitycznie znane, licze numerycznie."""
    def integ(r):
        # calkowanie po kacie dla ustalonego r
        def inner(ct):
            r2 = np.sqrt(r*r + d*d - 2*r*d*ct)
            r2 = max(r2, 1e-9)
            return np.exp(-m*r2)/r2
        ct, w = np.polynomial.legendre.leggauss(80)
        s = sum(wi*inner(ci) for ci, wi in zip(ct, w))
        return 2*np.pi * r*r * (np.exp(-m*r)/max(r,1e-9)) * s
    v, _ = quad(integ, 1e-6, 40.0, limit=300)
    return v

print("="*78)
print("[1] Skad trzy rezimy z DWOCH czlonow? (profil Yukawy)")
print("="*78)
print("  E_lin  ~ -K * e^{-m d}/d      (przyciaganie, zanik ~ e^{-m d})")
print("  E_beta ~ +b * I2(d)           (odpychanie,  I2 ~ e^{-2 m d} -- zanika 2x SZYBCIEJ)")
print()
print(f"  {'d':>7}{'E_lin ~ -e^-d/d':>18}{'I2(d)':>14}{'I2/|E_lin|':>14}")
for d in [0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0]:
    el = np.exp(-m*d)/d
    i2 = I2_yukawa(d)
    print(f"  {d:>7.2f}{-el:>18.5f}{i2:>14.5f}{i2/el:>14.5f}")

print("\n  Mechanizm:")
print("   - MALE d : E_lin ~ -1/d  DYWERGUJE, a I2(d) jest SKONCZONE -> przyciaganie")
print("   - SREDNIE d: I2 jeszcze duze, E_lin juz stlumione -> odpychanie")
print("   - DUZE d : I2 ~ e^{-2md} znika szybciej niz E_lin ~ e^{-md} -> przyciaganie")
print("\n  ==> 'Studnia' w tym skrypcie to PO PROSTU osobliwosc -1/d punktowego zrodla.")
print("      Nie ma w niej nic z czlonu kwartycznego ani z confinementu.")
print("      To ten sam -1/d, ktory na duzych d nazywamy grawitacja.")

print("\n" + "="*78)
print("[2] Czy E_gamma jest w tym rachunku perturbacyjne?")
print("="*78)
print("  Skrypt rdzenia, T4b:  V_bez_gamma(d=1) = +0.9057")
print("                        V_z_gamma (d=1) = -16.8694")
print(f"  zmiana: {abs(-16.8694-0.9057)/abs(0.9057):.1f}x wartosci wyjsciowej, ze zmiana znaku")
print("\n  ==> E_gamma nie jest poprawka -- CALKOWICIE DOMINUJE.")
print("      Rozwiniecie perturbacyjne, na ktorym stoi caly rachunek, tam nie zbiega.")
print("      Skrypt zalicza to jako PASS ('E_gamma poglabia studnie').")
print("      Dodatkowo test liczy I3 tylko dla d >= 0.8 -- male d sa POMINIETE.")

print("\n" + "="*78)
print("[3] TRZY niezgodne rachunki TEJ SAMEJ wielkosci")
print("="*78)
rows = [
    ("thm:three-regimes (eq:scales)", "r_well ~ (alpha/beta) r0", "r_rep ~ (beta/gamma)(qM/Phi0)",
     "r_well < r_rep", "zalezy od alpha (NIEOBLICZONE)"),
    ("prop:trzy-rezimy-beta-gamma",   "d_well = 2b + sqrt(4b^2-18bC)", "d_rep = 2b - sqrt(...)",
     "d_well > d_rep  [ODWROTNIE]", "brak alpha, brak r0"),
    ("two_source_potential.py",       "r_well ~ 0.32/m", "r_rep ~ 6.0/m",
     "r_well < r_rep", "trzy rezimy BEZ E_gamma"),
]
print(f"  {'zrodlo':<32}{'r_well':<32}{'r_rep':<32}")
for r in rows:
    print(f"  {r[0]:<32}{r[1]:<32}{r[2]:<32}")
    print(f"  {'':<32}porzadek: {r[3]:<24}{r[4]}")
print("\n  ==> Trzy rachunki tej samej wielkosci daja:")
print("      - dwa RozNE porzadki skal (r_well<r_rep vs d_well>d_rep),")
print("      - dwie rozne zaleznosci od masy,")
print("      - dwa rozne mechanizmy studni (kwartyczny vs osobliwosc -1/d),")
print("      - jeden zalezny od czlonu, ktory nigdy nie zostal obliczony.")
print("      ZADEN nie jest oznaczony w rdzeniu jako obalajacy pozostale.")
