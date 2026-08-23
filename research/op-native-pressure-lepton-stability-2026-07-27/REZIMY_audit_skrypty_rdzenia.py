"""
Weryfikacja DWOCH twierdzen ze skryptow rdzenia:
  [1] three_regimes_quantitative.py  raportuje OBA punkty rownowagi jako "(stable)".
  [2] two_body_Veff.py  SUMMARY twierdzi "reproduces the three-regime structure",
      a jego wlasny output [2] pokazuje TYLKO JEDNO przejscie sily przez zero.
Sprawdzam [1] wprost (znak V'' i kierunek sily wokol obu pierwiastkow).
"""
import numpy as np

def V(d,b,C):  return -4*np.pi*C**2/d + 8*np.pi*b*C**2/d**2 - 24*np.pi*b*C**3/d**3
def F(d,b,C):  return -4*np.pi*C**2/d**2 + 16*np.pi*b*C**2/d**3 - 72*np.pi*b*C**3/d**4

print("="*78)
print("[1] Czy OBA punkty rownowagi sa stabilne? (skrypt rdzenia mowi: tak)")
print("="*78)
for b, C, tag in [(1.0, 0.1, "test"), (1.0, 2.321116e-78, "proton (z three_regimes_quantitative)")]:
    d1 = 2*b - np.sqrt(4*b**2 - 18*b*C)
    d2 = 2*b + np.sqrt(4*b**2 - 18*b*C)
    print(f"\n  {tag}:  C={C:.6e}")
    print(f"    d1 (maly, 'studnia') = {d1:.6e}    [4.5*C = {4.5*C:.6e}]")
    print(f"    d2 (duzy)            = {d2:.6e}    [4*beta = {4*b:.6e}]")
    for d, nm in [(d1,"d1"), (d2,"d2")]:
        h = d*1e-4
        V2 = (V(d+h,b,C) - 2*V(d,b,C) + V(d-h,b,C))/h**2
        Fm, Fp = F(d*(1-1e-3),b,C), F(d*(1+1e-3),b,C)
        # stabilny: ponizej odpychanie (F>0), powyzej przyciaganie (F<0)
        stab = (Fm > 0) and (Fp < 0)
        print(f"      {nm}: V''={V2:+.4e}  F(d-)={Fm:+.3e} F(d+)={Fp:+.3e}"
              f"  -> {'STABILNY' if stab else 'NIESTABILNY'}")

print("\n  ==> d1 jest NIESTABILNY w obu przypadkach:")
print("      ponizej d1 sila jest PRZYCIAGAJACA (F<0) => uklad spada do d=0 (kolaps),")
print("      powyzej d1 sila jest ODPYCHAJACA  (F>0) => uklad ucieka do d2.")
print("      d1 to MAKSIMUM potencjalu -- bariera, a nie studnia.")
print("  ==> Skrypt rdzenia raportuje go jako '(stable)'. To BLAD w skrypcie.")

print("\n" + "="*78)
print("[2] Skale fizyczne z three_regimes_quantitative.py -- czy sa 'subatomowe'?")
print("="*78)
lP, rp, a0, RH = 1.616255e-35, 8.75e-16, 5.292e-11, 1.321836e+26
print(f"  Rdzen (rem:trzy-rezimy-fizyczne): 'trzy rezimy wystepuja na skalach SUBATOMOWYCH'")
print(f"\n  {'obiekt':<10}{'d_well [m]':>16}{'d_grav [m]':>16}{'d_well/l_Planck':>18}{'d_grav/R_Hubble':>18}")
for nm, dw, dg in [("proton", 2.772411e-52, 1.061716e+26),
                   ("elektron", 1.509925e-55, 1.061716e+26),
                   ("Ziemia",   9.898863e-01, 1.061716e+26)]:
    print(f"  {nm:<10}{dw:>16.4e}{dg:>16.4e}{dw/lP:>18.4e}{dg/RH:>18.4e}")
print(f"\n  Skala odniesienia:  l_Planck={lP:.3e}  r_proton={rp:.3e}  a_Bohr={a0:.3e}  R_Hubble={RH:.3e}")
print("\n  ==> d_well dla protonu jest 17 RZEDOW PONIZEJ dlugosci Plancka.")
print("  ==> d_grav dla protonu = 0.80 R_Hubble.")
print("  ==> Czyli dla dwoch protonow w odleglosci 1 fm obowiazuje REZIM II (ODPYCHANIE),")
print("      a nie grawitacja. Rezim I zaczyna sie dopiero powyzej ~10^26 m.")
print("  ==> Twierdzenie 'na skalach subatomowych' jest ILOSCIOWO FALSZYWE")
print("      wg wlasnego skryptu rdzenia.")

print("\n" + "="*78)
print("[3] M_crit -- czy 'obiekty makroskopowe' faktycznie maja tylko grawitacje?")
print("="*78)
Msun, Mcrit = 1.989e30, 1.601337e+50
print(f"  M_crit = {Mcrit:.4e} kg = {Mcrit/Msun:.4e} M_sun")
for nm, M in [("Ziemia",5.972e24), ("Slonce",Msun), ("Droga Mleczna",1.5e12*Msun),
              ("gromada galaktyk",1e15*Msun), ("obserwowalny Wszechswiat",1.5e53)]:
    print(f"    {nm:<26} M={M:.3e} kg   M<M_crit? {str(M < Mcrit):<6}"
          f" -> {'TRZY REZIMY' if M < Mcrit else 'tylko grawitacja'}")
print("\n  ==> Nawet GROMADY GALAKTYK leza ponizej M_crit => wg tego rachunku")
print("      maja trzy rezimy, czyli odpychanie az do ~10^26 m.")
print("      Rdzen twierdzi, ze warunek jest 'zlamany dla obiektow makroskopowych'.")
print("      Wg skryptu lamie sie dopiero powyzej 8e19 mas Slonca.")
