"""
AUDYT prop:trzy-rezimy-beta-gamma (sek03_rezimy, eq:Veff-beta-eq-gamma).
Sprawdzam TRZY rzeczy u zrodla, nie ufajac narracji:
  [A] czy algebra warunku beta > 9C/2 i pierwiastkow sie zgadza,
  [B] czy struktura znakow F rzeczywiscie daje przyciaganie-odpychanie-studnia,
  [C] czy "studnia" ma DNO (stabilne minimum), czy jest bezdenna (kolaps).
"""
import numpy as np

# V_eff(d) = -4pi C^2/d + 8pi beta C^2/d^2 - 24pi beta C^3/d^3
def V(d, b, C):  return -4*np.pi*C**2/d + 8*np.pi*b*C**2/d**2 - 24*np.pi*b*C**3/d**3
def F(d, b, C):  return -4*np.pi*C**2/d**2 + 16*np.pi*b*C**2/d**3 - 72*np.pi*b*C**3/d**4

print("="*76)
print("[A] Algebra: F=0  <=>  d^2 - 4*beta*d + 18*beta*C = 0")
print("="*76)
for b, C in [(1.0, 0.10), (1.0, 0.20), (5.0, 0.50), (1.0, 0.2222), (1.0, 0.30)]:
    disc = 16*b**2 - 72*b*C
    cond = b > 4.5*C
    print(f"  beta={b:<5} C={C:<7} beta>9C/2? {str(cond):<6} Delta={disc:+.4f}"
          f"  {'2 pierwiastki' if disc>0 else 'BRAK (1 rezim)'}")
    if disc > 0:
        r1 = 2*b - np.sqrt(4*b**2 - 18*b*C)   # mniejszy
        r2 = 2*b + np.sqrt(4*b**2 - 18*b*C)   # wiekszy
        # weryfikacja: F(r)=0 ?
        print(f"        pierwiastki: d1={r1:.6f} (F={F(r1,b,C):+.2e}), "
              f"d2={r2:.6f} (F={F(r2,b,C):+.2e})")

print("\n" + "="*76)
print("[B] Struktura znakow F  (F<0 = przyciaganie, F>0 = odpychanie)")
print("="*76)
b, C = 1.0, 0.10
r1 = 2*b - np.sqrt(4*b**2 - 18*b*C); r2 = 2*b + np.sqrt(4*b**2 - 18*b*C)
print(f"  beta={b}, C={C}:  d1={r1:.4f}, d2={r2:.4f}")
for d, opis in [(0.05,"d << d1"), (0.2,"d < d1"), (r1,"= d1"), (1.0,"d1 < d < d2"),
                (r2,"= d2"), (10.0,"d > d2"), (100.0,"d >> d2")]:
    f = F(d, b, C)
    print(f"    d={d:<8.4f} {opis:<14} F={f:+.4e}  -> "
          f"{'PRZYCIAGANIE' if f<0 else ('ODPYCHANIE' if f>0 else 'rownowaga')}")

print("\n  ==> struktura: przyciaganie / odpychanie / przyciaganie  (3 rezimy) POTWIERDZONA")
print("  ==> ALE: nazwy w rdzeniu. eq:d-rep = 2b - sqrt(...) = MNIEJSZY pierwiastek")
print("           eq:d-well = 2b + sqrt(...) = WIEKSZY pierwiastek")
print("      a tabela w ssec:diagram mowi r_well < r_rep.  => ETYKIETY ZAMIENIONE.")

print("\n" + "="*76)
print("[C] Czy 'studnia' ma DNO?  Zachowanie V(d) przy d -> 0")
print("="*76)
print(f"  {'d':>10}{'V(d)':>16}{'F(d)':>16}")
for d in [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001]:
    print(f"  {d:>10.4f}{V(d,b,C):>16.4e}{F(d,b,C):>16.4e}")
print("\n  V(d) ~ -24*pi*beta*C^3 / d^3  ->  -INF  gdy d -> 0")
print("  F(d) ~ -72*pi*beta*C^3 / d^4  ->  -INF  (przyciaganie rosnie bez ograniczenia)")
print("\n  ==> STUDNIA JEST BEZDENNA. Brak stabilnego minimum, brak 'plaskiego dna'.")
print("      To nie jest confinement -- to KOLAPS do d=0.")
print("      Rdzen (ssec:studnia) twierdzi: 'asymptotyczna swoboda: wewnatrz studni")
print("      gradienty sa male (plaskie dno)'. W eq:Veff-beta-eq-gamma TEGO NIE MA.")

print("\n" + "="*76)
print("[D] Gdzie sa ekstrema V i jakiego typu?")
print("="*76)
for d, nm in [(r1,"d1 (mniejszy)"), (r2,"d2 (wiekszy)")]:
    h = 1e-6
    V2 = (V(d+h,b,C) - 2*V(d,b,C) + V(d-h,b,C))/h**2
    typ = "MINIMUM (stabilne)" if V2 > 0 else "MAKSIMUM (niestabilne)"
    print(f"  {nm:<16} d={d:.6f}  V={V(d,b,C):+.4e}  V''={V2:+.4e}  -> {typ}")
print("\n  ==> jedyne stabilne minimum lezy przy WIEKSZYM pierwiastku (granica")
print("      odpychanie/grawitacja), a nie w 'studni'. Studnia jest zboczem do -INF.")

print("\n" + "="*76)
print("[E] Zakres waznosci: czy wolno uzywac tego rozwiniecia przy malych d?")
print("="*76)
print("  Szereg -C^2/d + beta*C^2/d^2 - beta*C^3/d^3 to rozwiniecie w potegach 1/d.")
print("  Kolejne czlony rosna wzgledem poprzedniego jak ~ beta/d oraz ~ C/d.")
for d in [10.0, 2.0, 1.0, 0.5, 0.2]:
    t1 = abs(4*np.pi*C**2/d); t2 = abs(8*np.pi*b*C**2/d**2); t3 = abs(24*np.pi*b*C**3/d**3)
    print(f"    d={d:<6.2f} |E_lin|={t1:.3e}  |E_beta|={t2:.3e} ({t2/t1:6.2f}x)  "
          f"|E_gamma|={t3:.3e} ({t3/t2:6.2f}x)")
print("\n  ==> przy d < ~2*beta kolejne czlony PRZEWYZSZAJA poprzednie => szereg")
print("      rozbiezny w rezimie II i III. Rezim III lezy POZA zakresem waznosci")
print("      wzoru, ktorym jest wyprowadzany.")
