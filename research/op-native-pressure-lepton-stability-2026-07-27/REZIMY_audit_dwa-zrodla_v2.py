"""
POPRAWKA do REZIMY_audit_dwa-zrodla.py.
Tam uzylem profilu Phi ~ A/(r+r0) BEZ ekranowania. Calka INT Phi1 Phi2 dla Phi~1/r
rozbiega sie w IR (integrand ~1/r^2, miara ~r^2 dr => INT dr ~ L), wiec mierzylem
OBJETOSC PUDLA, nie zaleznosc od d. Stad falszywe "I2 ~ d^-0.045".

UWAGA -- to samo w sobie jest wynikiem: rdzen (lin. 2813-2816) twierdzi
  "calka dywerguje LOGARYTMICZNIE, ale jest regularyzowana przez nieliniowa
   strukture profilu BLISKO ZRODEL"
Oba czlony tego zdania sa falszywe:
  - rozbieznosc jest LINIOWA w L, nie logarytmiczna,
  - i jest w IR (duze r), wiec regularyzacja BLISKO ZRODEL jej nie usuwa.
Sprawdzam oba twierdzenia jawnie, a potem licze z ekranowaniem Yukawy.
"""
import numpy as np

def I2_box(d, L, N=200, r0=1.0, mu=0.0):
    x = np.linspace(-L, L, N); h = x[1]-x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r1 = np.sqrt((X+d/2)**2 + Y**2 + Z**2); r2 = np.sqrt((X-d/2)**2 + Y**2 + Z**2)
    P1 = np.exp(-mu*r1)/(r1+r0); P2 = np.exp(-mu*r2)/(r2+r0)
    return np.sum(P1*P2)*h**3

print("="*78)
print("[1] Jakiego typu jest rozbieznosc INT Phi1 Phi2 ?  (bez ekranowania, mu=0)")
print("="*78)
print(f"{'L':>8}{'I2(d=2)':>16}{'I2/L':>12}{'I2/ln(L)':>12}")
for L in [10, 20, 40, 80]:
    v = I2_box(2.0, L, N=200)
    print(f"{L:>8}{v:>16.4f}{v/L:>12.4f}{v/np.log(L):>12.4f}")
print("\n  I2/L jest w przyblizeniu STALE => rozbieznosc LINIOWA w L (IR), nie logarytmiczna.")
print("  ==> Twierdzenie rdzenia 'dywerguje logarytmicznie' jest FALSZYWE.")
print("  ==> I regularyzacja 'blisko zrodel' nie moze tego naprawic (to rozbieznosc IR).")

print("\n" + "="*78)
print("[2] Z ekranowaniem Yukawy (mu>0) -- teraz calki sa skonczone. Zaleznosc od d?")
print("="*78)
mu, L, N = 0.5, 40.0, 240
def integrals(d):
    x = np.linspace(-L, L, N); h = x[1]-x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r1 = np.sqrt((X+d/2)**2 + Y**2 + Z**2); r2 = np.sqrt((X-d/2)**2 + Y**2 + Z**2)
    P1 = np.exp(-mu*r1)/(r1+1.0); P2 = np.exp(-mu*r2)/(r2+1.0)
    return np.sum(P1*P2)*h**3, np.sum(P1*P2*(P1+P2))*h**3

ds = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0])
print(f"{'d':>8}{'I2':>16}{'I3':>16}{'E_beta~+I2':>14}{'E_gamma~-I3':>14}")
I2s, I3s = [], []
for d in ds:
    a, b = integrals(d); I2s.append(a); I3s.append(b)
    print(f"{d:>8.2f}{a:>16.6f}{b:>16.6f}{'+'+format(a,'.4f'):>14}{'-'+format(b,'.4f'):>14}")
I2s, I3s = np.array(I2s), np.array(I3s)

print("\n  Wykladniki (fit na d>=2, gdzie asymptotyka rdzenia ma obowiazywac):")
m = ds >= 2.0
p2 = np.polyfit(np.log(ds[m]), np.log(I2s[m]), 1)[0]
p3 = np.polyfit(np.log(ds[m]), np.log(I3s[m]), 1)[0]
print(f"    I2 ~ d^({p2:+.2f})    rdzen twierdzi E_beta  ~ (1/d)ln(d/r0)  [~ -1]")
print(f"    I3 ~ d^({p3:+.2f})    rdzen twierdzi E_gamma ~ -1/d^3         [~ -3]")
print("\n  (Przy ekranowaniu spadek jest WYKLADNICZY, nie potegowy -- fit potegowy")
print("   jest tu tylko wskaznikiem, ze zadna z form '1/d ln d' ani '1/d^3' nie pasuje.)")

print("\n" + "="*78)
print("[3] Kluczowe: stosunek E_gamma/E_beta -- czy gamma MOZE zdominowac przy malych d?")
print("="*78)
print("  E_gamma/E_beta = -(3 gamma/Phi0^2 * I3) / (2 beta/Phi0 * I2)")
print("                 = -(3/2)(gamma/beta)(1/Phi0) * (I3/I2)")
print(f"\n  Przy beta=gamma (warunek prozniowy):  ratio = -(3/2)(1/Phi0)*(I3/I2)")
print(f"{'d':>8}{'I3/I2':>14}{'|ratio| przy Phi0=25':>24}")
for d, a, b in zip(ds, I2s, I3s):
    print(f"{d:>8.2f}{b/a:>14.6f}{1.5*(b/a)/25.0:>24.6f}")
print("\n  ==> I3/I2 MALEJE gdy d maleje (nie rosnie!). A przy Phi0=25 (prop:Lambda-eff)")
print("      stosunek jest rzedu 1e-3 -- czlon kwartyczny jest o 3 rzedy SLABSZY.")
print("      Rdzen (lin. 2839-2843) twierdzi, ze E_gamma DOMINUJE dla d<r_well.")
print("      W tej rachubie NIE dominuje na zadnym d.")
print("\n  Zastrzezenie: to zalezy od profilu. Przy profilu silnie spietrzonym blisko")
print("  zrodla (Phi >> Phi0) stosunek moze sie odwrocic -- ale wtedy jestesmy")
print("  poza zadeklarowanym zakresem waznosci |dPhi/Phi0| <~ 0.2 (lin. 2660-2662).")
print("  Czyli: albo E_gamma nie dominuje, albo funkcjonal nie obowiazuje. Nie da sie mieć obu.")
