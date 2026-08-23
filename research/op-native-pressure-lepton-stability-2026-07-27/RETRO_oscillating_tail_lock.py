"""
RETROSPEKCJA 2026-08-16: czy oscylujacy ogon ODE (status_map O-L5, v47b:
'linearyzacja wokol prozni daje r. Bessela sferycznego z omega=1,
oscylacje, nie zanik wykladniczy'; uniwersalne, niezalezne od alpha)
implikuje BLOKOWANIE PAR na dyskretnych odleglosciach?

Mechanizm (teoria liniowa, wazna w ogonie |dPhi| << Phi0 -- czyli DOKLADNIE
tam, gdzie funkcjonal wg audytu ssec:dwa-zrodla JEST wiarygodny):
  E_int(d) = - q^2 * Re G(d),  G(d) ~ e^{-kappa d} cos(d + phi) / (4 pi d)
Jesli ogon oscyluje -> E_int oscyluje -> nieskonczenie wiele MINIMOW
(stabilnych rownowag) w odstepie ~pi (w jednostkach rdzenia r0!).

TEST KTORY MOZE NIE PRZEJSC:
  [T1] ogon oscylacyjny: E_int(d) ma >=3 minima z V''>0 w d in (2, 30)
  [T2] kontrola negatywna, Yukawa (bez cos): ZERO minimow w (2, 30)
       -- jesli Yukawa tez da minima, mechanizm jest artefaktem procedury
  [T3] odstep kolejnych minimow == pi z dokladnoscia <2% (dla kappa<=0.3)
  [T4] lancuch 3 zrodel (suma parowa, teoria liniowa): istnieje stabilna
       konfiguracja rownoodlegla (Hessian 2x2 dodatnio okreslony)
"""
import numpy as np

def E_pair(d, kappa, phi=0.0, osc=True):
    if osc:
        return -np.exp(-kappa*d)*np.cos(d+phi)/d
    return -np.exp(-kappa*d)/d

def find_minima(kappa, osc=True, lo=2.0, hi=30.0, n=200001):
    d = np.linspace(lo, hi, n)
    E = E_pair(d, kappa, osc=osc)
    mins = []
    for i in range(1, n-1):
        if E[i] < E[i-1] and E[i] < E[i+1]:
            # V'' numerycznie
            h = d[1]-d[0]
            Vpp = (E[i+1]-2*E[i]+E[i-1])/h**2
            if Vpp > 0:
                mins.append((d[i], E[i], Vpp))
    return mins

print("="*70)
print("[T1] ogon oscylacyjny e^{-kd} cos(d)/d")
ok_T1 = True
for kappa in (0.0, 0.1, 0.3, 0.5):
    m = find_minima(kappa, osc=True)
    print(f"  kappa={kappa}: {len(m)} minimow; pierwsze 4 d* = "
          f"{[round(x[0],3) for x in m[:4]]}")
    if len(m) < 3:
        ok_T1 = False
print("  T1:", "PASS" if ok_T1 else "FAIL")

print("="*70)
print("[T2] KONTROLA NEGATYWNA: Yukawa e^{-kd}/d (musi byc 0 minimow)")
ok_T2 = True
for kappa in (0.0, 0.1, 0.3, 0.5):
    m = find_minima(kappa, osc=False)
    print(f"  kappa={kappa}: {len(m)} minimow")
    if len(m) != 0:
        ok_T2 = False
print("  T2:", "PASS (Yukawa nie blokuje -> mechanizm = oscylacja)"
      if ok_T2 else "FAIL (artefakt procedury!)")

print("="*70)
print("[T3] odstep minimow vs 2*pi")
# KOREKTA (2026-08-16): pierwotna specyfikacja mowila 'odstep pi' -- BLAD
# autora testu, nie mechanizmu. Minima -cos(d)/d leza tam, gdzie cos(d)~1,
# czyli co 2*pi; pomiedzy nimi sa MAKSIMA (cos=-1). Pierwszy przebieg
# FAIL byl artefaktem zlej predykcji; oczekiwanie poprawione na 2*pi
# PRZED jakimkolwiek uzyciem wyniku.
ok_T3 = True
for kappa in (0.0, 0.1, 0.3):
    m = find_minima(kappa, osc=True)
    gaps = np.diff([x[0] for x in m[:6]])
    rel = np.abs(gaps - 2*np.pi)/(2*np.pi)
    print(f"  kappa={kappa}: odstepy={np.round(gaps,4)}, "
          f"max odch. od 2pi = {100*rel.max():.2f}%")
    if rel.max() > 0.02:
        ok_T3 = False
print("  T3:", "PASS" if ok_T3 else "FAIL")

print("="*70)
print("[T4] lancuch 3 zrodel na osi: x=(-a, 0, +b); E = suma parowa")
kappa = 0.1
def E3(a, b):
    return (E_pair(a, kappa) + E_pair(b, kappa) + E_pair(a+b, kappa))
# start w poblizu pierwszego minimum pary
m1 = find_minima(kappa)[0][0]
from scipy.optimize import minimize
res = minimize(lambda v: E3(v[0], v[1]), x0=[m1, m1], method='Nelder-Mead',
               options=dict(xatol=1e-10, fatol=1e-14))
a, b = res.x
h = 1e-5
H = np.zeros((2,2))
f0 = E3(a,b)
H[0,0] = (E3(a+h,b)-2*f0+E3(a-h,b))/h**2
H[1,1] = (E3(a,b+h)-2*f0+E3(a,b-h))/h**2
H[0,1] = H[1,0] = (E3(a+h,b+h)-E3(a+h,b-h)-E3(a-h,b+h)+E3(a-h,b-h))/(4*h**2)
ev = np.linalg.eigvalsh(H)
ok_T4 = (ev > 0).all()
print(f"  rownowaga: a={a:.4f}, b={b:.4f} (para: d*={m1:.4f}); "
      f"wart. wlasne Hessianu = {np.round(ev,6)}")
print(f"  symetria a==b: |a-b| = {abs(a-b):.2e}")
print("  T4:", "PASS (stabilny lancuch)" if ok_T4 else "FAIL")

print("="*70)
print("UWAGA o skali: d* w JEDNOSTKACH RDZENIA SOLITONU (r0 ODE), bo omega=1")
print("w tych jednostkach -- NIE w jednostkach beta kalibrowanych kosmologia.")
print("WYNIK:", "4/4 PASS" if (ok_T1 and ok_T2 and ok_T3 and ok_T4)
      else "SA FAILE -- mechanizm do odrzucenia lub poprawy")
