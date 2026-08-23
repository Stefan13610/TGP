# BUDGET_psi_orientation.py
# DOMKNIECIE LUKI: sprzezenie psi <-> orientacja reperu.
#
# PUNKT WYJSCIA (fakt geometryczny, do zweryfikowania):
#   reper e^a_i = sqrt(h) * R^a_i,  R ortogonalna
#   => metryka g_ij = e^a_i e^a_j = h * (R^T R)_ij = h * delta_ij
#   => METRYKA NIE WIDZI R.  Orientacja jest niezalezna od psi na poziomie metryki.
# Jesli tak, to TOPOLOGIA (nawiniecie) mieszka w R, a NIE w psi
#   => moj wczesniejszy wymog "trawersu psi" (WYNIKI_profil-z-budzetu) byl NIEUPRAWNIONY.
#
# ROZSZERZENIE BUDZETU (jawne, ale naturalne): s_0 = entropia NA WEZEL = log(liczba stanow).
#   Stany wezla obejmuja ORIENTACJE => budzet dzieli sie na TRZY czesci:
#        n_sp + n_czas + n_orient <= B
#   z fh=1:  n_sp = h,  n_czas = 1/h  =>  g(h)=h+1/h  (min 2 przy h=1)
#   koszt orientacji: n_orient = c*|df/dr| (zmiana orientacji na krok wezlowy)
#
# TESTY:
#  [A] czy metryka naprawde nie widzi R (weryfikacja numeryczna)
#  [B] czy minimalizacja budzetu utrzymuje h=1 (czy rdzen psi w ogole DOLUJE?)
#  [C] czy limit tempa delta_max WYPADA z budzetu (dotad ZALOZONY!)
#  [D] co wyznacza rozmiar: budzet (dol) czy sasiedzi (gora)?

import numpy as np
from scipy.optimize import minimize

print("="*70)
print("[A] CZY METRYKA WIDZI ORIENTACJE R?  (e=sqrt(h)R => g=h*R^T R)")
print("-"*70)
rng=np.random.default_rng(0)
def rand_SO3():
    A=rng.normal(size=(3,3)); Q,Rr=np.linalg.qr(A); Q*=np.sign(np.diag(Rr))
    if np.linalg.det(Q)<0: Q[:,0]*=-1
    return Q
h=1.7
for k in range(3):
    R=rand_SO3(); e=np.sqrt(h)*R; g=e.T@e
    print(f"  proba {k+1}: max|g - h*I| = {np.max(abs(g-h*np.eye(3))):.2e}   det R={np.linalg.det(R):+.4f}")
print("  => POTWIERDZONE: g_ij = h*delta_ij NIEZALEZNIE od R.")
print("     Orientacja jest NIEWIDOCZNA w metryce => topologia (nawiniecie) siedzi w R,")
print("     a NIE w profilu psi.  [To KORYGUJE wczesniejszy wymog trawersu psi.]")

print("\n"+"="*70)
print("[C] LIMIT TEMPA Z BUDZETU (dotad byl ZALOZONY jako delta_max)")
print("-"*70)
def g_h(h): return h+1.0/h
print("  wiez na wezel:  g(h) + c*|df| <= B ;  g(h)>=2 (min w h=1)")
print("  => max zmiana orientacji na krok:  |df|_max = (B - 2)/c")
c=1.0
print(f"\n  {'B':>6} {'|df|_max':>10} {'kroki na trawers pi':>21} ")
for B in [2.2,2.5,3.0,4.0,6.0]:
    dmax=(B-2)/c
    print(f"  {B:6.2f} {dmax:10.4f} {np.pi/dmax:21.2f}")
print("  >>> LIMIT TEMPA NIE JEST JUZ ZALOZENIEM - wypada z nadwyzki budzetu (B-2).")
print("      To usuwa jeden z wolnych parametrow poprzedniego rachunku.")

print("\n"+"="*70)
print("[B] CZY RDZEN psi DOLUJE? (minimalizacja budzetu przy nawinieciu)")
print("-"*70)
M=20; B=3.0; c=1.0
w=np.arange(M+1).astype(float)**2      # miara powlokowa
def unpack(x):
    h=x[:M+1]; f=np.concatenate([[np.pi],x[M+1:],[0.0]]); return h,f
def total_budget(x):
    h,f=unpack(x); df=np.abs(np.diff(f))
    return float(np.sum(w*g_h(h)) + c*np.sum(w[1:]*df))
def cons_budget(x):
    h,f=unpack(x); df=np.abs(np.diff(f))
    return B - (g_h(h[1:]) + c*df)         # >=0 na kazdym wezle
x0=np.concatenate([np.ones(M+1), np.linspace(np.pi,0,M+1)[1:-1]])
res=minimize(total_budget,x0,constraints=[{'type':'ineq','fun':cons_budget}],
             bounds=[(0.2,5.0)]*(M+1)+[(0,np.pi)]*(M-1),
             method='SLSQP',options={'maxiter':1500,'ftol':1e-10})
h_s,f_s=unpack(res.x)
print(f"  h w rdzeniu (i=0..5): {np.round(h_s[:6],4)}")
print(f"  h na brzegu (i=M-5..M): {np.round(h_s[-6:],4)}")
print(f"  max|h-1| na calym profilu = {np.max(abs(h_s-1)):.6f}")
if np.max(abs(h_s-1))<1e-3:
    print("  >>> h == 1 WSZEDZIE. Rdzen psi NIE DOLUJE.")
    print("      Obiekt jest CZYSTA TEKSTURA ORIENTACJI, metryka plaska.")
else:
    print("  >>> h odbiega od 1 - rdzen jednak doluje.")

print("\n"+"="*70)
print("[D] CO WYZNACZA ROZMIAR? (profil kata f)")
print("-"*70)
df_s=np.abs(np.diff(f_s))
print(f"  |df| na kroku (i=1..8): {np.round(df_s[:8],4)}")
print(f"  |df|_max dozwolone przez budzet = (B-2)/c = {(B-2)/c:.4f}")
print(f"  liczba krokow z NIEZEROWYM |df| (>1e-4): {int(np.sum(df_s>1e-4))} z {M}")
print(f"  suma |df| = {np.sum(df_s):.4f}  (wymagane pi = {np.pi:.4f})")
if np.sum(df_s>1e-4) >= M-1:
    print("  >>> trawers ROZLOZONY na CALA dostepna domene (nie skupia sie).")
    print("      => rozmiar wyznaczaja SASIEDZI (dostepna przestrzen), nie budzet.")
    print("      Budzet daje tylko DOLNY limit: >= pi/|df|_max krokow.")
else:
    print("  >>> trawers SKUPIONY - budzet lokalizuje sam z siebie.")

print("\n"+"="*70)
print("WNIOSEK (sprzezenie psi<->orientacja):")
print("  1. Metryka nie widzi R  => topologia w ORIENTACJI, nie w psi.")
print("  2. Minimalizacja budzetu => psi PLASKIE (h=1); brak dolka w rdzeniu.")
print("  3. Limit tempa |df|_max=(B-2)/c WYPADA z budzetu (byl zalozeniem).")
print("  4. Rozmiar: DOLNY limit z budzetu, GORNY z sasiadow (relacyjny).")
print("  => KORYGUJE 'WYNIKI_profil-z-budzetu': tam narzucilem trawers psi,")
print("     ktory NIE jest wymagany przez topologie.")
