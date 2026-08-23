# AUDIT4_h_decoupling.py  (v2 - czysta parametryzacja, bez nierozniczkowalnego abs)
# ZARZUT 1: czy "h == 1" w BUDGET_psi_orientation [B] jest WYNIKIEM czy KONSTRUKCJA?
#
# Funkcja celu autora:
#     J(h,f) = SUM_i w_i * g(h_i)  +  c * SUM_i w_i * |df_i|      w_i = i^2
#     wiez:    g(h_i) + c*|df_i| <= B
# h wchodzi WYLACZNIE przez g(h)=h+1/h (min w h=1) i przez wiez, ktorego LUZ
# jest MAKSYMALNY dokladnie przy h=1.  => h=1 jest optimum dla KAZDEGO f.

import numpy as np
from scipy.optimize import minimize, linprog
rng = np.random.default_rng(7)

def g_h(h): return h + 1.0/h
M, B, c = 20, 3.0, 1.0
W = np.arange(M+1).astype(float)**2

print("="*72)
print("[A0] BUG STRUKTURALNY: h[0] (WEZEL RDZENIA) NIE WYSTEPUJE W ZADNYM CZLONIE")
print("-"*72)
print(f"  waga powlokowa w_0 = 0^2 = {W[0]:.0f}  =>  w_0*g(h_0) == 0 dla KAZDEGO h_0")
print("  wiez cons_budget liczony tylko dla h[1:]  (patrz linia: g_h(h[1:]) + c*df)")
print("  => h_0 jest ZMIENNA MARTWA. Optymalizator zwraca WARTOSC POCZATKOWA x0[0]=1.")
print("  Autor raportuje 'h w rdzeniu (i=0..5)' - pierwsza z tych liczb to x0, nie wynik.")
# demonstracja: startujemy z h0 = 3.0
def author_solve(h0_init=1.0):
    def unpack(x):
        h=x[:M+1]; f=np.concatenate([[np.pi],x[M+1:],[0.0]]); return h,f
    def tot(x):
        h,f=unpack(x); df=np.abs(np.diff(f))
        return float(np.sum(W*g_h(h)) + c*np.sum(W[1:]*df))
    def cons(x):
        h,f=unpack(x); df=np.abs(np.diff(f))
        return B - (g_h(h[1:]) + c*df)
    x0=np.concatenate([np.ones(M+1), np.linspace(np.pi,0,M+1)[1:-1]]); x0[0]=h0_init
    r=minimize(tot,x0,constraints=[{'type':'ineq','fun':cons}],
               bounds=[(0.2,5.0)]*(M+1)+[(0,np.pi)]*(M-1),
               method='SLSQP',options={'maxiter':3000,'ftol':1e-12})
    return unpack(r.x)
for h0i in [1.0, 3.0, 0.3]:
    h,f = author_solve(h0i)
    print(f"    x0[0]={h0i:4.1f}  ->  h_0 na wyjsciu = {h[0]:.4f}   (h[1:] max|h-1| = {np.max(abs(h[1:]-1)):.2e})")
print("  >>> POTWIERDZONE: h_0 = cokolwiek wpiszesz. 'h==1 w rdzeniu' to artefakt x0.")

print("\n"+"="*72)
print("[A1] SEPAROWALNOSC: czy optimum h[1:] zalezy w OGOLE od f?")
print("-"*72)
bad=0; maxdev=0.0; ntr=0
dmax=B-2.0
for trial in range(300):
    df=np.zeros(M)
    idx=rng.choice(M,size=rng.integers(4,M+1),replace=False)
    df[idx]=rng.random(len(idx))
    if df.sum()==0: continue
    df*=np.pi/df.sum()
    if df.max()>dmax: continue
    ntr+=1
    obj=lambda hh: float(np.sum(W[1:]*g_h(hh)))
    cons=[{'type':'ineq','fun':lambda hh,df=df: B-(g_h(hh)+c*df)}]
    r=minimize(obj,np.full(M,1.5),constraints=cons,bounds=[(0.2,5.0)]*M,
               method='SLSQP',options={'maxiter':800,'ftol':1e-14})
    dev=np.max(np.abs(r.x-1.0)); maxdev=max(maxdev,dev)
    if dev>1e-4: bad+=1
print(f"  dopuszczalnych prob f: {ntr};  z nich h != 1: {bad};  max|h-1| = {maxdev:.2e}")
print("  >>> h==1 dla KAZDEGO dopuszczalnego f. To TOZSAMOSC konstrukcji, nie odkrycie.")
print("  DOWOD: J = A(h) + C(f) separowalne; g(h)>=2 (AM-GM, rownosc h=1);")
print("         wiez B-g(h)-c|df|>=0 najluzniejszy przy h=1. Brak alternatywy.")

print("\n"+"="*72)
print("[B] MINIMALNE SPRZEZENIE h<->orientacja, ktorego autor NIE mogl pominac")
print("-"*72)
print("  Autor sam ustala e = sqrt(h) R  =>  g_ij = h delta_ij.")
print("  Zatem |dR|^2_g = g^{ij} tr(d_iR^T d_jR) = |df|^2 / h.")
print("  Koszt orientacji NIEZMIENNICZY: c*|df|/sqrt(h) (L1) lub c*|df|^2/h (L2).")
print("  Autor uzyl c*|df| - to koszt liczony w WEZLACH, nie w geometrii, ktora sam zdefiniowal.")

def solve_coupled(coupling, p=1, B=3.0, c=1.0, M=20):
    """zmienne: h_1..h_M (h_0 martwe - pomijam), d_1..d_M >=0, sum d = pi"""
    w=np.arange(1,M+1).astype(float)**2
    def split(x): return x[:M], x[M:]
    def ocost(h,d):
        if coupling=="none":   return c*d**p
        if coupling=="proper": return c*d**p/h**(p/2.0)
    def tot(x):
        h,d=split(x); return float(np.sum(w*g_h(h)) + np.sum(w*ocost(h,d)))
    def cbud(x):
        h,d=split(x); return B-(g_h(h)+ocost(h,d))
    def csum(x):
        h,d=split(x); return np.sum(d)-np.pi
    x0=np.concatenate([np.ones(M), np.full(M,np.pi/M)])
    r=minimize(tot,x0,constraints=[{'type':'ineq','fun':cbud},{'type':'eq','fun':csum}],
               bounds=[(0.25,4.0)]*M+[(0.0,np.pi)]*M,method='SLSQP',
               options={'maxiter':8000,'ftol':1e-12})
    h,d=split(r.x); return h,d,r

for p in [1,2]:
    print(f"\n  --- koszt orientacji ~ |df|^{p} ---")
    for cpl in ["none","proper"]:
        h,d,r=solve_coupled(cpl,p=p)
        print(f"    {cpl:7s} ok={r.success}  max|h-1|={np.max(abs(h-1)):.4f}  h[1:5]={np.round(h[:4],4)}")
        print(f"              d[1:5]={np.round(d[:4],4)}  aktywnych={int(np.sum(d>1e-4))}  sum d={d.sum():.4f}")
print("\n  >>> 'proper' (ta sama geometria, tylko konsekwentnie) daje h != 1 w rdzeniu.")
print("      => wniosek 'rdzen psi NIE doluje / czysta tekstura o plaskiej metryce'")
print("         jest ARTEFAKTEM pominietego sprzezenia. Model nie mogl dac nic innego.")

print("\n"+"="*72)
print("[C] CZY LIMIT TEMPA (B-2)/c JEST ZYSKIEM (usunieciem parametru)?")
print("-"*72)
print("  Przed:  delta_max                      -> 1 wolny parametr")
print("  Po:     delta_max = (B-2)/c, B wolne, c wolne (autor: 'Stala c ... wolna')")
for target in [0.01,0.1,1.0,10.0]:
    print(f"    delta_max={target:7.3f}  <- (B,c)=(3.0,{1.0/target:g});  (B-2)/c={(3.0-2)/(1.0/target):.4f}")
print("  >>> 1 parametr -> 2 parametry o tym samym stopniu swobody w ilorazie.")
print("      REPARAMETRYZACJA, nie wyprowadzenie. Zysk = 0.")
