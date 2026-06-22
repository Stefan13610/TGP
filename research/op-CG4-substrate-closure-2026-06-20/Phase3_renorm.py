#!/usr/bin/env python3
"""
Phase 3 — renormalizacja operatora zlozonego O_ab = d_a s d_b s (kanal TT spin-2 = sigma_ab).
Pytanie decydujace (C-D): czy wspolczynnik p^2 w <O^TT O^TT>(p) ma scheme-independent continuum,
czy rozbieznosc UV jest niezerowa (-> C_sigma = wolny parametr UV, nie predykcja).

Akcja: free scalar, propagator 1/(k^2+m^2) (3D euklidesowy).
<O_ab(p) O_cd(-p)>_conn = int_k W_ab W_cd /[(k^2+m^2)((p-k)^2+m^2)],  W_ab=1/2(k_a q_b+k_b q_a), q=p-k.
Projekcja TT (p||z): |W^TT|^2 = P_abcd W_ab W_cd.
"""
import sympy as sp
import numpy as np

# ============================================================
# CZESC A — STRUKTURA ANALITYCZNA (sympy)
# ============================================================
print("="*70); print("A. STRUKTURA ANALITYCZNA (sympy)"); print("="*70)
kx,ky,kz,pz,m,k,th = sp.symbols('kx ky kz pz m k theta', real=True)

# W_ab = 1/2(k_a q_b + k_b q_a), q = p - k = (-kx,-ky,pz-kz)
qx,qy,qz = -kx,-ky,pz-kz
def W(a,b):
    ka=[kx,ky,kz][a]; kb=[kx,ky,kz][b]; qa=[qx,qy,qz][a]; qb=[qx,qy,qz][b]
    return sp.Rational(1,2)*(ka*qb+kb*qa)
# TT (p||z): plaszczyzna xy, symetryczna bezsladowa. |W^TT|^2 = 1/2(Wxx-Wyy)^2 + 2 Wxy^2
Wxx,Wyy,Wxy = W(0,0),W(1,1),W(0,1)
WTT2 = sp.simplify(sp.Rational(1,2)*(Wxx-Wyy)**2 + 2*Wxy**2)
print("|W^TT|^2 =", WTT2, "  (UWAGA: niezalezne od pz!)")
# w sferycznych: kx=k sin th cos ph, ky=k sin th sin ph; |W^TT|^2 usrednione po ph
ph=sp.symbols('phi',real=True)
WTT2_sph = WTT2.subs({kx:k*sp.sin(th)*sp.cos(ph), ky:k*sp.sin(th)*sp.sin(ph)})
WTT2_phavg = sp.simplify(sp.integrate(WTT2_sph,(ph,0,2*sp.pi))/(2*sp.pi))
print("<|W^TT|^2>_phi =", WTT2_phavg, "  = 1/2 k^4 sin^4(theta)")

# p^2-coefficient of 1/((p-k)^2+m^2) at p=0 (p along z):
Dp = (pz-kz)**2 + kx**2+ky**2 + m**2
c_p2 = sp.simplify(sp.diff(Dp**-1, pz, 2).subs(pz,0)/2)
print("\nd^2/dpz^2 [1/D_p]|_0 / 1 (p^2 coeff) =", sp.simplify(c_p2))
print("   = (4 kz^2 - k^2 - m^2)/(k^2+m^2)^3   [k^2=kx^2+ky^2+kz^2]")

# Asymptotyka duzego k wspolczynnika p^2 calki:
#   integrand_p2 ~ (1/2) * <|W^TT|^2>_phi * (p^2coeff) / (k^2+m^2)
#   (pierwszy propagator 1/(k^2+m^2)); duze k:
print("\n--- stopien rozbieznosci wspolczynnika p^2 (C_sigma) ---")
print("[O_ab]=3 (d=3) => [<OO>(p)]=3 => [coeff p^2]=1 => LINIOWA rozbieznosc (~Lambda)")
# wspolczynnik rozbieznosci liniowej = kat. calka asymptotyczna:
mu=sp.symbols('mu',real=True) # mu=cos theta
# duze k: (1/2 k^4 sin^4)*(k^2(4mu^2-1))/k^8 = (1/2) sin^4 (4mu^2-1)/k^2 ; * 1/(k^2)->k^{-2}*k^{-2}? sprawdz:
# integrand = 1/2 * <|W^TT|^2> * (4kz^2-k^2-m^2)/(k^2+m^2)^3 * 1/(k^2+m^2)
#           = 1/2 * (1/2 k^4 sin^4) * (k^2(4mu^2-1)-m^2)/(k^2+m^2)^4
# duze k -> 1/2 * 1/2 k^4 sin^4 * k^2(4mu^2-1)/k^8 = 1/4 sin^4(4mu^2-1)/k^2
ang = sp.integrate( (1-mu**2)**2*(4*mu**2-1), (mu,-1,1))  # sin^4 = (1-mu^2)^2
print("kat. calka  int_{-1}^{1}(1-mu^2)^2(4mu^2-1)dmu =", ang, "= -16/35  (NIEZEROWA!)")
print("=> wspolczynnik rozbieznosci LINIOWEJ jest NIEZEROWY")
print("=> C_sigma z golego babla jest scheme-DEPENDENT (rosnie ~Lambda); brak continuum scheme-indep.")

# ============================================================
# CZESC B — WERYFIKACJA NUMERYCZNA: C_sigma(Lambda) vs Lambda
# ============================================================
print("\n"+"="*70); print("B. WERYFIKACJA NUMERYCZNA  C_sigma(Lambda)"); print("="*70)
mval=1.0
def Csigma_of_Lambda(Lam, Nk=4000, Nth=400):
    kk=np.linspace(1e-4,Lam,Nk); tt=np.linspace(0,np.pi,Nth)
    K,TH=np.meshgrid(kk,tt,indexing='ij')
    sin=np.sin(TH); mu=np.cos(TH)
    WTT2avg=0.5*K**4*sin**4
    p2coeff=(K**2*(4*mu**2-1)-mval**2)/(K**2+mval**2)**4  # = (4kz^2-k^2-m^2)/D0^4 z kz=k mu
    integ=0.5*WTT2avg*p2coeff                    # czynnik 1/2 z F_TT
    # miara: d^3k/(2pi)^3 = k^2 sin th dk dth dphi/(2pi)^3 ; phi->2pi
    meas=K**2*sin*(2*np.pi)/(2*np.pi)**3
    f=integ*meas
    dk=kk[1]-kk[0]; dth=tt[1]-tt[0]
    return np.sum(f)*dk*dth
print(f"{'Lambda':>8} {'C_sigma(Lam)':>16}")
Ls=[5,10,20,40,80,160]; Cs=[]
for L in Ls:
    c=Csigma_of_Lambda(L); Cs.append(c)
    print(f"{L:8.1f} {c:16.5f}")
Cs=np.array(Cs); Ls=np.array(Ls,float)
# fit C = a*Lambda + b
A=np.vstack([Ls,np.ones_like(Ls)]).T
(slope,intercept),*_=np.linalg.lstsq(A,Cs,rcond=None)
print(f"\nfit C_sigma(Lam) = {slope:.6f}*Lambda + {intercept:.5f}")
print(f"przewidziany wsp. liniowy = 1/2 * 1/2 * (-16/35) * 2pi/(2pi)^3 = "
      f"{0.5*0.5*(-16/35)*(2*np.pi)/(2*np.pi)**3:.6f}")
print(f"|slope| {'>> 0 => LINIOWO ROZBIEZNE (scheme-dependent)' if abs(slope)>1e-3 else '~0 => zbiezne'}")

print("\n"+"="*70)
print("WERDYKT Phase 3 (wyliczony):")
print(" C_sigma (wsp. p^2 kanalu TT sigma_ab) ma NIEZEROWA rozbieznosc LINIOWA w Lambda.")
print(" => NIE istnieje scheme-independent continuum golego operatora d_s d_s.")
print(" => C_sigma jest UV-czuly = wolny parametr (wymaga warunku renormalizacyjnego: residuum bieguna),")
print("    NIEustalany przez sam babel. C-D = GAP (nie PARTIAL-survival).")
print(" To TLUMACZY pasmo lattice-MC #31 [0.04,11.1]: wartosc skaluje sie z odcieciem sieci.")
print("="*70)
