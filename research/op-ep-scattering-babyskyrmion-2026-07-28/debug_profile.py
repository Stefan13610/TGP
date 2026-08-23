#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Find a good soliton size by 1-parameter energy minimization (robust, no collapse)."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.optimize import minimize_scalar
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1
Rmax, M = 24.0, 1200
r = np.linspace(0, Rmax, M + 1); dr = r[1] - r[0]

def energy_of(f):
    f = f.copy(); f[0] = np.pi; f[-1] = 0.0
    E = bs.profile_energy_1d(f, dr, KAPPA2, MU2)
    fp = np.gradient(f, dr); s = np.sin(f); rr = r.copy(); rr[0] = 1e-9
    Es = 2*np.pi*np.sum(rr[1:-1]*dr*0.5*KAPPA2*fp[1:-1]**2*s[1:-1]**2/rr[1:-1]**2)
    Ep = 2*np.pi*np.sum(rr[1:-1]*dr*MU2*(1-np.cos(f[1:-1])))
    return E, Es, Ep

def fam_exp(w):   return np.pi*np.exp(-r/w)
def fam_atan(w):  return 4*np.arctan(np.exp(-r/w))

for name, fam in (("exp", fam_exp), ("atan", fam_atan)):
    res = minimize_scalar(lambda w: energy_of(fam(w))[0], bounds=(0.5, 8.0), method='bounded')
    w = res.x; E, Es, Ep = energy_of(fam(w))
    print(f"{name:5s}: w*={w:.3f}  E={E:.4f}  Eskyrme={Es:.4f}  Epot={Ep:.4f}  virial={Es/Ep:.3f}")

# now full L-BFGS but STARTING from the best atan ansatz, WITH a size that already satisfies virial
from scipy.optimize import minimize
res = minimize_scalar(lambda w: energy_of(fam_atan(w))[0], bounds=(0.5, 8.0), method='bounded')
f0 = fam_atan(res.x); f0[0] = np.pi; f0[-1] = 0.0
def build(x):
    f = np.empty(M+1); f[0]=np.pi; f[-1]=0.0; f[1:M]=x; return f
r2 = minimize(lambda x: bs.profile_energy_1d(build(x), dr, KAPPA2, MU2),
              f0[1:M], jac=lambda x: bs.profile_grad_1d(build(x), dr, KAPPA2, MU2)[1:M],
              method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-12, gtol=1e-9))
E, Es, Ep = energy_of(build(r2.x))
print(f"full : E={E:.4f}  Eskyrme={Es:.4f}  Epot={Ep:.4f}  virial={Es/Ep:.3f}  success={r2.success}")
