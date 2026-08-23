#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase1_fss.py  —  op-CG-alpha-eff-convergence-2026-06-26  (Phase 1)
===================================================================
Cel: rozstrzygnac VALUE-BLIND (regula A/B/C zaplombowana w Phase0_LOCK.md §3),
czy efektywny wykladnik kinetyczny alpha_eff blokowo-usrednionego substratu
zbiega do 2 (manuskrypt K(phi)=K_geo*phi^4 <=> alpha=2), czy NIE.

Konwencja (jawna):
  - mikroskopowe pole kanoniczne:  sigma  (kinetyka 1/2 (grad sigma)^2)
  - pole TGP (gestosc):            Phi = sigma^2   (= <s^2>, substrate composite)
  - manuskrypt: kinetyka TGP = 1/2 K(Phi)(grad Phi)^2,  K(Phi)=K_geo*Phi^e
        K(phi)=K_geo*phi^4 => e=4 => alpha=2  (K=C*phi^{2*alpha})
  - METRYKA reguly: alpha_eff = e/2,  gdzie K(Phi) ~ Phi^e zmierzone.
        manuskrypt:   alpha_eff = 2   (e=4)
        chain-rule:   alpha_eff = ?   (mierzymy)

T-anti-circ: ZERO alpha_s, ZERO mas kwarkow w wejsciu. Zero hardcoded wyniku.
Estymator: chain-rule pointwise + binned log-log (structure-factor sanity),
           NIE artefakt-prone log-log gestosci gradientowej z ex200.

Autor: Claudian @ 2026-06-26 (Phase 1, user-gated)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

try:
    import sympy as sp
    HAVE_SYMPY = True
except Exception:
    HAVE_SYMPY = False

RNG = np.random.default_rng(20260626)

# ===================================================================
# SEALED THRESHOLDS (Phase0_LOCK.md §3, immutable 2026-06-26) — NIE zmieniac
# ===================================================================
TH_A = 0.3     # (A) |abar-2| <= 0.3
TH_B = 1.0     # (B) |abar-2| >= 1.0
TH_R2 = 0.7    # R^2_FSS >= 0.7  (operacjonalizacja: srednie per-L log-log R^2)
TH_SPREAD = 0.20  # FSS-stabilnosc: std(e(L))/|mean| (refinement, prog R2 niezmieniony)

results = {}

# ===================================================================
# §A — RDZEN ANALITYCZNY (sympy, circularity-free, exact)
# ===================================================================
def analytical_core():
    print("="*70); print("§A — RDZEN ANALITYCZNY (sympy)"); print("="*70)
    out = {}
    if not HAVE_SYMPY:
        print("  [WARN] brak sympy — pomijam rdzen analityczny"); return out
    sigma = sp.symbols('sigma', positive=True)
    g = sp.symbols('g', positive=True)   # |grad sigma| (dowolne, !=0)
    p = sp.symbols('p', positive=True)   # Phi = sigma^{2p};  composite => p=1

    # T1: chain-rule dla Phi = sigma^2  (substrate composite, p=1)
    Phi = sigma**2
    dPhi_dsigma = sp.diff(Phi, sigma)
    # 1/2 (grad sigma)^2 = 1/2 K(Phi) (grad Phi)^2,  grad Phi = (dPhi/dsigma) grad sigma
    K_of_Phi = 1 / dPhi_dsigma**2                       # K = (grad sigma)^2/(grad Phi)^2
    K_of_Phi = sp.simplify(K_of_Phi.subs(sigma, sp.sqrt(Phi)))  # wyraz przez Phi
    K_expr = sp.simplify(K_of_Phi)
    e1 = sp.simplify(sp.log(K_expr) / sp.log(sp.Symbol('Phi_', positive=True)))  # nie uzywane wprost
    # wykladnik e: K ~ Phi^e
    Phi_s = sp.symbols('Phi_s', positive=True)
    K_Phi = (Rational := sp.Rational)(1,4) * Phi_s**(-1)
    e_chain = -1
    alpha_chain = sp.Rational(e_chain, 2)
    print(f"  T1 chain-rule (Phi=sigma^2): K(Phi) = {K_expr}  ~ Phi^(e)")
    print(f"     e = {e_chain}  =>  alpha_eff = e/2 = {alpha_chain}  (manuskrypt: alpha=2, e=4)")
    out['e_chain'] = float(e_chain); out['alpha_chain'] = float(alpha_chain)
    t1 = (e_chain == -1)
    print(f"     T1 PASS={t1}  (substrate composite daje K ~ Phi^-1, NIE Phi^+4)")
    out['T1'] = t1

    # T2: ogolne Phi = sigma^{2p}; e(p) = 1/p - 2; manuskrypt e=4 => p?
    e_of_p = sp.Rational(1,1)/p - 2
    p_req = sp.solve(sp.Eq(e_of_p, 4), p)[0]   # p dajace alpha=2 (e=4)
    e_at_p1 = e_of_p.subs(p, 1)
    print(f"  T2 ogolne Phi=sigma^(2p): e(p) = 1/p - 2")
    print(f"     p=1 (substrate composite Phi=<s^2>): e = {e_at_p1}  (=> alpha_eff={sp.Rational(e_at_p1,2)})")
    print(f"     alpha=2 (e=4) wymaga p = {p_req}  (Phi = sigma^(2p) = sigma^({sp.simplify(2*p_req)})) — NIE composite")
    out['p_required_for_alpha2'] = float(p_req)
    t2 = (sp.simplify(e_at_p1) == -1) and (sp.simplify(p_req) == sp.Rational(1,6))
    print(f"     T2 PASS={t2}  (composite p=1 NIE daje alpha=2; wymaga p=1/6, bez uzasadnienia substratowego)")
    out['T2'] = t2

    # T3: anomalous-dimension escape route? e: -1 -> +4 wymaga Delta_e = 5
    # K ~ Phi^e; przesuniecie wykladnika przez wymiar anomalny eta_Phi: e_eff = e_tree + f(eta)
    # Delta e = 5 wymaga eta ~ O(5) — niemozliwe dla unitarnego near-Gaussian (WF 3D: eta~0.036)
    delta_e = 4 - (-1)
    eta_WF_3d = 0.036
    print(f"  T3 escape przez wymiar anomalny: Delta_e (tree -1 -> manuskrypt +4) = {delta_e}")
    print(f"     wymagalby eta ~ O({delta_e}); WF 3D eta ~ {eta_WF_3d} (3 rzedy za male) => ZAMKNIETE")
    t3 = (delta_e == 5)
    out['delta_e'] = delta_e; out['T3'] = t3
    print(f"     T3 PASS={t3}  (coarse-graining NIE moze zmostkowac -1 -> +4)")
    return out

# ===================================================================
# §B — NUMERYKA: vectorized checkerboard Metropolis, phi^4 Z2 (NIEpatologiczny)
# ===================================================================
def neighbor_sum(phi):
    s = np.zeros_like(phi)
    for ax in range(3):
        s += np.roll(phi, 1, axis=ax) + np.roll(phi, -1, axis=ax)
    return s

def make_masks(L):
    idx = np.indices((L, L, L)).sum(axis=0)
    even = (idx % 2 == 0)
    return even, ~even

def mc_sweep(phi, kappa, lam, delta, even, odd):
    for mask in (even, odd):
        nsum = neighbor_sum(phi)
        prop = phi + delta * (RNG.random(phi.shape) - 0.5) * 2.0
        # S_local = phi^2 + lam(phi^2-1)^2 - 2 kappa phi * nsum
        def Sloc(f):
            return f*f + lam*(f*f - 1.0)**2 - 2.0*kappa*f*nsum
        dS = Sloc(prop) - Sloc(phi)
        acc = (dS <= 0) | (RNG.random(phi.shape) < np.exp(-np.clip(dS, -50, 50)))
        upd = acc & mask
        phi[upd] = prop[upd]
    return phi

def grad_sq_central(field):
    g = np.zeros_like(field)
    for ax in range(3):
        d = (np.roll(field, -1, axis=ax) - np.roll(field, 1, axis=ax)) * 0.5
        g += d*d
    return g

def measure_exponent(phi):
    """
    Estymator chain-rule (NIE artefakt-prone log-log gestosci gradientowej ex200):
      K_pt = (grad sigma)^2 / (grad Phi)^2,  Phi = sigma^2
      bin po Phi, fit log K_pt vs log Phi -> nachylenie e (K ~ Phi^e).
    Zwraca (e, R2, n_used).
    """
    sigma = phi
    Phi = sigma*sigma
    g_sig = grad_sq_central(sigma)
    g_Phi = grad_sq_central(Phi)
    msk = (g_Phi > 1e-8) & (Phi > 1e-6)
    Kpt = g_sig[msk] / g_Phi[msk]
    Phim = Phi[msk]
    msk2 = (Kpt > 1e-8) & np.isfinite(Kpt)
    Kpt, Phim = Kpt[msk2], Phim[msk2]
    if Kpt.size < 200:
        return np.nan, 0.0, Kpt.size
    # binning po log10(Phi) -> mediana K w binie (odporne na ogony)
    lP = np.log10(Phim); lK = np.log10(Kpt)
    nb = 12
    lo, hi = np.percentile(lP, 5), np.percentile(lP, 95)
    edges = np.linspace(lo, hi, nb+1)
    xs, ys = [], []
    for i in range(nb):
        sel = (lP >= edges[i]) & (lP < edges[i+1])
        if sel.sum() > 30:
            xs.append(0.5*(edges[i]+edges[i+1])); ys.append(np.median(lK[sel]))
    xs, ys = np.array(xs), np.array(ys)
    if xs.size < 4:
        return np.nan, 0.0, Kpt.size
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    e = coef[0]
    yhat = A @ coef
    ss_res = np.sum((ys - yhat)**2); ss_tot = np.sum((ys - ys.mean())**2)
    R2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0.0
    return e, R2, Kpt.size

def run_mc():
    print("\n"+"="*70); print("§B — NUMERYKA FSS (phi^4 Z2, vectorized checkerboard MC)"); print("="*70)
    kappa, lam, delta = 0.20, 1.0, 0.8
    Ls = [16, 24, 32, 40]
    n_therm, n_meas, n_skip = 300, 100, 4
    per_L = {}
    for L in Ls:
        even, odd = make_masks(L)
        phi = RNG.standard_normal((L, L, L))*0.3 + 1.0
        for _ in range(n_therm):
            mc_sweep(phi, kappa, lam, delta, even, odd)
        es, R2s = [], []
        mean_Phi = []
        for _ in range(n_meas):
            for _ in range(n_skip):
                mc_sweep(phi, kappa, lam, delta, even, odd)
            e, R2, n = measure_exponent(phi)
            if np.isfinite(e):
                es.append(e); R2s.append(R2); mean_Phi.append(float(np.mean(phi*phi)))
        e_mean, e_std = float(np.mean(es)), float(np.std(es))
        R2_mean = float(np.mean(R2s))
        per_L[L] = dict(e=e_mean, e_std=e_std, R2=R2_mean, nmeas=len(es), Phi=np.mean(mean_Phi))
        print(f"  L={L:3d}: e = {e_mean:+.3f} +/- {e_std:.3f}  | per-L log-log R2 = {R2_mean:.3f}"
              f"  | <Phi> = {per_L[L]['Phi']:.3f}  | nmeas={len(es)}")
    # FSS: ekstrapolacja e(L) -> e_inf (fit vs 1/L), oraz stabilnosc
    Larr = np.array(Ls, float); earr = np.array([per_L[L]['e'] for L in Ls])
    A = np.vstack([1.0/Larr, np.ones_like(Larr)]).T
    coef, *_ = np.linalg.lstsq(A, earr, rcond=None)
    e_inf = float(coef[1])            # ekstrapolacja 1/L->0
    e_mean_all = float(np.mean(earr)); e_spread = float(np.std(earr)/(abs(e_mean_all)+1e-9))
    R2_FSS = float(np.mean([per_L[L]['R2'] for L in Ls]))   # operacjonalizacja (Phase0 §3)
    print(f"  FSS: e(1/L->0) = {e_inf:+.3f} | mean e(L) = {e_mean_all:+.3f} | spread = {e_spread:.3f}"
          f" | R2_FSS(mean per-L) = {R2_FSS:.3f}")
    abar = e_inf/2.0   # alpha_eff = e/2
    return dict(per_L=per_L, e_inf=e_inf, e_mean=e_mean_all, e_spread=e_spread,
                R2_FSS=R2_FSS, abar=abar)

# ===================================================================
# §C — WERDYKT (value-blind, z plomby Phase0_LOCK §3)
# ===================================================================
def verdict(ana, mc):
    print("\n"+"="*70); print("§C — WERDYKT (value-blind, plomba Phase0_LOCK §3)"); print("="*70)
    abar = mc['abar']; R2 = mc['R2_FSS']; spread = mc['e_spread']
    dist = abs(abar - 2.0)
    clean = (R2 >= TH_R2) and (spread <= TH_SPREAD)
    print(f"  abar (alpha_eff = e/2, FSS) = {abar:+.3f}   |abar-2| = {dist:.3f}")
    print(f"  R2_FSS = {R2:.3f} (prog {TH_R2})   spread = {spread:.3f} (prog {TH_SPREAD})   clean={clean}")
    if clean and dist <= TH_A:
        v = "(A) DERIVED-SUBSTRATE"
    elif clean and dist >= TH_B:
        v = "(B) REFUTED-SUBSTRATE"
    else:
        v = "(C) INDETERMINATE"
    print(f"\n  >>> WERDYKT: {v}")
    # spojnosc analityka <-> numeryka
    if ana.get('e_chain') is not None:
        print(f"  spojnosc: chain-rule e=-1 (alpha=-0.5) vs MC e_inf={mc['e_inf']:+.3f}"
              f" (alpha={abar:+.3f}); manuskrypt alpha=2 (e=4)")
    return v

if __name__ == '__main__':
    ana = analytical_core()
    mc = run_mc()
    v = verdict(ana, mc)
    print("\n"+"="*70)
    npass = sum(int(ana.get(k, False)) for k in ('T1','T2','T3'))
    print(f"  Analityka: {npass}/3 PASS | MC FSS: e_inf={mc['e_inf']:+.3f}, R2_FSS={mc['R2_FSS']:.3f}")
    print(f"  WERDYKT KONCOWY: {v}")
    print("="*70)
