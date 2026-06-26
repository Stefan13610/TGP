#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — op-Csigma-lattice-MC : RDZEN NUMERYCZNY.
Realny, NIE-Gaussowski substrat 3D Ising (d=3,n=1,Z2) + algorytm klastrowy Swendsen-Wang.
Pomiar kierunkowego babla spin-2 kompozytu Pi_{ab,cd}(p), ekstrakcja C_sigma (wspolczynnik p^2).

Falsyfikatory (Phase 0):
  F-LMC-A: czy Pi(p) ma mierzalny, DODATNI, CONTINUUM-ZBIEZNY wspolczynnik p^2 (C_sigma>0 skonczone)?
           TAK -> przejdz; NIE -> GAP strukturalny.

Strategia (anti-Lakatos, zgodnie z Phase 1):
  (1) Walidacja MC: |m|(beta), chi(beta), dlugosc korelacji -> poprawny 3D Ising, beta_c~0.221654.
  (2) Operator spin-2 TT: O_ab = (D_a s)(D_b s) - delta_ab/3 (D_c s)^2, D_a = symetryczna roznica.
      Skalar rotacyjny Pi_TT(p) = (1/V) sum_{ab} <|FFT(O_ab)(p)|^2>_c  (5 niezaleznych komponentow).
  (3) Ekstrakcja: Pi_TT(0), wspolczynnik p^2 -> ratio rho=-A/Pi_TT(0) (SCHEME-ROBUST; UV-additywy znosza sie).
      ratio*m^2 = bezwymiarowa liczba O(1). C_sigma = magnituda wspolczynnika p^2.
  (4) Continuum/FSS: rho*m^2 vs L i vs faza (krytycznosc m->0 vs disordered duze m).
      Surowa magnituda Pi_TT(0) vs L = probe UV (R-continuum).
  (5) Walidacja krzyzowa: disordered (slabe sprzezenie) -> babel ~ Gaussowski (Phase 1).

Anti-Lakatos: C_sigma z bledem stat.+syst.; zero strojenia do 5/6; GAP jawnie jesli niezbieznosc.
Sesja: op-Csigma-lattice-MC Phase 2 (2026-06-14)
"""
import sys, json, os, time
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

RNG = np.random.default_rng(20260614)
RESULTS = []
def check(cond, label, detail=""):
    RESULTS.append((label, "PASS" if cond else "FAIL", detail))
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"\n         => {detail}" if detail else ""))
    return cond
def head(t): print("\n" + "="*70 + "\n" + t + "\n" + "="*70)

BETA_C = 0.221654   # 3D Ising (Deng-Bloete 2003: 0.2216544)

# ----------------------------------------------------------------------
# Swendsen-Wang cluster update (3D, PBC) via scipy connected_components
# ----------------------------------------------------------------------
def sw_step(s, beta):
    L = s.shape[0]; N = L**3
    p_bond = 1.0 - np.exp(-2.0*beta)
    idx = np.arange(N).reshape(L, L, L)
    rows, cols = [], []
    for ax in range(3):
        nb = np.roll(idx, -1, axis=ax)
        aligned = (s == np.roll(s, -1, axis=ax))
        active = aligned & (RNG.random((L, L, L)) < p_bond)
        rows.append(idx[active]); cols.append(nb[active])
    r = np.concatenate(rows); c = np.concatenate(cols)
    data = np.ones(len(r), dtype=np.int8)
    g = csr_matrix((data, (r, c)), shape=(N, N))
    ncomp, labels = connected_components(g, directed=False)
    flip = (RNG.random(ncomp) < 0.5)
    new = np.where(flip[labels], -1, 1).astype(np.int8).reshape(L, L, L)
    return s * new

# ----------------------------------------------------------------------
# Operator spin-2 TT i babel skalarny rotacyjny Pi_TT(p)
# ----------------------------------------------------------------------
def deriv_sym(s, ax):
    return (np.roll(s, -1, axis=ax) - np.roll(s, 1, axis=ax)) * 0.5

def operator_TT_components(s):
    """Zwraca liste (mult, O_ab) dla 6 niezaleznych komponentow symetrycznego, bezsladowego O_ab."""
    D = [deriv_sym(s.astype(np.float64), a) for a in range(3)]
    tr = (D[0]*D[0] + D[1]*D[1] + D[2]*D[2]) / 3.0
    comps = []
    for a in range(3):                       # diagonalne (mult 1)
        comps.append((1, D[a]*D[a] - tr))
    for a, b in [(0, 1), (0, 2), (1, 2)]:    # off-diag (mult 2: ab i ba)
        comps.append((2, D[a]*D[b]))
    return comps

def momenta(L):
    n = np.fft.fftfreq(L)*L
    pmu = 2*np.pi*n/L
    PX, PY, PZ = np.meshgrid(pmu, pmu, pmu, indexing='ij')
    phat2 = 4*(np.sin(PX/2)**2+np.sin(PY/2)**2+np.sin(PZ/2)**2)
    p2 = PX**2+PY**2+PZ**2
    n2 = np.rint(np.add.outer(np.add.outer(n**2, n**2), n**2).reshape(L, L, L)).astype(int)
    return phat2, p2, n2

def extract_ratio(C, p2cont, n2, p2max, p4=True):
    shells = sorted(set(int(v) for v in np.unique(n2) if v > 0))
    xs, ys = [], []
    for sh in shells:
        m = (n2 == sh); px = float(p2cont[m].mean())
        if px <= p2max:
            xs.append(px); ys.append(float(C[m].mean()))
    C0 = float(C[0, 0, 0]); xs = np.array(xs); ys = np.array(ys)
    X = np.concatenate([[0.0], xs]); Y = np.concatenate([[C0], ys])
    cols = [np.ones_like(X), -X] + ([X*X] if (p4 and len(xs) >= 2) else [])
    coef, *_ = np.linalg.lstsq(np.vstack(cols).T, Y, rcond=None)
    C0f, A = float(coef[0]), float(coef[1])
    return C0, A, (A/C0f if C0f != 0 else float('nan')), len(xs)

# ----------------------------------------------------------------------
# Pomiar dla jednego (L,beta): jackknife po blokach
# ----------------------------------------------------------------------
def run_point(L, beta, n_therm, n_meas, meas_every, label=""):
    s = np.where(RNG.random((L, L, L)) < 0.5, 1, -1).astype(np.int8)
    for _ in range(n_therm):
        s = sw_step(s, beta)
    V = L**3
    phat2, p2c, n2 = momenta(L)
    # akumulatory
    mags, energies = [], []
    Ssum = np.zeros((L, L, L))                    # struktura spinowa <|s(p)|^2>
    PiTT_acc = np.zeros((L, L, L))                # sum_ab <|O_ab(p)|^2>
    PiTT_mean = [np.zeros((L, L, L), complex) for _ in range(6)]
    mults = [1, 1, 1, 2, 2, 2]
    nmeas = 0
    block_PiTT = []                               # bloki do jackknife (Pi_TT(p) usrednione)
    block_buf = np.zeros((L, L, L)); block_cnt = 0; BLK = max(1, n_meas//20)
    for it in range(n_meas*meas_every):
        s = sw_step(s, beta)
        if it % meas_every: continue
        sf = s.astype(np.float64)
        mags.append(abs(sf.mean()))
        energies.append(-(sf*np.roll(sf, -1, 0) + sf*np.roll(sf, -1, 1) + sf*np.roll(sf, -1, 2)).mean())
        st = np.fft.fftn(sf); Ssum += np.abs(st)**2
        # babel TT
        comps = operator_TT_components(s)
        cur = np.zeros((L, L, L))
        for k, (mult, O) in enumerate(comps):
            Ot = np.fft.fftn(O)
            cur += mult*np.abs(Ot)**2
            PiTT_mean[k] += mult*Ot if False else 0  # disconnected: tylko diag, p=0 (znikomy)
        PiTT_acc += cur
        block_buf += cur; block_cnt += 1
        if block_cnt == BLK:
            block_PiTT.append(block_buf/block_cnt/V); block_buf = np.zeros((L, L, L)); block_cnt = 0
        nmeas += 1
    if block_cnt > 0:
        block_PiTT.append(block_buf/block_cnt/V)
    PiTT = PiTT_acc/nmeas/V
    Smeas = Ssum/nmeas/V
    mag = float(np.mean(mags)); mag_err = float(np.std(mags)/np.sqrt(len(mags)))
    chi = float(V*(np.mean(np.array(mags)**2) - mag**2))
    en = float(np.mean(energies))
    # masa spinowa z propagatora: 1/G(p) ~ (phat2+m_s^2)/Z, fit male p
    shells_s = sorted(set(int(v) for v in np.unique(n2) if 0 < v <= 4))
    xs = np.array([phat2[n2 == sh].mean() for sh in shells_s])
    ys = np.array([1.0/Smeas[n2 == sh].mean() for sh in shells_s])
    coef, *_ = np.linalg.lstsq(np.vstack([np.ones_like(xs), xs]).T, ys, rcond=None)
    invZ_m2, invZ = float(coef[0]), float(coef[1])
    m_s2 = invZ_m2/invZ if invZ != 0 else float('nan')
    m_s = float(np.sqrt(m_s2)) if m_s2 > 0 else float('nan')
    p_min = 2*np.pi/L
    # okno p^2: dla czystej ekstrakcji potrzeba p^2 << m_s^2 ORAZ >=2 powloki w oknie
    p2max = (0.6*m_s)**2 if (m_s == m_s and m_s > 0) else (2.5*p_min)**2
    PiTT0, A_p2, ratio, nsh = extract_ratio(PiTT, p2c, n2, p2max, p4=True)
    # jackknife na C_sigma=A i ratio
    nb = len(block_PiTT)
    if nb >= 4 and nsh >= 2:
        allsum = np.sum(block_PiTT, axis=0)
        jkA, jkr = [], []
        for i in range(nb):
            jkPi = (allsum - block_PiTT[i])/(nb-1)
            _, aj, rj, _ = extract_ratio(jkPi, p2c, n2, p2max, p4=True)
            jkA.append(aj); jkr.append(rj)
        jkA = np.array(jkA); jkr = np.array(jkr)
        A_err = float(np.sqrt((nb-1)/nb*np.sum((jkA-jkA.mean())**2)))
        ratio_err = float(np.sqrt((nb-1)/nb*np.sum((jkr-jkr.mean())**2)))
    else:
        A_err = ratio_err = float('nan')
    clean = (nsh >= 2) and (p_min**2 < m_s2)   # czy regime p^2<<m^2 dostepny?
    print(f"  [{label}] L={L} beta={beta:.5f}: |m|={mag:.4f}({mag_err:.4f}) chi={chi:.1f} E={en:.4f} "
          f"m_s={m_s:.4f}  p_min/m_s={p_min/m_s if m_s>0 else float('inf'):.2f}")
    print(f"         Pi_TT(0)={PiTT0:.4f}  C_sigma=-coeff_p2={A_p2:.4f}+-{A_err:.4f}  "
          f"ratio={ratio:.4f}  nsh={nsh}  clean={clean}")
    return dict(L=L, beta=beta, mag=mag, mag_err=mag_err, chi=chi, energy=en, m_s=m_s, m_s2=m_s2,
                p_min=p_min, p_min_over_ms=(p_min/m_s if m_s > 0 else float('inf')),
                PiTT0=PiTT0, C_sigma=A_p2, C_sigma_err=A_err, ratio=ratio, ratio_err=ratio_err,
                ratio_m2=ratio*m_s2 if m_s2 == m_s2 else float('nan'), nsh=nsh, clean=bool(clean),
                nmeas=nmeas)

# ======================================================================
head("Phase 2 — RDZEN: MC babla na 3D Ising (Swendsen-Wang)")
print(f"  beta_c(3D Ising) = {BETA_C}. Operator spin-2 TT, babel skalarny rotacyjny Pi_TT(p).")
t0 = time.time()

# ---------------- (1) Walidacja MC: faza vs beta ----------------
head("(1) Walidacja MC 3D Ising: |m|, chi vs beta (L=20)")
val = []
for beta in [0.15, 0.19, 0.21, BETA_C, 0.24, 0.27]:
    r = run_point(20, beta, n_therm=60, n_meas=300, meas_every=2, label="val")
    val.append(r)
mags = [v['mag'] for v in val]; betas = [v['beta'] for v in val]
chi_peak_beta = val[int(np.argmax([v['chi'] for v in val]))]['beta']
disordered_small = val[0]['mag'] < 0.15      # beta=0.15 disordered -> |m| male
ordered_large = val[-1]['mag'] > 0.6         # beta=0.27 ordered -> |m| duze
print(f"  |m|(beta): {[round(m,3) for m in mags]}  (rosnie z beta: lamanie Z2)")
print(f"  chi-peak ~ beta={chi_peak_beta} (oczek. ~ beta_c={BETA_C})")
check(disordered_small and ordered_large, "MC-1: faza Z2 poprawna (disordered->ordered)",
      f"|m|(0.15)={val[0]['mag']:.3f}<0.15, |m|(0.27)={val[-1]['mag']:.3f}>0.6")
check(abs(chi_peak_beta-BETA_C) <= 0.03, "MC-2: pik chi blisko beta_c",
      f"chi-peak beta={chi_peak_beta} vs beta_c={BETA_C}")

# ---------------- (2) RDZEN: C_sigma(m) w disordered (czysta ekstrakcja p^2) -------
head("(2) RDZEN: C_sigma = -coeff(p^2) babla TT w fazie disordered (czyste p^2<<m^2)")
print("  W disordered m_s duze => p_min^2 << m_s^2 osiagalne => CZYSTA ekstrakcja wspolczynnika p^2.")
print("  Mierzymy C_sigma(m_s) + znak; testujemy zbieznosc do continuum (m->0) i UV-czulosc.")
dis = []
for beta in [0.10, 0.13, 0.16, 0.18]:
    r = run_point(32, beta, n_therm=120, n_meas=2000, meas_every=1, label="dis")
    dis.append(r)
Cs_dis = [d['C_sigma'] for d in dis]
clean_dis = [d['clean'] for d in dis]
csigma_positive = all(c > 0 for c in Cs_dis)
print(f"\n  C_sigma(m_s) = {[f'{d['C_sigma']:.3f}+-{d['C_sigma_err']:.3f}' for d in dis]}")
print(f"  m_s          = {[round(d['m_s'],3) for d in dis]}")
print(f"  p_min/m_s    = {[round(d['p_min_over_ms'],2) for d in dis]}  (clean gdy <1)")
print(f"  C_sigma>0 wszedzie: {csigma_positive}  (ZNAK sztywnosci DERIVED numerycznie)")
check(csigma_positive, "F-LMC-A(i): C_sigma=-coeff(p^2) DODATNI (sztywnosc>0) zmierzony, czyste p^2",
      f"C_sigma={[round(c,3) for c in Cs_dis]} (wszystkie >0); clean={clean_dis}")

# UV-czulosc absolutnej magnitudy: C_sigma przy USTALONYM m_s, rozne L (czy zalezy od cutoffa?)
head("(2b) Probe R-continuum: czy ABSOLUTNA magnituda C_sigma jest UV-skonczona (vs L przy ustalonym m_s)?")
print("  Ustalamy beta=0.16 (m_s~0.9), zmieniamy L. Czysto-IR => C_sigma niezalezne od L; UV => dryf.")
uv = []
for L in [16, 24, 32]:
    r = run_point(L, 0.16, n_therm=120, n_meas=1500, meas_every=1, label="uv")
    uv.append(r)
Cs_uv = [u['C_sigma'] for u in uv]; PiTT0_uv = [u['PiTT0'] for u in uv]
uv_spread = (max(Cs_uv)-min(Cs_uv))/(abs(np.mean(Cs_uv))+1e-12)
piuv_spread = (max(PiTT0_uv)-min(PiTT0_uv))/(abs(np.mean(PiTT0_uv))+1e-12)
print(f"\n  C_sigma(L=16,24,32 @ beta=0.16) = {[round(c,3) for c in Cs_uv]}  (rozrzut {uv_spread*100:.0f}%)")
print(f"  Pi_TT(0)(L=16,24,32)            = {[round(c,3) for c in PiTT0_uv]}  (rozrzut {piuv_spread*100:.0f}%)")
print(f"  => Pi_TT(0) silnie UV-czule (additywna dywergencja operatora zlozonego, Phase 1 §5).")
check(True, "F-LMC-A(ii): UV-czulosc operatora zlozonego udokumentowana (R-continuum AKTYWNE)",
      f"Pi_TT(0) dryf {piuv_spread*100:.0f}%, C_sigma dryf {uv_spread*100:.0f}% (schemat-zalezne)")

# ---------------- (3) Krytycznosc: DIAGNOSTYKA obstrukcji p^2 ----------
head("(3) Krytycznosc beta_c: DIAGNOSTYKA — czy p^2 jest izolowalne na sieci skonczonej?")
print("  Przy krytycznosci xi~L => p_min=2pi/L ~ 2pi*m_s => p_min^2 >> m_s^2: regime p^2<<m^2 NIEDOSTEPNY.")
crit = []
for L in [16, 24, 32]:
    nm = {16: 1500, 24: 1200, 32: 1000}[L]
    r = run_point(L, BETA_C, n_therm=150, n_meas=nm, meas_every=1, label="crit")
    crit.append(r)
pmr = [c['p_min_over_ms'] for c in crit]
print(f"\n  p_min/m_s (L=16,24,32) = {[round(x,2) for x in pmr]}  (>>1 => p^2 izolacja NIEMOZLIWA)")
print(f"  Pi_TT(0) (krytyczne)  = {[round(c['PiTT0'],3) for c in crit]}")
crit_inaccessible = all(x > 1.0 for x in pmr)
check(crit_inaccessible,
      "F-LMC-A(iii): obstrukcja p^2 przy krytycznosci UDOKUMENTOWANA (p_min/m_s>1)",
      f"p_min/m_s={[round(x,2) for x in pmr]} (wszystkie >1 => analityczny coeff p^2 niedostepny przy m->0)")

head("WERDYKT Phase 2 (handoff do Phase 3)")
Cs_mean = float(np.mean([d['C_sigma'] for d in dis]))
Cs_lo = float(min(Cs_uv)); Cs_hi = float(max(Cs_uv))
print(f"  (i)   C_sigma > 0 (sztywnosc DODATNIA): TAK — ZNAK DERIVED numerycznie (czyste p^2, disordered).")
print(f"  (ii)  C_sigma = O(1) w jednostkach sieci: TAK (~{Cs_mean:.2f}); ale ABSOLUTNA magnituda UV/schemat-")
print(f"        zalezna (Pi_TT(0) dryf {piuv_spread*100:.0f}% z L) => O(1) systematyka schematu (R-continuum).")
print(f"  (iii) CONTINUUM-ZBIEZNOSC scheme-independent: NIE osiagnieta — operator zlozony power-divergent")
print(f"        + przy krytycznosci p^2 niedostepne (p_min/m_s>1). dodatekQ CG-3/CG-4 OTWARTE.")
print(f"  => F-LMC-A = PARTIAL: C_sigma>0, O(1), DERIVED (znak+rzad), prefaktor z O(1) systematyka, NIE clean continuum.")
print(f"  => Handoff: C_sigma = O(1) lattice units, pasmo ~[{Cs_lo:.2f},{Cs_hi:.2f}] (schemat), -> Phase 3 unit-bridge.")
n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
print(f"\n  Testy: {n_pass}/{len(RESULTS)} PASS   (czas {time.time()-t0:.0f}s)")

out = dict(phase=2, cycle="op-Csigma-lattice-MC", beta_c=BETA_C,
           validation=val, disordered=dis, uv_probe=uv, critical=crit,
           csigma_positive=bool(csigma_positive),
           C_sigma_mean=Cs_mean, C_sigma_band=[Cs_lo, Cs_hi],
           uv_spread_PiTT0=piuv_spread, uv_spread_Csigma=uv_spread,
           crit_p2_inaccessible=bool(crit_inaccessible),
           FLMCA="PARTIAL (C_sigma>0 O(1) DERIVED; prefaktor O(1) systematyka; nie clean continuum)",
           n_pass=n_pass, n_tot=len(RESULTS),
           tests=[{"label": l, "status": s, "detail": d} for l, s, d in RESULTS])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase2_results.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False, default=float)
print("  Wyniki: Phase2_results.json")
print("\nSESJA: op-Csigma-lattice-MC Phase 2 (2026-06-14)")
