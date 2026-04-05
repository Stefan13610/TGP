#!/usr/bin/env python3
"""p115 Phase 2: Fixed-g0* robustness test."""
import sys, numpy as np, time
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, 'scripts/advanced')
from p115_mc_kink_profiles import (
    make_ode, solve_atail, run_single_realization,
    PHI, R21_PDG, R31_PDG, M_E, M_MU,
    SIGMA_EPS_V, SIGMA_EPS_K, SIGMA_DBG, SIGMA_LAM6,
    G0_MIN, G0_MAX, N_SCAN
)

np.random.seed(42)
g0_arr = np.linspace(G0_MIN, G0_MAX, N_SCAN)

print('Bazowa...')
base = run_single_realization(0, 0, 0, 0, g0_arr)
g0_e = base['g0_star']
g0_mu = PHI * g0_e
g0_tau = base.get('g0_tau', PHI**2 * g0_e)
print('g0_e=%.6f g0_mu=%.6f g0_tau=%.6f' % (g0_e, g0_mu, g0_tau))
print('r21=%.4f r31=%.2f' % (base['r21'], base['r31']))

N = 500
fix_r21, fix_r31, fix_koide = [], [], []
n_fail = 0
t0 = time.time()
for i in range(N):
    eV = np.random.normal(0, SIGMA_EPS_V)
    eK = np.random.normal(0, SIGMA_EPS_K)
    dbg = np.random.normal(0, SIGMA_DBG)
    l6 = abs(np.random.normal(0, SIGMA_LAM6))
    ode_f, ic_f = make_ode(eV, eK, dbg, l6)
    Ae = solve_atail(g0_e, ode_f, ic_f)
    Am = solve_atail(g0_mu, ode_f, ic_f)
    At = solve_atail(g0_tau, ode_f, ic_f)
    if np.isnan(Ae) or np.isnan(Am) or np.isnan(At) or Ae < 1e-12:
        n_fail += 1
        continue
    r21 = (Am / Ae)**4
    r31 = (At / Ae)**4
    mt = M_E * r31
    kn = M_E + M_MU + mt
    kd = (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(mt))**2
    fix_r21.append(r21)
    fix_r31.append(r31)
    fix_koide.append(kn / kd)
    if (i + 1) % 100 == 0:
        dt = time.time() - t0
        print('  [%d/%d] ok=%d fail=%d t=%.0fs' % (i+1, N, len(fix_r21), n_fail, dt))
        sys.stdout.flush()

fix_r21 = np.array(fix_r21)
fix_r31 = np.array(fix_r31)
fix_koide = np.array(fix_koide)
dt_total = time.time() - t0

print('')
print('Udanych: %d/%d, fail: %d, t=%.0fs' % (len(fix_r21), N, n_fail, dt_total))

print('')
print('--- r_21 (fixed g0*) ---')
print('  mean=%.4f std=%.4f spread=%.3f%%' % (
    np.mean(fix_r21), np.std(fix_r21), 100*np.std(fix_r21)/R21_PDG))
print('  delta=%+.4f%%' % (100*(np.mean(fix_r21)-R21_PDG)/R21_PDG))
print('  min/max = %.3f / %.3f' % (np.min(fix_r21), np.max(fix_r21)))

print('')
print('--- r_31 (fixed g0*) ---')
print('  mean=%.2f std=%.2f spread=%.3f%%' % (
    np.mean(fix_r31), np.std(fix_r31), 100*np.std(fix_r31)/R31_PDG))
print('  delta=%+.4f%%' % (100*(np.mean(fix_r31)-R31_PDG)/R31_PDG))
print('  min/max = %.1f / %.1f' % (np.min(fix_r31), np.max(fix_r31)))

print('')
print('--- Koide ---')
print('  mean=%.7f std=%.7f' % (np.mean(fix_koide), np.std(fix_koide)))

sp21 = 100 * np.std(fix_r21) / R21_PDG
sp31 = 100 * np.std(fix_r31) / R31_PDG
d21 = abs(100 * (np.mean(fix_r21) - R21_PDG) / R21_PDG)
d31 = abs(100 * (np.mean(fix_r31) - R31_PDG) / R31_PDG)

print('')
print('--- TESTY Fazy 2 ---')
tests = [
    ('F2-T1', 'spread(r21)=%.3f%% < 5%%' % sp21, sp21 < 5),
    ('F2-T2', 'spread(r31)=%.3f%% < 10%%' % sp31, sp31 < 10),
    ('F2-T3', '|delta(r21)|=%.4f%% < 1%%' % d21, d21 < 1),
    ('F2-T4', '|delta(r31)|=%.4f%% < 3%%' % d31, d31 < 3),
]
n_pass = 0
for tid, desc, p in tests:
    tag = 'PASS' if p else 'FAIL'
    if p:
        n_pass += 1
    print('  [%s] %s: %s' % (tag, tid, desc))

print('')
print('Wynik Fazy 2: %d/4 PASS' % n_pass)
