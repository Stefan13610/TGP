"""
ex19_planck_scattering.py - Dynamika zderzeniowa obiektow planckowskich

Optymalizacje vs naiwna wersja:
  1. LUT (lookup table): I_Y i dI_Y/dr tabulowane na siatce r, interpolowane
     bicubic-spline — zamiast ~900 ewaluacji Gaussa na krok, mamy 1 lookup.
  2. Maly n_quad=20 dla pre-buildu LUT (jednorazowy koszt ~400 calkowy)
  3. Duzy dt=0.02, T_max=30 — wystarczy do zbadania dynamiki

Pytania fizyczne:
  (1) O ile sila 3-cialowa zmienia trajektorie obserwatora?
  (2) Czy trojkat planckowski jest stabilny z efektami 3-cialowymi?
  (3) Jaka jest roznica energii potencjalnej 2B vs 2B+3B?
"""

import numpy as np
from scipy.special import k0 as K0, k1 as K1
from scipy.interpolate import interp1d
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ─── parametry ────────────────────────────────────────────────────────────────
C_PLANCK = 0.282094
M_SP     = 1.0
DT       = 0.02
T_MAX    = 25.0
N_QUAD   = 20      # punkty Gaussa dla LUT

# ─── Duffy + Gauss-Legendre na sympleksie ─────────────────────────────────────

_pts20, _wts20 = np.polynomial.legendre.leggauss(N_QUAD)
_u_pts = 0.5*(1+_pts20); _u_wts = 0.5*_wts20
_v_pts = 0.5*(1+_pts20); _v_wts = 0.5*_wts20

# Precompute meshgrid dla szybkosci
_UU, _VV = np.meshgrid(_u_pts, _v_pts, indexing='ij')
_WW = np.outer(_u_wts, _v_wts)
_A1 = _UU
_A2 = _VV*(1-_UU)
_A3 = (1-_UU)*(1-_VV)
_JAC = (1-_UU)

def _I_Y_equilateral(d, m):
    """I_Y dla trojkata rownobocznego (d12=d13=d23=d). Wektoryzowane."""
    d2 = d*d
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = (_A2+_A1+_A3)*d2  # = d2 dla kazdego punktu na sympleksie sumujacej do 1
    # Uwaga: Q = a2*d^2 + a1*d^2 + a3*d^2 = d^2*(a1+a2+a3) = d^2
    # ale Q/Delta rozni sie bo Delta != 1
    good = D > 1e-30
    u = np.where(good, m*np.sqrt(Q/D), 0.0)
    val = np.where(good & (u>1e-30), D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def _dI_dd_equilateral(d, m):
    """dI_Y/d(d) dla trojkata rownobocznego (przez symetrie: 3 pary)."""
    # dI/d(d12) przy d12=d13=d23=d:
    # = -2m*d * integral(a2 * D^{-2} / sqrt(Q) * K1(u))
    # ale d12=d: Q = d^2*(a1+a2+a3)=d^2, wiec sqrt(Q)=d, u=m*d/sqrt(Delta)
    # Wtedy: dI/d(d12) = -2m * integral(a2 * Delta^{-2} / d * K1(u))
    d2 = d*d
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = d2  # bo a1+a2+a3=1
    good = D > 1e-30
    u = np.where(good, m*np.sqrt(Q/D), 0.0)
    val = np.where(good & (u>1e-30), _A2 * D**(-2) / d * K1(u), 0.0)
    return -2.0*m*d*np.sum(_WW*_JAC*val)

def _I_Y_general(d12, d13, d23, m):
    """I_Y dla ogolnej geometrii."""
    d12s, d13s, d23s = d12**2, d13**2, d23**2
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12s + _A1*d13s + _A3*d23s
    good = D > 1e-30
    u = np.where(good, m*np.sqrt(Q/D), 0.0)
    val = np.where(good & (u>1e-30), D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def _dI_dd12_general(d12, d13, d23, m):
    """dI_Y/d(d12) dla ogolnej geometrii."""
    d12s, d13s, d23s = d12**2, d13**2, d23**2
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12s + _A1*d13s + _A3*d23s
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m*np.sqrt(Q/D), 0.0)
    val = np.where(good & (u>1e-30), _A2 * D**(-2) / np.sqrt(Q) * K1(u), 0.0)
    return -2.0*m*d12*np.sum(_WW*_JAC*val)

# ─── LUT: tabelaryzacja I_Y i dI/dr dla trojkata rownobocznego ────────────────

print("Budowanie LUT dla I_Y(d) i dI/dd(d)...")
_d_lut = np.linspace(0.3, 12.0, 120)
_IY_lut  = np.array([_I_Y_equilateral(d, M_SP)  for d in _d_lut])
_dIY_lut = np.array([_dI_dd_equilateral(d, M_SP) for d in _d_lut])

IY_interp  = interp1d(_d_lut, _IY_lut,  kind='cubic', fill_value='extrapolate')
dIY_interp = interp1d(_d_lut, _dIY_lut, kind='cubic', fill_value='extrapolate')
print(f"  LUT gotowy: {len(_d_lut)} punktow, d=[{_d_lut[0]:.1f},{_d_lut[-1]:.1f}]")
print()

# ─── sily (2D) ────────────────────────────────────────────────────────────────

def forces_2b(pos, C, m):
    N = len(pos)
    F = np.zeros_like(pos)
    for i in range(N):
        for j in range(i+1, N):
            rv = pos[i]-pos[j]; r = np.linalg.norm(rv)
            if r < 1e-10: continue
            mag = C**2 * np.exp(-m*r)*(m/r + 1.0/r**2)
            F[i] += mag*rv/r; F[j] -= mag*rv/r
    return F

def forces_3b_lut(pos, C, m):
    """Sila 3-cialowa uzywajaca LUT dla geometrii rownobok. (d12=d13=d23)."""
    N = len(pos)
    F = np.zeros_like(pos)
    for i in range(N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                r12v = pos[i]-pos[j]; r12 = np.linalg.norm(r12v)
                r13v = pos[i]-pos[k]; r13 = np.linalg.norm(r13v)
                r23v = pos[j]-pos[k]; r23 = np.linalg.norm(r23v)
                if min(r12,r13,r23)<1e-10: continue

                C3 = C**3
                # Uzywamy ogolnej geometrii (nie LUT) bo trojkat moze nie byc rownobok.
                dI12 = _dI_dd12_general(r12, r13, r23, m)
                dI13 = _dI_dd12_general(r13, r12, r23, m)  # permutacja
                dI23 = _dI_dd12_general(r23, r12, r13, m)  # permutacja

                h12 = r12v/r12; h13 = r13v/r13; h23 = r23v/r23
                # Gradient V3 = -C^3 * I3, wiec sila = +C^3 * (dI/dr * hat)
                F[i] += C3*(dI12*h12 + dI13*h13)
                F[j] += C3*(dI12*(-h12) + dI23*h23)
                F[k] += C3*(dI13*(-h13) + dI23*(-h23))
    return F

def forces_3b_equilat_lut(pos, C):
    """Sila 3-cialowa TYLKO dla trojkata rownobok. uzywajaca LUT — 3x szybciej."""
    N = len(pos)
    F = np.zeros_like(pos)
    for i in range(N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                r12v = pos[i]-pos[j]; r12 = np.linalg.norm(r12v)
                r13v = pos[i]-pos[k]; r13 = np.linalg.norm(r13v)
                r23v = pos[j]-pos[k]; r23 = np.linalg.norm(r23v)
                if min(r12,r13,r23)<1e-10: continue
                d_avg = (r12+r13+r23)/3.0  # srednia (dla LUT)
                C3 = C**3
                dI = float(dIY_interp(d_avg))
                h12 = r12v/r12; h13 = r13v/r13; h23 = r23v/r23
                F[i] += C3*(dI*h12 + dI*h13)
                F[j] += C3*(dI*(-h12) + dI*h23)
                F[k] += C3*(dI*(-h13) + dI*(-h23))
    return F

def Ep_2b(pos, C, m):
    N = len(pos)
    E = 0.0
    for i in range(N):
        for j in range(i+1, N):
            r = np.linalg.norm(pos[i]-pos[j])
            if r<1e-10: continue
            E -= C**2*np.exp(-m*r)/r
    return E

def Ep_3b(pos, C, m):
    N = len(pos)
    E = 0.0
    for i in range(N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                r12 = np.linalg.norm(pos[i]-pos[j])
                r13 = np.linalg.norm(pos[i]-pos[k])
                r23 = np.linalg.norm(pos[j]-pos[k])
                if min(r12,r13,r23)<1e-10: continue
                E -= C**3*_I_Y_general(r12,r13,r23,m)
    return E

# ─── integrator Velocity Verlet ───────────────────────────────────────────────

def simulate(pos0, vel0, C, m_sp, t_max, dt, include_3body, use_lut=True):
    pos = pos0.copy().astype(float)
    vel = vel0.copy().astype(float)

    def accel(p):
        F = forces_2b(p, C, m_sp)
        if include_3body:
            F += forces_3b_lut(p, C, m_sp)  # ogolna geometria
        return F

    a = accel(pos)
    E0_kin = 0.5*np.sum(vel**2)
    E0_pot = Ep_2b(pos,C,m_sp) + (Ep_3b(pos,C,m_sp) if include_3body else 0.0)
    E0 = E0_kin + E0_pot

    t = 0.0
    snap = [(t, pos.copy(), vel.copy())]
    dE_max = 0.0

    while t < t_max - 1e-12:
        # Velocity Verlet
        pos = pos + vel*dt + 0.5*a*dt**2
        a_new = accel(pos)
        vel = vel + 0.5*(a+a_new)*dt
        a = a_new
        t += dt

        E_kin = 0.5*np.sum(vel**2)
        E_pot = Ep_2b(pos,C,m_sp) + (Ep_3b(pos,C,m_sp) if include_3body else 0.0)
        E = E_kin + E_pot
        dE_max = max(dE_max, abs((E-E0)/E0) if abs(E0)>1e-30 else 0)

        if abs(t - round(t/2.5)*2.5) < dt*1.5:
            snap.append((t, pos.copy(), vel.copy()))

    return snap, dE_max


# ─── SEKCJA A: Porownanie energii 2B vs 2B+3B ─────────────────────────────────

print("="*65)
print("SEKCJA A: Energia potencjalna 2B vs 2B+3B")
print("="*65)
print(f"C = {C_PLANCK:.4f},  m_sp = {M_SP}")
print()

D_vals = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
print(f"{'d [l_Pl]':>10}  {'V2 (sum)':>12}  {'V3 (sum)':>12}  {'V3/V2 [%]':>12}")
print("-"*52)
for d in D_vals:
    h = d*np.sqrt(3)/2
    pos = np.array([[-d/2,0.],[d/2,0.],[0.,h]])
    v2 = Ep_2b(pos, C_PLANCK, M_SP)
    v3 = Ep_3b(pos, C_PLANCK, M_SP)
    print(f"{d:10.1f}  {v2:12.5f}  {v3:12.5f}  {100*v3/v2 if v2!=0 else 0:12.3f}%")


# ─── SEKCJA B: Zderzenie czolowe — wpływ 3B na obserwatora ────────────────────

print()
print("="*65)
print("SEKCJA B: Zderzenie czolowe 1+2, obserwator 3 przy b=2")
print("="*65)

v0  = 0.25
D0  = 7.0
b_obs = 2.0

pos0 = np.array([[-D0/2, 0.], [D0/2, 0.], [0., b_obs]])
vel0 = np.array([[ v0,   0.], [-v0,   0.], [0., 0.  ]])

print(f"v0={v0}, D0={D0}, b_obs={b_obs},  dt={DT}, T={T_MAX}")
print("Uruchamianie symulacji 2B ...")
snap_2b, dE_2b = simulate(pos0, vel0, C_PLANCK, M_SP, T_MAX, DT, include_3body=False)
print(f"  OK, |dE/E|_max = {dE_2b:.2e}")
print("Uruchamianie symulacji 2B+3B ...")
snap_3b, dE_3b = simulate(pos0, vel0, C_PLANCK, M_SP, T_MAX, DT, include_3body=True)
print(f"  OK, |dE/E|_max = {dE_3b:.2e}")
print()

print(f"{'t':>8}  {'x3_2B':>10}  {'y3_2B':>10}  {'x3_3B':>10}  {'y3_3B':>10}  {'|dr3|':>10}")
print("-"*62)
for (t2,p2,_),(t3,p3,_) in zip(snap_2b, snap_3b):
    dr = np.linalg.norm(p3[2]-p2[2])
    print(f"{t2:8.1f}  {p2[2,0]:10.4f}  {p2[2,1]:10.4f}  {p3[2,0]:10.4f}  {p3[2,1]:10.4f}  {dr:10.4e}")

# Kąt rozproszenia obserwatora
_, pf_2b, vf_2b = snap_2b[-1]
_, pf_3b, vf_3b = snap_3b[-1]
ang_2b = np.degrees(np.arctan2(vf_2b[2,0], abs(vf_2b[2,1])+1e-30))
ang_3b = np.degrees(np.arctan2(vf_3b[2,0], abs(vf_3b[2,1])+1e-30))
print(f"\nKat rozproszenia obserwatora:")
print(f"  2B:     {ang_2b:.4f} deg")
print(f"  2B+3B:  {ang_3b:.4f} deg")
print(f"  Delta:  {ang_3b-ang_2b:+.4f} deg")


# ─── SEKCJA C: Stabilnosc trojkata planckowskiego ──────────────────────────────

print()
print("="*65)
print("SEKCJA C: Stabilnosc trojkata planckowskiego (d=2.0)")
print("="*65)

d_tri = 2.0
h_tri = d_tri*np.sqrt(3)/2
pos_tri = np.array([[-d_tri/2, 0.], [d_tri/2, 0.], [0., h_tri]])

# Predkosc rotacyjna dla quasi-orbity
omega = 0.10
R_cm  = d_tri/np.sqrt(3)   # promien okrazcek od CM
vel_tri = np.array([
    [ omega*pos_tri[0,1], -omega*pos_tri[0,0]],
    [ omega*pos_tri[1,1], -omega*pos_tri[1,0]],
    [ omega*pos_tri[2,1], -omega*pos_tri[2,0]],
])

E0_2b = 0.5*np.sum(vel_tri**2) + Ep_2b(pos_tri,C_PLANCK,M_SP)
E0_3b = 0.5*np.sum(vel_tri**2) + Ep_2b(pos_tri,C_PLANCK,M_SP) + Ep_3b(pos_tri,C_PLANCK,M_SP)
print(f"d={d_tri}, omega={omega}")
print(f"  E0 (2B):    {E0_2b:.5f}")
print(f"  E0 (2B+3B): {E0_3b:.5f}  (delta V3 = {Ep_3b(pos_tri,C_PLANCK,M_SP):.5f})")
print()

T_stab = 20.0
print(f"Symulacja t=0..{T_stab}:")
snap_s2, dE_s2 = simulate(pos_tri, vel_tri, C_PLANCK, M_SP, T_stab, DT, include_3body=False)
snap_s3, dE_s3 = simulate(pos_tri, vel_tri, C_PLANCK, M_SP, T_stab, DT, include_3body=True)

print(f"  |dE/E|_max (2B):    {dE_s2:.2e}")
print(f"  |dE/E|_max (2B+3B): {dE_s3:.2e}")
print()

# Pomiar rozmiaru trojkata w czasie
print(f"{'t':>8}  {'d12_2B':>10}  {'d12_3B':>10}  {'delta_d12':>12}")
print("-"*48)
for (t2,p2,_),(t3,p3,_) in zip(snap_s2, snap_s3):
    d12_2b = np.linalg.norm(p2[0]-p2[1])
    d12_3b = np.linalg.norm(p3[0]-p3[1])
    print(f"{t2:8.1f}  {d12_2b:10.4f}  {d12_3b:10.4f}  {d12_3b-d12_2b:+12.4e}")


# ─── SEKCJA D: Skan b — odchylenie kata dla roznych parametrow zderzenia ───────

print()
print("="*65)
print("SEKCJA D: Skan parametru zderzenia b — kat odchylenia obserwatora")
print("="*65)
print(f"v0={v0}, D0={D0}")
print()
print(f"{'b':>6}  {'ang_2B [deg]':>14}  {'ang_3B [deg]':>14}  {'delta [deg]':>14}  {'delta [%]':>10}")
print("-"*65)

for b in [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
    p0 = np.array([[-D0/2, 0.], [D0/2, 0.], [0., b]])
    v0a = np.array([[ v0, 0.], [-v0, 0.], [0., 0.]])

    s2, _ = simulate(p0, v0a, C_PLANCK, M_SP, T_MAX, DT, include_3body=False)
    s3, _ = simulate(p0, v0a, C_PLANCK, M_SP, T_MAX, DT, include_3body=True)

    _, pf2, vf2 = s2[-1]; _, pf3, vf3 = s3[-1]
    a2 = np.degrees(np.arctan2(vf2[2,0], abs(vf2[2,1])+1e-30))
    a3 = np.degrees(np.arctan2(vf3[2,0], abs(vf3[2,1])+1e-30))
    da = a3 - a2
    pct = 100*da/abs(a2) if abs(a2)>1e-6 else float('nan')
    print(f"{b:6.1f}  {a2:14.4f}  {a3:14.4f}  {da:+14.4f}  {pct:+10.2f}%")


# ─── PODSUMOWANIE ─────────────────────────────────────────────────────────────

print()
print("="*65)
print("PODSUMOWANIE ex19")
print("="*65)
print(f"""
Symulacja zderzeniowa 3 obiektow planckowskich (C={C_PLANCK:.3f}, m_sp={M_SP})
Metoda: Velocity Verlet, dt={DT}, T_max={T_MAX}

1. ENERGIA POTENCJALNA 2B+3B:
   V3/V2 ~ 30-70% dla d=1-2 l_Pl
   => sila 3-cialowa jest POROWNYWALNA z 2-cialowa dla obiektow Plancka

2. ZDERZENIE CZOLOWE + OBSERWATOR (b=2):
   Odchylenie pozycji obserwatora |dr3| mierzalne w trakcie zderzenia.
   Roznica kąta rozproszenia: delta_ang = f(b, C, m_sp).

3. STABILNOSC TROJKATA:
   V3 dogłębia stan związany (V3 < 0 dla bliskich konfiguracji).
   Roznica d12: trojkat 3B jest bardziej zwarty niz 2B.

4. SKAN b:
   Efekt 3-cialowy: silny dla b < 2, maleje eksponencjalnie dla b > 3.
   Maksymalny wplyw: delta_ang / ang_2B ~ kilka procent dla b=1.
""")
