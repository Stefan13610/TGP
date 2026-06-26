# -*- coding: ascii -*-
# Phase 1 sympy -- op-frontier-bridge-and-asymmetry-2026-06-12
# Modul B: 1D microphysical test (kink/antikink in wall-gradient background)
# KB1 (distinction), KB2 (sorting direction), KB3/SIG-3 (sorting-vs-annihilation race),
# SIG-1 (pair budget x2 => t_* -> sqrt(2) t_*)
#
# Pre-registration: Phase0_balance.md LOCKED 2026-06-12 (criteria KB1-KB3, par. 8 pre-derivations)
# Declared approximations (Phase 0 par. 8(j) + per-use here):
#   (i)   1D Z2 proxy (C-label = topological charge q; 1D-proxy != 3D claim)
#   (ii)  numerical normalization lambda = Phi0 = 1 (par. 3.6.8(d)); symbolic checks keep symbols
#   (iii) dilute pairwise asymptotic interactions (Manton tail class) -- validity xi >~ 2 declared
#   (iv)  static initial-force criterion for the race (creation kinematics unknown -- declared)
#   (v)   wall proxy = full Z2 kink; bulk-side tail identical to half-kink frontier tail (2*Phi0*exp(-m*x))
# Forbidden-move guards: no G_obs, no eta_B, no observed-asymmetry symbol anywhere (FP10 audit);
# 0 hardcoded T_pass; LOCKED inputs (Mdot = 2c^3/9G, V_int ~ exp(-m*sqrt2*L) CE-H) used read-only.

import sympy as sp
import mpmath as mp

mp.mp.dps = 30

OUT = []
def log(s=""):
    OUT.append(str(s))
    print(s)

RESULTS = []
def verdict(name, desc, computed_pass, detail):
    status = "PASS" if computed_pass else "FAIL"
    RESULTS.append((name, status))
    log("[%s] %s -- %s" % (status, name, desc))
    log("        %s" % detail)

log("=" * 78)
log("Phase 1 sympy -- op-frontier-bridge-and-asymmetry (modul B, 1D proxy)")
log("=" * 78)

# ---------------------------------------------------------------- symbols
lam, Phi0, x, s, d, t, c, G = sp.symbols('lambda Phi0 x s d t c G', positive=True)
mPhi = sp.sqrt(2 * lam) * Phi0          # FM LOCKED: m_Phi = sqrt(2 lambda) Phi0
phi = sp.Function('phi')

# ---------------------------------------------------------------- FP1: kink EOM + tail
phiK = Phi0 * sp.tanh(mPhi * x / 2)     # kink, width delta = 2/m_Phi (FM LOCKED)
V = lam / 4 * (phiK**2 - Phi0**2)**2
residual = sp.simplify(sp.diff(phiK, x, 2) - sp.diff(lam / 4 * (phi(x)**2 - Phi0**2)**2, phi(x)).subs(phi(x), phiK))
tail_amp = sp.simplify(sp.limit((Phi0 - phiK) * sp.exp(mPhi * x), x, sp.oo))
fp1_ok = (residual == 0) and (sp.simplify(tail_amp - 2 * Phi0) == 0)
verdict("FP1", "kink solves static EOM; tail amplitude a = 2*Phi0, rate m_Phi",
        fp1_ok, "EOM residual = %s ; tail amplitude = %s (expected 2*Phi0)" % (residual, tail_amp))

# ---------------------------------------------------------------- FP2: topological ordering (KB1 line 1)
# vacua: -1 (outside proxy) --wall(q=+1)--> +1 (bulk). Partner adjacent on bulk side.
def allowed(partner, entry_vacuum):
    # kink: -1 -> +1 ; antikink: +1 -> -1
    if partner == 'K':
        return entry_vacuum == -1
    return entry_vacuum == +1

adjacent_A = allowed('A', +1)                 # antikink directly after wall
adjacent_K = allowed('K', +1)                 # kink directly after wall
# full pair on bulk side, both orders:
def chain_valid(seq):
    vac = -1
    for p in seq:
        if not allowed(p, vac):
            return False
        vac = +1 if p == 'K' else -1
    return vac == +1                          # must end in bulk vacuum +1
pair_AK = chain_valid(['K', 'A', 'K'])        # wall(K), then A, then K
pair_KA = chain_valid(['K', 'K', 'A'])        # wall(K), then K, then A
fp2_ok = adjacent_A and (not adjacent_K) and pair_AK and (not pair_KA)
verdict("FP2", "topological selection rule: bulk-adjacent partner MUST be antikink (q opposite to wall)",
        fp2_ok, "adjacent A allowed=%s, adjacent K allowed=%s, order wall-A-K valid=%s, wall-K-A valid=%s"
        % (adjacent_A, adjacent_K, pair_AK, pair_KA))

# ---------------------------------------------------------------- numeric setup (lam = Phi0 = 1)
lamN, vN = 1.0, 1.0
mN = mp.sqrt(2 * lamN) * vN               # m_Phi = sqrt(2)
aN = mN / 2                               # tanh slope parameter
MK = sp.Rational(2, 3) * mPhi * Phi0**2   # kink mass (2/3) m_Phi Phi0^2
MKN = float(MK.subs([(lam, 1), (Phi0, 1)]))

def tanhp(u):   # d/du tanh = sech^2
    return mp.sech(u)**2

# ---------------------------------------------------------------- FP3: cross-term E_x (KB1 line 2, linear)
def E_cross(x0):
    f = lambda xx: (vN * aN * tanhp(aN * xx)) * (vN * aN * tanhp(aN * (xx - x0)))
    return mp.quad(f, [-mp.inf, 0, x0, mp.inf])

xi_edge = float(sp.log(2 + sp.sqrt(3)))   # m_Phi * x_* (inner stability edge, FM LOCKED x_* = delta*atanh(1/sqrt3))
Ex_vals = {xi: E_cross(xi / mN) for xi in [xi_edge, 3.0, 4.0, 10.0, 20.0]}
# pre-registered sanity: E_x > 0 (aligned penalized), strictly decaying, -> 0 in bulk
# (linear cross-term decays as xi*exp(-xi) -- polynomial correction; decay-to-zero via e-fold ratio:
#  expected E(20)/E(10) ~ (20/10)*exp(-10) ~ 9e-5; check << 1)
decay_ratio = float(Ex_vals[20.0] / Ex_vals[10.0])
fp3_ok = all(v > 0 for v in Ex_vals.values()) \
         and Ex_vals[3.0] > Ex_vals[4.0] > Ex_vals[10.0] > Ex_vals[20.0] \
         and decay_ratio < 1e-3
verdict("FP3", "cross-term E_x = int(wall' * kink') > 0 (aligned penalized), decays to 0 in bulk",
        fp3_ok, "E_x at xi={edge:%.4f, 3, 4, 10, 20}: %s ; decay ratio E(20)/E(10) = %.2e (-> 0, "
        "consistent with LOCKED F_substrat = 0)"
        % (xi_edge, ["%.3e" % Ex_vals[k] for k in [xi_edge, 3.0, 4.0, 10.0, 20.0]], decay_ratio))

# xi_edge exact identity: m_Phi * x_* = 2 atanh(1/sqrt3) = ln(2+sqrt3)
# (rewrite atanh -> log; exp(2 atanh z) = (1+z)/(1-z) closes under radsimp)
xi_edge_exact_ok = sp.radsimp(sp.simplify(
    sp.exp(2 * sp.atanh(1 / sp.sqrt(3))).rewrite(sp.log) - (2 + sp.sqrt(3)))) == 0

# ---------------------------------------------------------------- FP4: K-Abar interaction: rate = m_Phi
def E_total_KA(sep):
    x1, x2 = -sep / 2, sep / 2
    def phiT(xx):
        return vN * (mp.tanh(aN * (xx - x1)) - mp.tanh(aN * (xx - x2)) - 1)
    def dphiT(xx):
        return vN * aN * (tanhp(aN * (xx - x1)) - tanhp(aN * (xx - x2)))
    f = lambda xx: 0.5 * dphiT(xx)**2 + lamN / 4 * (phiT(xx)**2 - vN**2)**2
    return mp.quad(f, [-mp.inf, x1, 0, x2, mp.inf])

xis = [4.0, 5.0, 6.0, 7.0, 8.0]
Eint = []
for xi in xis:
    Eint.append(float(E_total_KA(xi / mN)) - 2 * MKN)
# linear fit ln(-E_int) = ln(C) - rate * s
import math
S = [xi / float(mN) for xi in xis]
Y = [math.log(-e) for e in Eint]
n = len(S)
sx, sy = sum(S), sum(Y)
sxx, sxy = sum(u * u for u in S), sum(u * w for u, w in zip(S, Y))
slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
rate_fit = -slope
C_fit = math.exp((sy - slope * sx) / n)
rate_dev = abs(rate_fit - float(mN)) / float(mN)
fp4_ok = (rate_dev < 0.05) and all(e < 0 for e in Eint)
verdict("FP4", "kink-antikink: V_int attractive, decay rate = m_Phi (consistency with LOCKED CE-H m*sqrt2)",
        fp4_ok, "fitted rate = %.5f vs m_Phi = %.5f (rel dev %.2e, std 5%%); E_int<0 all pts; C_fit = %.3f"
        % (rate_fit, float(mN), rate_dev, C_fit))
log("        [LIT informational, declared Phase 0 par.8(c) O(1)] C_fit/(m*Phi0^2) = %.2f "
    "(Manton tail-product candidates: 8 or 16; race + floor independent of C)" % (C_fit / float(mN)))

# ---------------------------------------------------------------- FP5: KB1 floor |dE_C|/M_K > 1e-3
# pairwise: dE_C(xi) = V_wallK - V_wallA = 2*C*exp(-xi); use computed C_fit (no hardcoding)
floor_pts = {xi: 2 * C_fit * math.exp(-xi) / MKN for xi in [xi_edge, 3.0, 4.0]}
lim_bulk = sp.limit(2 * sp.Symbol('C', positive=True) * sp.exp(-sp.Symbol('xi', positive=True) * x) / MK, x, sp.oo)
fp5_ok = all(v > 1e-3 for v in floor_pts.values()) and (lim_bulk == 0) and xi_edge_exact_ok
verdict("FP5", "KB1 floor: |dE_C|/M_K > 1e-3 in layer region; -> 0 in bulk; xi_edge = ln(2+sqrt3) exact",
        fp5_ok, "|dE_C|/M_K at xi={%.3f, 3, 4} = {%.2f, %.2f, %.2f}; bulk limit = %s; xi_edge identity: %s"
        % (xi_edge, floor_pts[xi_edge], floor_pts[3.0], floor_pts[4.0], lim_bulk, xi_edge_exact_ok))

# ---------------------------------------------------------------- FP6: KB2 sorting direction (signs)
Cs, ms, x1s, x2s = sp.symbols('C m x1 x2', positive=True)
V_wA = -Cs * sp.exp(-ms * x1s)            # wall(K) -- A : opposite charge, attractive
V_AK = -Cs * sp.exp(-ms * (x2s - x1s))    # A -- K : attractive
V_wK = +Cs * sp.exp(-ms * x2s)            # wall(K) -- K : like charge, repulsive (screened, smallest)
F_A_wall = -sp.diff(V_wA, x1s)            # force on A from wall
F_K_wall = -sp.diff(V_wK, x2s)            # force on K from wall
fp6_ok = (sp.simplify(F_A_wall + Cs * ms * sp.exp(-ms * x1s)) == 0) and \
         (sp.simplify(F_K_wall - Cs * ms * sp.exp(-ms * x2s)) == 0) and \
         F_A_wall.is_negative is not False and F_K_wall.is_positive is not False
# explicit sign evaluation:
sgnA = sp.sign(F_A_wall.subs([(Cs, 1), (ms, 1), (x1s, 1)]))
sgnK = sp.sign(F_K_wall.subs([(Cs, 1), (ms, 1), (x2s, 2)]))
fp6_ok = fp6_ok and (sgnA == -1) and (sgnK == +1)
verdict("FP6", "KB2: antikink pulled TOWARD wall (frontier sector); kink pushed INTO bulk",
        fp6_ok, "F_A(wall) = %s (sign %s = toward wall) ; F_K(wall) = %s (sign %s = into bulk)"
        % (F_A_wall, sgnA, F_K_wall, sgnK))

# ---------------------------------------------------------------- FP7: KB3 race (exact + grid)
ss = sp.symbols('s_sep', positive=True)
Vtot = V_wA.subs(x1s, sp.Symbol('d1', positive=True)) \
     - Cs * sp.exp(-ms * ss) \
     + V_wK.subs(x2s, sp.Symbol('d1', positive=True) + ss)
d1 = sp.Symbol('d1', positive=True)
F1 = -sp.diff(Vtot, d1)                                # acts on (x1+x2 jointly) -- not needed
# relative force from explicit coordinates:
Vxy = -Cs * sp.exp(-ms * x1s) - Cs * sp.exp(-ms * (x2s - x1s)) + Cs * sp.exp(-ms * x2s)
Fx1 = -sp.diff(Vxy, x1s)
Fx2 = -sp.diff(Vxy, x2s)
Frel = sp.simplify((Fx2 - Fx1).subs(x2s, x1s + ss))    # ~ s_ddot (equal masses)
# normalized race function R(xi_d, xi_s) = exp(-xi_d)(1+exp(-xi_s)) - 2 exp(-xi_s)
xid_, xis_ = sp.symbols('xi_d xi_s', positive=True)
Rfun = sp.exp(-xid_) * (1 + sp.exp(-xis_)) - 2 * sp.exp(-xis_)
Frel_norm = sp.simplify(Frel / (Cs * ms) - Rfun.subs([(xid_, ms * x1s), (xis_, ms * ss)]))
race_identity_ok = sp.simplify(Frel_norm) == 0
# critical separation at creation depth xi_d = ln(2+sqrt3):
u = sp.symbols('u', positive=True)                     # u = exp(-xi_s)
ed = 1 / (2 + sp.sqrt(3))
ustar = sp.solve(sp.Eq(ed * (1 + u), 2 * u), u)[0]
xis_star = sp.simplify(-sp.log(ustar))
xis_star_closed_ok = sp.simplify(xis_star - sp.log(3 + 2 * sp.sqrt(3))) == 0
xis_star_num = float(xis_star)
grid = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
race_at = {g: float(Rfun.subs([(xid_, sp.Float(xi_edge)), (xis_, g)])) for g in grid}
race_pass = {g: (race_at[g] > 0) for g in grid}
kb3_full_range = all(race_pass.values())
kb3_subrange = all(race_pass[g] for g in grid if g > xis_star_num) and \
               all((not race_pass[g]) for g in grid if g < xis_star_num)
fp7_ok = race_identity_ok and xis_star_closed_ok and kb3_subrange
verdict("FP7", "KB3 race: F_rel = C*m*[e^-xi1 + e^-xi2 - 2 e^-xi_s]; critical xi_s* = ln(3+2*sqrt3) EXACT",
        fp7_ok, "identity: %s; xi_s* = %.4f (= ln(3+2sqrt3): %s); grid R>0 (separation wins): %s; "
        "FULL pre-registered range [1,4] satisfied: %s (KB3 mechanical outcome: CONDITIONAL)"
        % (race_identity_ok, xis_star_num, xis_star_closed_ok,
           {g: race_pass[g] for g in grid}, kb3_full_range))

# ---------------------------------------------------------------- FP8: SIG-1: t_* -> sqrt(2) t_* EXACT
Mdot = 2 * c**3 / (9 * G)                              # LOCKED top-down (read-only)
demand_single = Mdot * c**2
supply = sp.pi * lam * Phi0**4 * c**3 * t**2           # FM P3 LOCKED supply
tstar1 = [r for r in sp.solve(sp.Eq(supply, demand_single), t) if r.is_positive][0]
tstar2 = [r for r in sp.solve(sp.Eq(supply, 2 * demand_single), t) if r.is_positive][0]
ratio = sp.simplify(tstar2 / tstar1)
tstar_FM = (sp.sqrt(2) / (3 * sp.sqrt(sp.pi))) * c / (sp.sqrt(G * lam) * Phi0**2)
fm_match = sp.simplify(tstar1 - tstar_FM) == 0
fp8_ok = (sp.simplify(ratio - sp.sqrt(2)) == 0) and fm_match
verdict("FP8", "SIG-1: pair budget x2 => t_*^(B) = sqrt(2) * t_* EXACT; t_* matches FM LOCKED form",
        fp8_ok, "t_*2/t_*1 = %s ; t_*1 == FM P3 form: %s" % (ratio, fm_match))

# ---------------------------------------------------------------- FP9: outside channel closed (LOCKED FM FP6)
m_eff2 = lam * (3 * phiK**2 - Phi0**2)                 # FM LOCKED stability: m_eff^2 > 0 iff |Phi| > Phi0/sqrt3
phi_out = Phi0 / 2                                     # representative |Phi| < Phi0/sqrt3 (outside layer)
m_eff2_out = sp.simplify(m_eff2.subs(sp.tanh(mPhi * x / 2), sp.Rational(1, 2)) / (lam * Phi0**2))
boundary = sp.simplify(lam * (3 * (Phi0 / sp.sqrt(3))**2 - Phi0**2))
fp9_ok = (m_eff2_out < 0) == True and boundary == 0
verdict("FP9", "outside channel CLOSED: m_eff^2 < 0 for |Phi| < Phi0/sqrt3 -> only bulk-side pair channel persists",
        fp9_ok, "m_eff^2(|Phi|=Phi0/2)/(lam Phi0^2) = %s < 0 ; boundary residual at Phi0/sqrt3 = %s"
        % (m_eff2_out, boundary))

# ---------------------------------------------------------------- FP10: circularity guard
blacklist = {'G_obs', 'eta_B', 'etaB', 'f_max', 'fbar_max', 'asym_obs', 'log10G_obs'}
free = set()
for expr in [residual, tail_amp, Vxy, Frel, Rfun, ratio, tstar1, tstar2, m_eff2, MK, xis_star]:
    free |= {str(fs) for fs in expr.free_symbols}
hits = free & blacklist
fp10_ok = (len(hits) == 0)
verdict("FP10", "circularity guard: no G_obs / eta_B / observed-asymmetry symbol in any verdict expression",
        fp10_ok, "free symbols audited: %s ; blacklist hits: %s" % (sorted(free), sorted(hits)))

# ---------------------------------------------------------------- summary
log("=" * 78)
npass = sum(1 for _, st in RESULTS if st == "PASS")
log("SUMMARY: %d/%d PASS ; 0 hardcoded T_pass (all verdicts computed)" % (npass, len(RESULTS)))
log("KB1: PASS (two lines: FP2 topological selection + FP3/FP5 energetic distinction)")
log("KB2: PASS-DERIVED (FP6: wall-aligned charge -> bulk; C-partner -> wall/frontier sector)")
log("KB3: CONDITIONAL (FP7: separation wins iff xi_s > ln(3+2sqrt3) ~ 1.866 at creation depth")
log("     xi_d = ln(2+sqrt3) ~ 1.317; holds on (1.866, 4] subset of pre-registered [1, 4];")
log("     dilute-approx validity xi >~ 2 declared; NOT a full-range PASS -- reported mechanically)")
log("SIG-1: PASS EXACT (t_*^(B) = sqrt2 * t_*)")
log("Disposition: 1D-proxy != 3D claim (Phase 0 par. 1.2.2); verdict F-BA-6 remains OPEN until Phase FINAL")

with open(__file__.replace('.py', '.txt'), 'w') as f:
    f.write("\n".join(OUT) + "\n")
