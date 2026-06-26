# -*- coding: utf-8 -*-
# op-PSR-spin-substrate-channel-2026-06-13 / correlation signs (EXPLORATORY)
# Author's challenge (2026-06-13): the source contributions should NOT just add;
# when you measure correlations between sources you get SUBTRACTION, not addition.
#
# Test: compute the SIGN of the cross-correlation 2<X_1 : X_2> for each multipole
# (monopole, dipole, quadrupole) of the J0737 two-body system in the CM frame.
# Does correlation subtract (destructive) or add (constructive), and for WHICH
# multipole? Is the subtraction the author senses already used by TGP?
# Run: python correlation_signs.py

import sympy as sp
line = "=" * 78
def hdr(s): print(line); print(s); print(line)

t, w, m1, m2, d, q = sp.symbols('t omega m1 m2 d q', positive=True)
M = m1 + m2
e = sp.Matrix([sp.cos(w*t), sp.sin(w*t), 0])     # unit orbital direction e(t)
x1 = (m2/M)*d*e
x2 = -(m1/M)*d*e                                  # opposite side of CM

def tavg(expr):
    T = 2*sp.pi/w
    return sp.simplify(sp.integrate(sp.expand_trig(expr), (t, 0, T))/T)

# ----------------------------------------------------------------------
hdr("MONOPOLE (charge q*m): conserved, no radiation")
Q1, Q2 = q*m1, q*m2
print(f"Q1 = q*m1, Q2 = q*m2 ; Qdot = 0 (charge conserved) -> no monopole radiation")

# ----------------------------------------------------------------------
hdr("DIPOLE  D_i = q * m * x_i  -- does the correlation SUBTRACT?")
D1 = q*m1*x1
D2 = q*m2*x2
print(f"body-1 dipole  q m1 x1 = {sp.simplify(D1.T)} ")
print(f"body-2 dipole  q m2 x2 = {sp.simplify(D2.T)} ")
print(f"are they anti-parallel?  D1 + D2 = {sp.simplify((D1+D2).T)}")
# cross-correlation of the 2nd derivative (dipole radiation ~ <Ddot^2>)
D1dd, D2dd = sp.diff(D1, t, 2), sp.diff(D2, t, 2)
self1 = tavg((D1dd.T*D1dd)[0]); self2 = tavg((D2dd.T*D2dd)[0])
cross = tavg((D1dd.T*D2dd)[0])
print(f"<Ddot1^2> = {self1}")
print(f"<Ddot2^2> = {self2}")
print(f"2*<Ddot1:Ddot2> (CROSS) = {sp.simplify(2*cross)}   <-- NEGATIVE")
print(f"total <(Ddot1+Ddot2)^2> = {sp.simplify(self1+2*cross+self2)}")
print("=> DIPOLE: correlation SUBTRACTS exactly (total = 0). The author is RIGHT.")
print("   TGP ALREADY banks this: scalar dipole D = q*Sum(m_i x_i) = 0 (T0, Pdot).")

# ----------------------------------------------------------------------
hdr("QUADRUPOLE  M_ij = m * x_i x_j  -- can the correlation also subtract?")
M1 = m1*(x1*x1.T)
M2 = m2*(x2*x2.T)
# are the two quadrupole tensors parallel (collinear) in tensor space?
M1d3, M2d3 = sp.diff(M1, t, 3), sp.diff(M2, t, 3)
def inner(A, B): return sp.simplify(sum(A[i, j]*B[i, j] for i in range(3) for j in range(3)))
i11 = tavg(inner(M1d3, M1d3))
i22 = tavg(inner(M2d3, M2d3))
i12 = tavg(inner(M1d3, M2d3))
align = sp.simplify(i12/sp.sqrt(i11*i22))         # normalized correlation (cos angle)
print(f"M1_ij = m1 x1 x1  ~  +e_i e_j  (sign of x1 SQUARED away)")
print(f"M2_ij = m2 x2 x2  ~  +e_i e_j  (sign of x2 SQUARED away)")
print(f"2*<Mddd1:Mddd2> (CROSS) = {sp.simplify(2*i12)}   <-- POSITIVE")
print(f"normalized correlation <Mddd1:Mddd2>/sqrt(<1><2>) = {align}")
print("=> QUADRUPOLE: the two tensors are PERFECTLY ALIGNED (correlation = +1).")
print("   They CANNOT subtract: M ~ x*x is EVEN in x, so opposite CM positions")
print("   give the SAME tensor direction. Constructive is geometrically forced.")

# ----------------------------------------------------------------------
hdr("WHY dipole subtracts but quadrupole adds (the parity argument)")
print("Position-odd moments (dipole ~ x, momentum ~ x-dot): opposite CM sides give")
print("OPPOSITE signs -> anti-parallel -> SUBTRACT -> cancel (= conservation laws:")
print("center of mass fixed, total momentum zero).")
print("Position-EVEN moments (quadrupole ~ x x): opposite sides give the SAME sign")
print("(minus times minus) -> parallel -> ADD. No cancellation possible.")
print("This is not a modelling choice; it is the parity of the multipole in x.")

# ----------------------------------------------------------------------
hdr("CONSEQUENCE for the falsification")
print("- The 'subtraction by correlation' the author invokes is REAL and is the")
print("  DIPOLE. TGP already exploits it fully: scalar dipole = 0 (T0).")
print("- It does NOT extend to the quadrupole (parity-even -> forced additive).")
print("- The surviving quadrupole is exactly the (1/6) P_GR scalar channel that")
print("  falsifies. Its cross term is +, not -, by geometry, so it cannot be")
print("  cancelled by any redecomposition into P1/P2/P3.")
print("- To make the quadrupole subtract you would need a SECOND quadrupole source")
print("  of OPPOSITE tensor phase carrying ~P_GR -> a spin-2 (helicity-2) field")
print("  anti-phased to the scalar -> again a graviton (forbidden axiom).")
print("=> The author's correlation-subtraction is correct physics, already used")
print("   where parity allows it (dipole=0). It cannot reach the parity-even")
print("   quadrupole. Verdict unchanged; R = 1/6 stands.")
