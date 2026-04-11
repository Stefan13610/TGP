#!/usr/bin/env python3
"""
koide_window_analytic_v47b.py -- Analytical estimate of alpha window for Koide.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2

print("KOIDE AVAILABILITY WINDOW (ANALYTICAL ESTIMATE)")
print("=" * 65)
print()
print("g0_crit = (2a+4)/(2a+1),  g0_mu approx phi * g0_e approx 1.41")
print("Tau needs: g0_mu < g0_tau < g0_crit")
print()

header = "%-6s %-10s %-10s %-10s %-10s"
print(header % ("alpha", "g0_crit", "g0_mu_est", "room", "tau?"))

g0_e_est = 0.868
g0_mu_est = PHI * g0_e_est

for alpha in np.arange(0.5, 4.1, 0.25):
    gc = (2 * alpha + 4) / (2 * alpha + 1)
    room = gc - g0_mu_est
    tau_fits = "YES" if room > 0.03 else "MARGINAL" if room > 0 else "NO"
    print("%-6.2f %-10.4f %-10.4f %-10.4f %-10s" %
          (alpha, gc, g0_mu_est, room, tau_fits))

print()
print("CRITICAL ALPHA VALUES:")

# g0_crit = g0_mu  =>  (2a+4)/(2a+1) = phi*g0_e
# (2a+4) = phi*g0_e*(2a+1)
# 2a(1 - phi*g0_e) = phi*g0_e - 4
# a = (phi*g0_e - 4) / (2*(1 - phi*g0_e))
a_muon_max = (g0_mu_est - 4) / (2 * (1 - g0_mu_est))
print("  alpha_max (muon hits collapse) = %.3f" % a_muon_max)

# For tau Koide: need room > ~0.03
# (2a+4)/(2a+1) = g0_mu_est + 0.03
target = g0_mu_est + 0.03
a_tau_max = (target - 4) / (2 * (1 - target))
print("  alpha_tau (tau Koide limit)    = %.3f" % a_tau_max)
print()

# TGP canonical
print("CANONICAL TGP (alpha=2):")
print("  g0_crit = 8/5 = 1.6")
print("  g0_mu   = %.4f" % g0_mu_est)
print("  room    = %.4f" % (1.6 - g0_mu_est))
print("  alpha/alpha_muon_max = %.2f" % (2.0 / a_muon_max))
print("  alpha/alpha_tau_max  = %.2f" % (2.0 / a_tau_max))
print()

# Fibonacci
print("FIBONACCI CONNECTION:")
print("  g0_crit = 8/5 = F(6)/F(5) (consecutive Fibonacci)")
print("  phi = lim F(n+1)/F(n) = %.6f" % PHI)
print("  g0_crit = %.6f" % 1.6)
print("  |g0_crit - phi| = %.6f (diff = %.2f%%)" %
      (abs(1.6 - PHI), abs(1.6 - PHI) / PHI * 100))
print()

# Check: 2*alpha+4 = 8 = F(6), 2*alpha+1 = 5 = F(5)
# This is because alpha=2, and 2*2+4=8, 2*2+1=5
# 8 and 5 are Fibonacci: F(6)=8, F(5)=5
# The ratio F(6)/F(5) = 8/5 converges to phi from below
print("  The Fibonacci connection is NUMERICAL COINCIDENCE:")
print("  2*alpha+4 = 8 = F(6) and 2*alpha+1 = 5 = F(5)")
print("  because alpha=2 happens to give these values.")
print("  For alpha=3: 10/7, which are NOT Fibonacci.")
print()

# Summary
print("=" * 65)
print("SUMMARY")
print("=" * 65)
print()
print("  Koide Q_K=3/2 exists for alpha in [0, approx %.1f]" % a_tau_max)
print("  Canonical TGP alpha=2 is at %.0f%% of this window." %
      (200 / a_tau_max))
print("  alpha=2 is NOT uniquely selected by Koide alone.")
print("  It IS selected by the TGP Lagrangian (K=g^4 => alpha=2).")
print("  The coincidence 8/5 = F(6)/F(5) is numerical, not structural.")
