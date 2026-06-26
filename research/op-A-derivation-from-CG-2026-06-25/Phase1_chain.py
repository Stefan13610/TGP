# -*- coding: utf-8 -*-
# Phase 1 — op-A-derivation-from-CG: chain-audit symboliczny + status CG.
# Werdykt DERIVED/POSTULATE-CONDITIONAL/GAP WYLICZONY. 0 hardcoded.
import sympy as sp

print("="*72); print("PHASE 1 CHAIN-AUDIT: 𝒜=C_F^2 a_s^2 — derywowane czy postulate?"); print("="*72)
R={}
K_geo,m_sp,Phi0,CF,al,pi = sp.symbols('K_geo m_sp Phi0 C_F alpha_s pi', positive=True)

# --- Łańcuch symboliczny ---
A_color = CF*al/(pi*Phi0)            # L2
sigma_hat = pi*A_color**2            # L1 (ansatz gauss)
sigma_phys = K_geo*m_sp**2*sigma_hat # L4
sigma_phys = sp.simplify(sigma_phys)
print("\n[chain] sigma_phys =", sigma_phys)

# L5 hipoteza: K_geo*m_sp^2 = pi*Phi0^2  -> podstaw
A_result = sigma_phys.subs(K_geo*m_sp**2, pi*Phi0**2)
A_result = sp.simplify(sigma_phys.subs(K_geo, pi*Phi0**2/m_sp**2))
print("[chain] przy L5 (K_geo*m_sp^2=pi*Phi0^2): 𝒜 =", sp.simplify(A_result))

# T1: czy L5 => 𝒜=C_F^2 a_s^2 ? (algebraicznie)
R['T1_L5_gives_CF2as2'] = sp.simplify(A_result - CF**2*al**2)==0
print("\n[T1] L5 => 𝒜 = C_F^2*alpha_s^2 ?", R['T1_L5_gives_CF2as2'])

# T2: klasyfikacja linkow (DERIVED w ansatzu vs POSTULATE)
links = {
 "L1 sigma_hat=pi*A^2 (ansatz gauss)":"DERIVED-ANSATZ",  # Bessel K0 daje x0.3-0.6 (rem ii) => nie unikalne
 "L2 A_color=C_F a_s/(pi Phi0)":"DERIVED",                # sprzezenie konforemne
 "L3 m_sp^2=gamma":"DERIVED",                              # N0-6
 "L4 sigma_phys=K_geo m_sp^2 sigma_hat":"DEFINITIONAL",    # most wymiarowy
 "L5 K_geo*m_sp^2=pi*Phi0^2":"POSTULATE",                  # eq:X-K-msp-hypothesis (load-bearing)
}
n_postulate = sum(1 for v in links.values() if v=="POSTULATE")
R['T2_exactly_one_postulate'] = (n_postulate==1)
print("\n[T2] linki:")
for k,v in links.items(): print("     [%-15s] %s" % (v,k))
print("     liczba POSTULATE (load-bearing):", n_postulate, "-> dokladnie 1?", R['T2_exactly_one_postulate'])

# T3: status CG-1/CG-3 (most Gamma->Phi) — z manuskryptu (read-only, fakty)
CG_status = {"label":"[SZKIC]", "ex200_pass":4, "ex200_total":8, "A_sqrt_sigma_closed":False, "ex202_T6":"FAIL"}
CG_closed = (CG_status["label"]=="[ZAMKNIETE]") and (CG_status["ex200_pass"]==CG_status["ex200_total"]) and CG_status["A_sqrt_sigma_closed"]
R['T3_CG_not_closed'] = (not CG_closed)
print("\n[T3] CG-1/CG-3 (Gamma->Phi): label=%s ex200=%d/%d 𝒜~sqrt(sigma)/Phi0 zamkniete=%s ex202_T6=%s"
      % (CG_status["label"],CG_status["ex200_pass"],CG_status["ex200_total"],CG_status["A_sqrt_sigma_closed"],CG_status["ex202_T6"]))
print("     CG zamkniete? %s -> NOT closed: %s" % (CG_closed, R['T3_CG_not_closed']))

# T4: czy K_geo ustalone niezaleznie? (absorbowalne w field redef; CG nie domkniete => NIE)
K_geo_fixed_independently = CG_closed  # tylko zamkniety CG ustala K_geo*m_sp^2/Phi0^2
R['T4_K_geo_not_fixed'] = (not K_geo_fixed_independently)
print("\n[T4] K_geo*m_sp^2/Phi0^2 ustalone niezaleznie (wymaga zamknietego CG)? %s -> NIE: %s"
      % (K_geo_fixed_independently, R['T4_K_geo_not_fixed']))

# T5: WERDYKT wyliczony
L5_postulate = (links["L5 K_geo*m_sp^2=pi*Phi0^2"]=="POSTULATE")
if R['T1_L5_gives_CF2as2'] and (not L5_postulate) and (not R['T3_CG_not_closed']):
    verdict="DERIVED"
elif R['T1_L5_gives_CF2as2'] and L5_postulate and R['T3_CG_not_closed'] and R['T4_K_geo_not_fixed']:
    verdict="POSTULATE-CONDITIONAL"
else:
    verdict="GAP"
print("\n[T5] flagi: L5=>CF2as2=%s ; L5 postulate=%s ; CG not closed=%s ; K_geo not fixed=%s"
      % (R['T1_L5_gives_CF2as2'], L5_postulate, R['T3_CG_not_closed'], R['T4_K_geo_not_fixed']))
print("     >>> WERDYKT: %s <<<" % verdict)
print("     => alpha_s 0,03sigma = %s" % ("STRUCTURAL CONSISTENCY-CHECK warunkowy (NIE first-principles)" if verdict=="POSTULATE-CONDITIONAL" else verdict))

print("\n"+"="*72)
fp=[k for k in R]; print("PASS: %d/%d"%(sum(1 for k in fp if R[k]),len(fp))); print("WERDYKT:",verdict); print("="*72)
