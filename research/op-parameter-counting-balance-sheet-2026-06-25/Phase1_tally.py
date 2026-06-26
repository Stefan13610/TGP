# -*- coding: utf-8 -*-
# Phase 1 — op-parameter-counting-balance-sheet-2026-06-25
# Uczciwy tally: N_free (ciagle wolne param) + N_axiom (dyskretne) vs N_pred (genuine).
# Werdykt HEADLINE-* WYLICZONY z progow 4/12. 0 hardcoded werdyktu.
# Konwencja symetryczna z SM: symetrie/grupy NIE liczone jako parametr po obu stronach.

print("="*72)
print("PHASE 1 BALANCE-SHEET: uczciwy licznik inputow TGP vs headline '3 inputy'")
print("="*72)

# ---- LEDGER (klasyfikacja per kryteria §0.2; cytat zrodla w .md) ----
# (nazwa, klasa, zrodlo)
ledger = [
    # FREE-NUMERICAL (ciagle wolne parametry)
    ("g0_e = 0.86941",            "FREE",   "README; dodatekJ (baza leptonowa, 1->3)"),
    ("Phi_0",                     "FREE",   "FOUNDATIONS §3.5.3: 'EFT scale-dependent FREE parameter'"),
    ("c0 / C_sigma",              "FREE",   "#37: wolny parametr UV (GW amplituda)"),
    ("quark anchor m_d",          "FREE",   "#41: down sektor m_1"),
    ("quark anchor m_s",          "FREE",   "#41: down sektor m_2"),
    ("quark anchor m_u",          "FREE",   "#41: up sektor m_1"),
    ("quark anchor m_c",          "FREE",   "#41: up sektor m_2"),
    ("N_e (e-folds inflacji)",    "FREE",   "sek00:338 n_s=1-2/N_e, r; N_e~60 input"),
    # TRADED (pozorny win = wymiana input<->postulate; NIE redukcja)
    ("Omega_Lambda / g_tilde",    "TRADED", "FOUNDATIONS §3.5.3.1: g~0.98 postulate; sek00 g~N_f=5 (sporne)"),
    ("alpha_s(M_Z)",              "TRADED", "#41: =sqrt(A)/C_F warunkowe (CG nie domkniete); input SM lub predykcja"),
    # STRUCTURAL-AXIOM (dyskretne; osobno, NIE w N_free; konwencja symetryczna z SM)
    ("N = 3 (generacje)",         "AXIOM",  "README; liczba calkowita"),
    ("Z_2 symetria",             "AXIOM",  "FOUNDATIONS §1 (jak grupa cechowania SM - nie param)"),
    ("alpha = 2 (klasa C1-C3)",   "AXIOM",  "#36/#38/#39: nieredukowalnie aksjomatyczne"),
    ("beta = gamma (vacuum)",     "AXIOM",  "FOUNDATIONS §3.5.3.1: fine-tuning LUB RG-FP (otwarte)"),
    ("phi-FP (golden ladder)",    "AXIOM",  "dodatekJ: hipoteza topologiczna"),
    ("Koide B=sqrt2 / K=2/3",     "AXIOM",  "dodatekX: dla leptonow algebraiczne z ODE"),
    # DERIVED (genuine predykcje - wins; N_out > N_in w bloku)
    ("leptony m_e,m_mu,m_tau",    "DERIVED","g0_e ->3 (phi-FP r21 + Koide r31); 0.006%  [BLOK 1->3]"),
    ("Koide K=2/3",               "DERIVED","tozsamosc algebraiczna"),
    ("kwarki m_b,m_t",            "DERIVED","#41: r_31 z r_21+A (4->2; m_t scheme-zal.)"),
    ("n_s = 1-2/N_e",             "DERIVED","Starobinsky (sek00:338)"),
    ("r (tensor-scalar)",         "DERIVED","Mukhanov-Sasaki"),
    ("c_GW = c",                  "DERIVED","strukturalne (exact)"),
    ("proton lifetime = inf",     "DERIVED","strukturalne"),
    ("m_W = 80.354",              "DERIVED","sek00:337 (0.01 sigma)"),
    ("alpha_s z A=C_F^2 a_s^2",   "DERIVED","#41: 0.03 sigma (warunkowe CG)"),
    # WEAK / FIT (NIE liczone jako genuine DERIVED)
    ("PMNS angles (zeta.1)",      "FIT",    "sek00:323 zeroth-order drifty 8.6/12.6/15.6%; iota.1+mu.1 WITHDRAWN"),
    ("m_H = v*57/112",            "FIT",    "sek00:341 uzywa v + algebraiczny 57/112"),
]

from collections import Counter
cnt = Counter(c for _,c,_ in ledger)
N_free  = cnt["FREE"]
N_traded= cnt["TRADED"]
N_axiom = cnt["AXIOM"]
N_deriv = cnt["DERIVED"]
N_fit   = cnt["FIT"]

# N_free uczciwy: FREE + TRADED (TRADED nie redukuje, wiec liczy sie jako utrzymany input)
N_free_honest = N_free + N_traded

print("\n--- LEDGER ---")
for n,c,s in ledger:
    print("  [%-7s] %-26s | %s" % (c, n, s))

print("\n--- PODSUMOWANIE ---")
print("  N_free (ciagle wolne)        = %d" % N_free)
print("  N_traded (input<->postulate) = %d" % N_traded)
print("  N_free_honest (FREE+TRADED)  = %d   <-- vs headline '3'" % N_free_honest)
print("  N_axiom (dyskretne selekcje) = %d" % N_axiom)
print("  N_pred genuine DERIVED       = %d" % N_deriv)
print("  N_fit (slabe/wycofane)       = %d" % N_fit)

# SM baseline (symetryczna konwencja: grupy/symetrie nie liczone)
N_free_SM = 19  # 9 mas + 4 CKM + 3 cechowania + 2 Higgs + theta_QCD (bez neutrin)
print("\n  N_free^SM (symetrycznie)     = %d (+7 z masywnymi neutrinami)" % N_free_SM)

# ---- WERDYKT (wyliczony z progow 4/12) ----
lepton_block_derived = any("leptony" in n and c=="DERIVED" for n,c,_ in ledger)  # 1->3 genuine
mostly_fits = (N_fit > N_deriv)

if N_free_honest <= 4:
    verdict = "HEADLINE-HONEST"
elif (5 <= N_free_honest <= 12) and (N_free_honest < N_free_SM) and lepton_block_derived and (not mostly_fits):
    verdict = "HEADLINE-OPTIMISTIC"
else:
    verdict = "HEADLINE-MISLEADING"

print("\n--- WERDYKT (wyliczony) ---")
print("  N_free_honest=%d (prog HONEST<=4, OPTIMISTIC 5-12, MISLEADING>12)" % N_free_honest)
print("  N_free_honest << SM(%d)? %s ; blok leptonowy genuine 1->3? %s ; wiekszosc fity? %s"
      % (N_free_SM, N_free_honest < N_free_SM, lepton_block_derived, mostly_fits))
print("  >>> %s <<<" % verdict)

print("\n  Uczciwy headline (rekomendacja):")
print("  'Ok. %d genuine predykcji z ~%d wolnych parametrow + %d aksjomatow selekcji"
      % (N_deriv, N_free_honest, N_axiom))
print("   (vs SM ~%d); blok leptonowy najmocniejszy: 3 masy z 1 sprzezenia.'" % N_free_SM)
print("="*72)
