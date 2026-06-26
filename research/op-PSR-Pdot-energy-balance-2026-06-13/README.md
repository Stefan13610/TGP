---
title: "op-PSR-Pdot-energy-balance-2026-06-13"
folder_status: closed-resolved
claim_status: "TRIGGERED-FALSIFIED (pending user ratification)"
date_opened: 2026-06-13
date_closed: 2026-06-13
pr_candidate: PR-025
---

# op-PSR-Pdot-energy-balance — bilans energii sektora radiacyjnego vs Ṗ_b pulsarów

**Pytanie:** czy LIVE TGP (skalar Φ + kanał σ_ab) reprodukuje obserwowany zanik
orbity PSR J0737−3039 na poziomie bilansu energii?

**Wynik (21/21 sympy PASS, single-session sprint):**
- P_φ/P_GR = **1/6 EXACT** (ab initio; wzór autora TGP potwierdzony),
- kanał σ przy LOCKED m_σ = 0.71 meV **nie propaguje** przy ω pulsarowych (15 rzędów poniżej progu),
- amplitude-lock h_TT^σ = h_TT^GR **nie pinuje strumienia energii** (κ_E niezależne — znalezisko R1 #23),
- Gałąź A: R = 1/6 → **13 227σ**; Gałąź B: R = 7/6 → **2 646σ** vs R_obs = 0.999963 ± 0.000063 (Kramer 2021).

**Werdykt:** TRIGGERED-FALSIFIED — sektor radiacyjny LIVE TGP wykluczony przez
istniejące dane pulsarowe przy każdym nie-strojonym przypisaniu energii.

Pliki: [[Phase0_balance.md]] (LOCK) → [[Phase1_sympy.py]]/[[Phase1_sympy.txt]] →
[[Phase1_results.md]] → [[Phase_FINAL_close.md]].

Geneza: external expert review 2026-06-13 (dyskusja trylemat energetyczny);
autoryzacja usera "ok, sprawdźmy to ;)".
