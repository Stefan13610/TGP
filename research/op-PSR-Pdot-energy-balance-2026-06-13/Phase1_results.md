---
title: "Phase 1 results — energy balance derivation + observational verdict"
date: 2026-06-13
type: phase-results
phase: 1
status: 🔴 EXECUTED — 21/21 sympy PASS — both LIVE branches TRIGGERED at >5σ
parent: "[[./README.md]]"
predecessor: "[[./Phase0_balance.md]] (LOCKED before computation)"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
tags:
  - pulsar-Pdot
  - energy-balance
  - scalar-quadrupole
  - sigma-channel
  - J0737-3039
---

# Phase 1 — wyniki

## §0 — Executive summary

**21/21 sympy PASS, 0 hardcoded.** Wszystkie cztery falsyfikatory rozstrzygnięte:

| Falsyfikator | Werdykt |
|---|---|
| F-PDOT-A | P_φ = (16/15)Gμ²d⁴ω⁶/c⁵ ab initio ⟹ **P_φ/P_GR = 1/6 EXACT** (claim usera POTWIERDZONY) |
| F-PDOT-B | ℏω_GW = 9.36×10⁻¹⁹ eV vs m_σ = 7.1×10⁻⁴ eV (ratio 1.3×10⁻¹⁵) ⟹ **NON_PROPAGATING** — Gałąź A obowiązuje |
| F-PDOT-C | **Gałąź A: R = 1/6 ⟹ 13 227σ TRIGGERED. Gałąź B (κ_E=1): R = 7/6 ⟹ 2 646σ TRIGGERED.** Cross-check B1913+16: 520σ / 105σ — te same werdykty |
| F-PDOT-D | Bound na shift stosunku z e = 0.0878: ~5%; luki 83%/17% ⟹ werdykt ROBUST |

## §1 — Łańcuch wyprowadzenia (wszystko sympy-exact)

1. **T0:** dipol skalarny D = qΣm_ix_i ≡ 0 dla orbity w układzie CM (potwierdzenie usera).
2. **T2a:** kalibracja statyczna: U_int = −Gm₁m₂/d ⟺ **q² = 4πG** — wyprowadzone, nie założone. Skalar niesie CAŁĄ statykę Newtona (σ przy 0.71 meV ma zasięg 0.278 mm — T4.2).
3. **T1+T3:** pole dalekie ψ̇ = (q/8πr)n_in_jM⃛_ij; strumień ∮⟨ψ̇²⟩r²dΩ z jawnymi całkami kątowymi i uśrednieniem po okresie:
   $$P_φ = \frac{16}{15}\frac{G\mu^2 d^4\omega^6}{c^5} = \frac{G}{30c^5}\langle \dddot I_{ij}\dddot I_{ij}\rangle = \frac{1}{6}P_{GR}\ \text{EXACT}$$
   Wzór usera z dyskusji 2026-06-13 zreprodukowany w 100%.
4. **T4:** bramka propagacji: pole masywne nie unosi energii przy ω < m (k² < 0, mod ewanescentny). Kanał σ przy LOCKED m_σ = 0.71 meV jest **martwy radiacyjnie** na częstościach pulsarowych (15 rzędów wielkości poniżej progu).
5. **T5 — KLUCZOWE ZNALEZISKO STRUKTURALNE (R1-kandydat):**
   warunek amplitudowy h_TT^σ = h_TT^GR ⟺ λ·ξ_eff = 16πG (JEDEN warunek);
   strumień energii κ_E = ξ_eff²/(16πG) — **NIEZALEŻNA druga kombinacja**.
   Amplitude-lock T3.4 **nie pinuje energii**. κ_E = 1 wymaga dodatkowo ξ_eff = λ,
   czego żaden LOCK w repo nie asertuje. Zamknięty wynik "h_TT^σ/h_TT^GR = 1"
   był niemym założeniem, że amplituda ⟹ energia — non sequitur.
6. **T6:** werdykty obserwacyjne (J0737−3039, Kramer et al. 2021, R_obs = 0.999963 ± 0.000063):
   - Gałąź A (LIVE, m_σ = 0.71 meV): R = 1/6 → **13 227σ**
   - Gałąź B (mechanizm-v, κ_E = 1 kanoniczne): R = 7/6 → **2 646σ**
   - Ratunek Gałęzi B wymagałby κ_E = 0.833296 ± 0.000315 (precyzja 3.8×10⁻⁴) — NOWA kotwica strojona daną ⟹ forbidden move §6.2 Phase 0.
   - Gałąź C: brak mechanizmu kasowania kanału skalarnego przy q² = 4πG (statyka).
   - Gałąź D (screening): samo-destrukcyjna — ekranowanie radiacji skalarnej ekranuje statykę ⟹ brak wiązania orbity.
7. **T7:** poprawki ekscentryczności (e = 0.0878): bound ~5% na stosunek kanałów; nie mogą zmostkować luk 0.833 / 0.167. Werdykt ROBUST.

## §2 — Zgodność z pre-rejestracją

Phase 0 §5 pre-deklarowane oczekiwanie zrealizowane DOKŁADNIE we wszystkich
czterech punktach (1/6 exact; NON_PROPAGATING; 1.3×10⁴σ; κ_E nieokreślone +
7/6 przy wyborze kanonicznym). Zero modyfikacji progów. Zero nowych stałych.

## §3 — Honest notes

- P_GR = (G/5c⁵)⟨I⃛²⟩ jest kotwicą zewnętrzną (textbook), nie wyprowadzeniem TGP — zgodnie z LOCK §2 Phase 0.
- Rachunek na orbicie kołowej; e-korekta oszacowana boundem, nie policzona kanał-po-kanale (wystarczające przy lukach >10³σ; F-PDOT-D PASS).
- Sensitivities ciał zwartych (à la Damour–Esposito-Farèse) NIE uwzględnione — mogłyby jedynie WŁĄCZYĆ dipol (−1PN, silniej wzmocniony), tj. pogorszyć werdykt. Pominięcie jest konserwatywne na korzyść TGP.
- Wynik dotyczy LIVE TGP (Gałąź A) oraz jedynego nie-strojonego wariantu mechanizmu-v (Gałąź B kanoniczna). Nie dotyczy hipotetycznych rozszerzeń poza LOCKED inventory.
