---
title: "Phase2_correction_note — udokumentowana korekta błędu implementacji (detektor breakdown), opisana PRZED użyciem wyników (LOCK §4.1)"
date: 2026-08-23
type: method-note
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: CORRECTION-BEFORE-USE
related:
  - "[[Phase_method_decisions.md]]"
  - "[[Phase0_balance.md]]"
---

# Korekta 1: detektor breakdown w `evolve()`

**Objaw (pierwszy przebieg Phase2, zachowany w Phase2_output.txt):**
zalockowany baseline P2a (#63 V3 = soliton μ) reprodukowany czysto
(λ_min(w1)=−1.3896 vs −1.389; t*=3.62 przy dt=0.004 i 0.002;
gate ≤1.64e−08 PASS). DODANY przeze mnie komparator e (g₀=1.24915)
dał gate=3.56e+95 FAIL przy a=+0.01, mimo że dynamika do t≈10 była
regularna.

**Diagnoza (błąd implementacji, nie fizyka):** detektor wyjścia
z dziedziny (`min(g) ≤ 0`) był sprawdzany tylko w punktach próbkowania
(co 0.02 t = co 5 kroków RK4). Wolny runaway e przechodzi przez g≈0
między próbkami; przy g→0⁺ człony log/1/g² eksplodują i licznik energii
zatruwa się overflowem PRZED detekcją. Dla μ (szybki, stromy runaway)
artefakt nie występował — dlatego #63 go nie widział.

**Korekta (przed użyciem wyników):** detekcja wyjścia z dziedziny
w KAŻDYM kroku RK4, z progiem podłogi modelu `g_floor = 0.01`
(profil statyczny #62/#63 i tak clampuje g≥1e−3; g=0.01 leży głęboko
pod każdym tłem — to nadal „g→0" z #63, wykryte zanim integrator
opuści reżim stosowalności). Energia liczona wyłącznie w stanach
przed-breakdown. Kryteria LOCKa (gate 1e−6, t*, ω², punkty d,
konfiguracje) — NIETKNIĘTE.

**Zakres gate'u P2a:** zalockowany gate P2a dotyczy reprodukcji #63 V3
(μ) — ta PRZESZŁA w obu przebiegach. Komparator e pozostaje dodatkiem
(dozwolonym), raportowanym z własnym gate po korekcie.

---

# Korekta 2 (Phase 3): próg Newtona poniżej podłogi arytmetyki

**Objaw (pierwszy przebieg Phase3, zachowany):** tła d=4 (N=800), d=2
i kontrola q=0.3 raportowane „Newton NIE ZBIEGL" przy ‖R‖∞ = 1.4–6.3e−11
(próg: 1e−11), podczas gdy d=8/6 i d=4 (N=400) zbiegły normalnie.

**Diagnoza (błąd implementacji, nie fizyka):** residuum zawiera
(ψ₊−2ψ+ψ₋)/h²; przy h=0.005 i ψ≈1.28 szum kasowania w double wynosi
~4·ε·ψ/h² ≈ 5e−11 — mój próg 1e−11 leży PONIŻEJ osiągalnej podłogi.
Iteracja osiągnęła plateau na podłodze maszynowej, czyli tło JEST
samouzgodnione do granic arytmetyki.

**Korekta (przed użyciem wyników):** próg zbieżności skalowany podłogą:
tol = max(1e−10, 20·ε·max(ψ)/h²). Kryteria Q2 (znak ω²_min, punkty d,
kontrola d=∞, klasyfikacja odpowiedzi) — NIETKNIĘTE.
