---
title: "Phase0_balance — LOCK (mini-cykl analityczno-decyzyjny, wzorzec L04): klasa dynamiki substratu a znak W — spektrum symplektyczne JL̂ na istniejących tłach łańcucha"
date: 2026-08-31
type: phase0-lock
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-31: 'ok, załuż nowy mini cykl' (po weryfikacji konwencji W_source/V_energy: relabeling jest widmowo inwariantny — statyka maszynerii 2 pinuje ω²=k²−γ w KAŻDEJ dynamice II rzędu i w gradient flow; jedyna nieprzebadana klasa = symplektyczna I rzędu)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase2_backgrounds.npz]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/README.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Phase 0 — LOCK mini-cyklu `op-symplectic-Jspectrum`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria,
kontrole, progi i forbidden moves zamknięte PRZED kodem.**

---

## 0. Pytanie binarne (Q) i status ontologiczny

**Czy tła sektora tachionowego (łańcuchy z cyklu bloch: d ∈ {3π, 4π, 6π},
Hessian z ujemnym kierunkiem λ_min ≈ −1.22), niestabilne w dynamice
II rzędu (Q-FAIL), są SPEKTRALNIE STABILNE w dynamice symplektycznej
I rzędu (NLS/GP-podobnej): σ(JL̂) ⊂ iℝ?**

- TAK ⟹ „dwa sektory znaku W" mogą być „dwiema dynamikami jednej akcji":
  solitonowy znak siedzi w statyce, a fluktuacje są stabilne orbitalnie —
  bez ducha i bez ręcznej zmiany znaku. Pierwszy pozytywny nośnik
  stabilności maszynerii 2 (materiał do decyzji aksjomatycznej autora:
  klasa dynamiki substratu).
- NIE ⟹ także klasa symplektyczna nie ratuje sektora tachionowego w 1D —
  hipoteza stabilizacji przez zmianę dynamiki zamknięta negatywnie
  w klasie zbadanej; pozostaje 3D i/lub decyzja aksjomatyczna bez tej
  podpory.

**Status: MODEL-EXTENSION / analytical-decision (wzorzec L04).** Wynik
NIE zmienia core; zasila decyzję ontologiczną (NEEDS user-gate).
Delimitacja od #63 V2 (jawna): tam Q-balle ω>0 w modelu f_ε z kryterium
VK (slope-positive) — tu STATYCZNE (ω=0) tła łańcucha w akcji
kanonicznej i pełne spektrum JL̂; obiekt i pytanie inne, brak dublowania.

## 1. Ustalenia dziedziczone (dowody na wejściu, NIE do powtarzania)

- Relabeling W_source/V_energy jest widmowo inwariantny: każde zanurzenie
  II rzędu odtwarzające statykę ∇²g = g²(1−g) daje ω² = k² − γ (rachunek
  konwencyjny 2026-08-31, sesja główna; zgodny z P1a cyklu bloch do
  8.3e−5). Gradient flow: stabilność rządzi ten sam Hessian ⟹ też
  niestabilne. Formalizacja tych dwóch faktów = Phase 1 (sympy).
- Tła: `../op-bloch-chain-stability-2026-08-31/Phase2_backgrounds.npz`
  (3π; 4π; 6π 2-garb i 1-garb; residua ≤ 5.7e−12) — UŻYĆ, nie liczyć
  od nowa. Hessian L₊ = operator z Phase 3 tamtego cyklu (kotwica:
  λ_min = −1.222191 dla 3π — gate reprodukcji ≤ 1e−4).

## 2. Model ZAMKNIĘTY

- **Zanurzenie symplektyczne (M-J):** pole zespolone u(x,t),
  dynamika i∂_t u = δE/δu*, energia kanoniczna rozszerzona
  E[u] = ∫ [½ K(|u|) |∂_x u|² + U(|u|)] dx, K(g)=g⁴,
  U(g)=g⁷/7−g⁸/8 (akcja kanoniczna — ta sama, której statyka = M-C
  cyklu bloch; tło rzeczywiste u₀ = g_d(x), ω=0).
- **Linearyzacja:** u = (g_d + a + ib)e^{i·0}; blokowo
  ∂_t(a,b) = J·diag-sprzężenie (L₊, L₋): L₊ = Hessian amplitudowy
  (MUSI odtworzyć operator Phase 3 cyklu bloch — gate), L₋ = operator
  fazowy wyprowadzony z E[u] (derywacja sympy w Phase 1; własność
  kontrolna: L₋ g_d = 0 exact — mod fazowy/cechowania, osiągalny FAIL).
- **Kryterium spektralne:** wartości własne λ problemu JL̂; stabilność
  spektralna ⟺ max Re λ ≤ tol, tol := max(1e−6, 1e−3·max|Im λ|)
  na siatce bazowej, zbieżnie (patrz §3).
- **Warunki brzegowe:** periodyczne na komórce d (jak w cyklu bloch);
  Bloch k=0 PRIMARY (tam siedział argmin Hessianu), pełny skan k
  (≥8 punktów w [0, π/d]) jako rozszerzenie zalockowane.
- **Rejestr WEJŚĆ:** wybór złożenia zespolonego u z g (moduł pola) =
  decyzja modelowa [INPUT-ONTO], flagować w każdym wyniku; β=γ=1.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — formalizacja analityczna (sympy, zero numeryki siatkowej)
- P1a: dowód inwariantności widmowej zanurzeń II rzędu (energia vs
  źródło) — obie dyspersje ω²=k²−γ symbolicznie.
- P1b: gradient flow — stabilność ⟺ Hessian ≥ 0 (symbolicznie).
- P1c: derywacja L₊ i L₋ z E[u] (sympy); tożsamość L₋g_d = 0 on-shell
  (statyczne EL) — exact. FAIL tożsamości ⟹ STOP (derywacja błędna).

### Phase 2 — bramka maszynerii numerycznej (osiągalne FAIL-e w OBIE strony)
- **C1 (znany STABILNY):** soliton NLS kubicznego 1D
  (iu_t + u_xx + |u|²u = 0): maszyneria MUSI dać σ(JL̂) ⊂ iℝ
  (max Re λ ≤ tol) mimo ujemnego kierunku L₊ — dokładnie mechanizm,
  o który pytamy.
- **C2 (znany NIESTABILNY):** soliton NLS z nieliniowością |u|⁶u
  (nadkrytyczna, p>5): maszyneria MUSI znaleźć λ z Re λ > 0.
- **C3 (próżnia, analitycznie):** dyspersja JL̂ wokół u=1 sektora
  tachionowego policzona zamknięcie w Phase 1, numerycznie odtworzona
  ≤ 1e−3.
- Gate reprodukcji L₊: λ_min(3π) = −1.222191 ± 1e−4 (kotwica bloch).
- Gate Phase 2 FAIL ⟹ STOP (kod nieważny).

### Phase 3 — RACHUNEK CENTRALNY
- σ(JL̂) dla 4 teł (3π; 4π; 6π×2), k=0 PRIMARY + skan k; ≥2 siatki
  (N i 2N); zbieżność: |Δ max Re λ| ≤ 0.01·max(max Re λ, 0.01).
- **Q-PASS:** max Re λ ≤ tol zbieżnie dla WSZYSTKICH 4 teł i wszystkich k,
  przy czystych C1/C2/C3.
- **Q-FAIL:** istnieje tło z Re λ > próg zbieżnie (raportować per tło;
  jeśli część teł stabilna — wynik MIXED, raportować per tło bez
  uśredniania).
- **Q-INCONCLUSIVE:** niezbieżność — jako niezbieżność.

### Phase 4 (warunkowa, TYLKO przy Q-PASS) — sanity nieliniowe
- Ewolucja split-step i∂_t u = δE/δu* z perturbacją 0.01 wzdłuż dawnego
  modu runaway, do t = 3·t*_ref (3.62·3); kryterium: brak ucieczki
  (‖|u|−g_d‖∞ ≤ 3× amplituda początkowa). Gate zachowania Q=∫|u|² i E
  ≤ 1e−8. (Rozbieżność Phase 3↔4 = incydent do raportu, nie do ukrycia.)

## 4. Forbidden moves
1. Zero zmian kryteriów/progów/tol/list teł po starcie; korekty tylko
   dla udokumentowanego błędu implementacji (correction_note PRZED
   użyciem wyniku; pierwotne outputy zachować).
2. Kontrole C1/C2/C3 i tożsamość L₋g_d=0 nieusuwalne; każdy test
   z osiągalnym FAIL.
3. Tła WYŁĄCZNIE z Phase2_backgrounds.npz (zakaz relaksowania nowych —
   to nie ten cykl); katalogów innych cykli nie dotykać.
4. Wyniki negatywne/MIXED wprost; zakaz interpretacji ontologicznej
   w close (wynik = wejście do decyzji autora, nie decyzja).
5. Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowany; git nieużywany
   (commit robi sesja główna).
6. Higiena ścieżek: bez `cd`; weryfikacja materializacji po każdym zapisie.

## 5. Deliverables
`Phase1_analytic.py`+output, `Phase2_gate_nls.py`+output,
`Phase3_Jspectrum.py`+output(+json), warunkowo `Phase4_evolution.py`+output,
`Phase_FINAL_close.md`, `NEEDS.md` (user-gated), `README.md` (log faz).

## 6. Drzewo decyzyjne
```text
P1c FAIL → STOP (derywacja L₋ błędna; raport)
P2 gate FAIL → STOP (kod nieważny)
Q-PASS → NEEDS: materiał do decyzji aksjomatycznej „klasa dynamiki
   substratu" (dwie dynamiki jednej akcji); kandydat: VK/continuation
   ω≠0 dla łańcuchów + wersja 3D
Q-FAIL → NEEDS: klasa symplektyczna nie ratuje 1D; decyzja o znaku W
   wraca do autora bez tej podpory; pozostaje 3D
MIXED → raport per tło; NEEDS: charakteryzacja które tła i czemu
```

---

**LOCK ZAMKNIĘTY 2026-08-31. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
