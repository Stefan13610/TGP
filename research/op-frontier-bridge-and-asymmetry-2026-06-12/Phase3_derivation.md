---
title: "Phase 3 — moduł A: F-BA-1 (GAP-1) DERIVED — most EQ-5 na poziomie pola (j₀ = (3/8)m_Φ(t_*/t)² w jednostkach pola) · F-BA-4 (GAP-4) DERIVED_IN_CLASS — jedyność (U, ρ₀, k) = (2/3, 1/6π, 0); jednorodność WYMUSZONA · F-BA-5 (GAP-5) DERIVED — sferyczny atraktor a_l ∝ t^(−1/2) dla wszystkich l ≥ 2"
type: phase_result
status: PHASE3_COMPLETE
phase: 3
cycle: op-frontier-bridge-and-asymmetry-2026-06-12
created_date: 2026-06-12
authorization: "User 2026-06-12: 'Phase 3'"
sympy_script: "[[./Phase3_sympy.py]]"
sympy_output: "[[./Phase3_sympy.txt]]"
sympy: "13/13 PASS; 0 hardcoded T_pass; circularity guard FP13"
falsifier_resolved: "F-BA-1 = DERIVED (tożsamość wymiany energii z Lagrangianu EXACT; amplituda źródła w jednostkach pola; domknięcie ledgera) · F-BA-4 = DERIVED_IN_CLASS (jedyność w rozszerzonej klasie self-similar; jednorodność wymuszona Eulerem; naddeterminacja domyka się EXACT) · F-BA-5 = DERIVED (toy nierelatywistyczny zadeklarowany: a_l ∝ t^(−1/2) → 0 ∀ l ≥ 2)"
anti_lakatos_lock: PRESERVED
---

# Phase 3 — moduł A: most polowy EQ-5 / jedyność / sferyczność

## §0 — Verdict at a glance

| Falsyfikator | Werdykt | Esencja |
|---|---|---|
| **F-BA-1 (GAP-1)** | **DERIVED** | Tożsamość wymiany energii w jednostkach pola EXACT (z Lagrangianu); transfer zlokalizowany na ścianie T_area = c·j₀·σ; **amplituda źródła EQ-5 wyprowadzona: j₀(t) = (3/8)·m_Φ·(t_*/t)²**; ledger domknięty T = Ṁc² EXACT; wymiary EXACT; η wynika z regulatora Phase 2 (nie wstawione) |
| **F-BA-4 (GAP-4)** | **DERIVED_IN_CLASS** | W rozszerzonej klasie self-similar {u = Ux/t, ρ = ρ₀ξ^k/(Gt²)}: **jedyne rozwiązanie (U, ρ₀, k) = (2/3, 1/(6π), 0)**; jednorodność (k = 0) WYMUSZONA przez skalowanie Eulera; marginalność = trzeci, naddeterminujący warunek — domyka się EXACT |
| **F-BA-5 (GAP-5)** | **DERIVED** (toy zadeklarowany) | Wszystkie mody kształtu l ≥ 2: **a_l ∝ t^(−1/2)·osc → 0 — jednolite tempo zaniku** (sferyczny atraktor); l = 1 = dryf translacyjny (nie kształt); napęd ΔV shape-neutral |

## §1 — F-BA-1: most EQ-5 w jednostkach pola (GAP-1)

### §1.1 Tożsamość wymiany energii (z Lagrangianu, EXACT — FP1)

Dla pola z sourced EOM (φ̈ = φ″ − V′(φ) + J — struktura EQ-5/concept §3.3):

```
∂_t[½φ̇² + ½φ′² + V(φ)] − ∂_x[φ̇φ′] = φ̇·J     (residual = 0 symbolicznie)
```

**Gęstość transferu S_Φ → S_matter = φ̇·J w jednostkach pola** — to jest brakująca w koncepcie
(§11.2: „schematyczne, no sympy") księgowość polowa: lewa strona = bilans sektora Φ, prawa =
praca źródła konwersji. Kierunek znaku: J odbiera energię polu tam, gdzie φ̇ ≠ 0 — czyli
wyłącznie w poruszającej się warstwie (bulk: φ̇ = 0 ⇒ transfer 0 — zgodność z A-i LOCKED).

### §1.2 Ewaluacja na ścianie: transfer ∝ σ (FP2)

Na poruszającym się profilu φ̇ = −c·w′ i źródle zlokalizowanym na gradiencie J = j₀·w′:
T_area = −∫φ̇·J dx = c·j₀·∫(w′)²dx, a tożsamość BPS daje **∫(w′)²dx = ∫√(2V)dΦ = σ EXACT**
(zgodność z LOCKED FM P1 wartością σ = (2/3)√(λ/2)Φ₀³). Transfer na jednostkę powierzchni
frontu: **T_area = c·j₀·σ** — pierwsza polowa postać prawej strony EQ-5.

### §1.3 Amplituda źródła wyprowadzona (FP3-FP4)

Stałe ledgera połączone w jednostkach pola: **ΔV/σ = (3/8)·m_Φ EXACT**. Domknięcie ledgera
(popyt Ṁc² z marginalności, read-only; η = (t_*/t)² z regulatora Phase 2 — wynika,
nie wstawione):

> **j₀(t) = η(t)·ΔV/σ = (3/8)·m_Φ·(t_*/t)²** — amplituda S_creation w EQ-5
> wyrażona w jednostkach pola ([j₀] = masa = [m_Φ] ✓), bezparametrowo w (λ, Φ₀, G, c).

Weryfikacje: 4πR²·T_area − Ṁc² = 0 EXACT; Ė_wall = supply·(1 − η) = 0 przy t_*, > 0 dla
t > t_* (nadwyżka → kinetyka ściany — spójność z FM bez modyfikacji); wymiary
[ΔV·4πR²Ṙ] = [Ṁc²] = E/T EXACT (FP5).

### §1.4 Kwalifikacja (jawna)

Wyprowadzony jest most **bilansowy na poziomie pola** (tożsamość lokalna + lokalizacja +
amplituda). Pozostały rezyduał (ten sam co P2 §1.4, flagowany, poza pre-rejestrowanym
zakresem): mikroskopowa postać operatorowa J[Φ] (jaki funkcjonał pola kreuje solitony —
poziom fluktuacyjno-statystyczny, Appendix E territory). **F-BA-1 = DERIVED.**

## §2 — F-BA-4: jedyność w klasie (GAP-4)

Klasa zadeklarowana Phase 0 §8(g), tu ROZSZERZONA (nadzbiór — wzmocnienie, nie goalpost):
u = U·x/t, ρ = ρ₀·ξ^k/(Gt²), ξ = x/ct; warunki: ciągłość source-free + Euler z samograwitacją
+ depozycja marginalna flow-matched.

1. **Ciągłość (FP6):** jedyne U(k) = (k+2)/(k+3).
2. **Euler (FP7):** człony skalują się jak ξ¹ vs ξ^(k+1) ⇒ residual znika tożsamościowo
   tylko dla **k = 0** (ρ₀ > 0) — **jednorodność WYMUSZONA, nie narzucona** (zamyka caveat
   FM COR-1 „uniqueness nie wykazana" w obrębie klasy; jednocześnie wzmacnia samo A-ii).
3. **Para (U, ρ₀) (FP8):** k = 0 ⇒ U = 2/3 (ciągłość) ⇒ ρ₀ = 1/(6π) (Euler) — jedyne.
4. **Naddeterminacja domyka się EXACT (FP8):** marginalność (½v_c² = GM/ct) jest TRZECIM
   warunkiem na DWIE niewiadome — residual = 0 tożsamościowo. Spójność układu
   naddeterminowanego = nietrywialna treść (mogła obalić konfigurację — nie obaliła).
5. **Audyt odrzuceń (FP9):** ρ₀(U) = 3U(1−U)/(4π): U = 1 → próżnia; U > 1 → ρ₀ < 0
   niefizyczne; gałąź U = 1/3 (to samo ρ₀) odrzucona ciągłością (wymaga k = −3/2 ≠ 0).

**F-BA-4 = DERIVED_IN_CLASS.** Deklaracja niezmienna: jedyność GLOBALNA (poza klasą
self-similar) = poza zakresem (Phase 0); klasa jest jednak naturalna (wymuszona symetrią
sferyczną GAP-5 + brakiem skali zewnętrznej w R = ct).

## §3 — F-BA-5: sferyczność jako atraktor (GAP-5)

Toy zadeklarowany (Phase 0 §8(h) + per-use): nierelatywistyczna dynamika poprzeczna
thin-wall, masa powierzchniowa μ = σ/c², napęd ΔV jednorodny; operator średniej krzywizny
zlinearyzowany (standardowa geometria — assumption).

1. **Sztywność modów (FP10):** δ(2H) = (l−1)(l+2)·δr/R₀² EXACT — mody kształtu l ≥ 2 mają
   dodatnią siłę przywracającą; l = 1 zerowa (translacja); napęd ΔV shape-neutral
   (FP12: brak zależności kątowej — jednorodność per-area λ, Φ₀; FM COR-1).
2. **Dynamika na R₀ = ct (FP11):** b̈ = −(l−1)(l+2)·c²b/R₀² ⇒ równanie Eulera-Cauchy'ego;
   indicial p(p−1) + (l−1)(l+2) = 0; dyskryminanta 9 − 4l(l+1) < 0 dla wszystkich l ≥ 2
   (przy l = 2: −15; malejąca w l) ⇒ Re p = ½ ⇒ |b_l| ∝ t^(1/2)·oscylacje ⇒
   **a_l = b_l/R₀ ∝ t^(−1/2) → 0 — JEDNOLITE tempo zaniku dla wszystkich modów kształtu**
   (czysty wynik: szybsze mody nie wybierają skali — front gładnieje samoczynnie).
3. **l = 1 (FP12):** pierwiastki {0, 1} ⇒ a_1 ~ const = dryf środka (nie asferyczność
   kształtu w rzędzie liniowym) — deklarowane.

**Caveats (jawne):** (i) nierelatywistycznie — rosnący γ ściany (FM P1 FP5) zamraża dynamikę
poprzeczną w czasie własnym, co działa W TYM SAMYM kierunku (szybszy względny zanik) —
kierunkowo zgodne, NIE obliczone (INFORMATIONAL, nie użyte do oceny); (ii) rząd liniowy;
(iii) zewnętrzna strona warstwy pre-metryczna (m_eff² < 0) — zaburzenia nie mają tam kanału
wzmocnienia (brak modów propagujących; SCOPING Q3, INFORMATIONAL).

**F-BA-5 = DERIVED** (w zadeklarowanym przybliżeniu — dokładnie tym, które pre-rejestrowała
Phase 0 §8(h); żaden mod nie rośnie ⇒ klasa PARTIAL/mode-restricted niepotrzebna).

## §4 — GAP REGISTER po Phase 3 (moduł A komplet)

| GAP | Werdykt | Faza |
|---|---|---|
| GAP-1 EQ-5 field-level | **DERIVED** (residual operatorowy flagowany) | P3 |
| GAP-2 bottom-up J_source | **DERIVED** (regulator marginalnościowy) | P2 |
| GAP-3 A-iv | **SUPPORTED_PARTIAL** (kanał odrzutu; atraktor pyłowy chroni machinery) | P2 |
| GAP-4 uniqueness | **DERIVED_IN_CLASS** | P3 |
| GAP-5 sphericity | **DERIVED** (toy zadeklarowany) | P3 |

**F-BA-D strict reading:** 4/5 DERIVED(-IN-CLASS) + GAP-3 SUPPORTED_PARTIAL ⇒ formalnie
**BRIDGE_PARTIAL** (nie BRIDGE_COMPLETE). Decyzja progowa — czy SUPPORTED_PARTIAL z atraktorem
pyłowym spełnia próg „domknięcia luk A2" (FM Phase_FINAL §3) — należy WYŁĄCZNIE do użytkownika
w Phase FINAL (precedens FM (iv): wtedy strict reading utrzymano). **NO PR-022 w tej fazie**
(forbidden move #8). Moduł B: bez zmian (KB4/SIG-2 → Phase 4).

## §5 — Anti-Lakatos (Phase 3): COMPLIANT ✓

0/13 hardcoded; rozszerzenie klasy GAP-4 = nadzbiór (wzmocnienie jedyności, nie osłabienie
progu); F-BA-5 oceniony w przybliżeniu pre-deklarowanym Phase 0 §8(h) (γ-freezing kierunkowy
NIE użyty do oceny); F-BA-1 rezyduał operatorowy flagowany (nie ukryty); naprawa plumbingu
FP2 (jawna gałąź dodatnia √(2V) na dziedzinie, z weryfikacją branch² = 2V — bez zmiany progu);
LOCKED inputs read-only; circularity guard FP13 czysty; λ, Φ₀ symboliczne; 0 nowych stałych;
0 predecessor verdicts modified.
