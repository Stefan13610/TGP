---
title: "Phase 3 — F-FM-COR: corollaries — A-ii DERIVED_SELF_CONSISTENT (shell-map exact closure); C-2 DERIVED (dissolved); A2 PARTIAL (energy-ledger bridge + 3 gaps declared)"
type: phase_result
status: PHASE3_COMPLETE
phase: 3
cycle: op-frontier-microphysics-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'FM Phase 3'"
sympy_script: "[[./Phase3_sympy.py]]"
sympy_output: "[[./Phase3_sympy.txt]]"
sympy: "8/8 PASS; 0 hardcoded T_pass"
falsifier_resolved: "F-FM-COR = {COR-1 DERIVED_SELF_CONSISTENT, COR-2 DERIVED (dissolved), COR-3 PARTIAL}; PR-022 conditions: (i) SATISFIED (residual: A-iv + uniqueness caveat), (ii) DERIVED_SELF_CONSISTENT, (iii) SATISFIED, (iv) PARTIAL → user decision at Phase FINAL; NO PR-022 this phase"
anti_lakatos_lock: PRESERVED
---

# Phase 3 — korolaria: A-ii / C-2 / A2 (warunki ii-iv PR-022)

## §0 — Verdict at a glance

| Korolarium | Werdykt |
|---|---|
| **COR-1 (A-ii, jednorodność)** | **DERIVED_SELF_CONSISTENT** ⭐⭐ — jednorodność NIE jest narzucona: jest dokładnym wynikiem mapy powłok (marginalność + depozycja 2c/3 + transport przepływem), z pełnym domknięciem dynamicznym |
| **COR-2 (C-2, siła substratowa)** | **DERIVED (dissolved)** — wymagana siła ≡ 0; balans tożsamościowy przy ε = 2/3 |
| **COR-3 (A2, most Φ→materia)** | **PARTIAL** — most wyspecyfikowany EXACT na poziomie ledgera energii; 3 luki rezydualne zadeklarowane |

## §1 — COR-1: jednorodność WYPROWADZONA (mapa powłok Lagrange'a)

Powłoka kreowana w epoce t₀: depozycja na x₀ = ct₀ z prędkością 2c/3 (Phase 2), dalszy ruch
po wyprowadzonym przepływie: **x(t; t₀) = c·t₀^(1/3)·t^(2/3)**. Księgowość gęstości:

```
ρ(x,t) = Ṁ / (4πx² ∂x/∂t₀) = 1/(6πGt²)   — WOLNE od t₀ (więc od x)  [FP1 EXACT]
```

**Jednorodność jest dokładną konsekwencją** trzech wyprowadzonych elementów: (1) marginalność
⇒ Ṁ = 2c³/9G = const (FCR P3 + Phase 2); (2) depozycja flow-matched v_c = 2c/3 (Phase 2);
(3) transport wyprowadzonym przepływem (FCR P2). Izotropia per-area: parametry ściany (σ, ΔV)
jednorodne na froncie (te same λ, Φ₀) ⇒ kreacja na jednostkę powierzchni jednorodna.

**Domknięcie dynamiczne (FP2) — kluczowe wzmocnienie:** trajektorie powłok spełniają PEŁNE
samograwitujące równanie ruchu ẍ = −GM_enc/x² EXACT (residual 0); M_enc = M(t₀) zachowane
wzdłuż trajektorii; **∂x/∂t₀ > 0 — brak przecinania powłok (no caustics) ⇒ single-stream**
(retroaktywne wsparcie A-iv); domknięcie masy Σ powłok = M(t) = 2c³t/9G; powłoki wypełniają
przestrzeń do x → 0. Konfiguracja {u = (2/3)x/t, v_c = 2c/3, Ṁ const, ρ = 1/(6πGt²)} jest
**dokładnym rozwiązaniem pełnego układu**, nie tylko kinematyką.

**Caveats (FP3, jawne):** (a) jedyność rozwiązania samouzgodnionego nie wykazana (istnienie
+ dokładność TAK); (b) sferyczność locusu = założenie odziedziczone Phase 0 §8(e)(iv);
(c) jednorodność nie jest ścisłym atraktorem — zaburzenia rosną potęgowo t^(2/3) (pierwiastki
{2/3, −1}; klasa no-runaway R17) — dokładnie tak, jak rama sama używa (wzrost liniowy na
jednorodnym tle); bez sprzeczności wewnętrznej.

**Nota cyrkularności (uczciwa):** Phase 2 K2 używało A-ii do wyboru przepływu; Phase 3
wyprowadza jednorodność z depozycji + przepływu. To pętla punktu stałego: status A-ii =
**zweryfikowany samouzgodniony punkt stały (existence + exactness)**, nie liniowa dedukcja.
Upgrade: IMPOSED → DERIVED_SELF_CONSISTENT.

## §2 — COR-2: C-2 rozpuszczone (formalne zaksięgowanie)

Z EQ-1 (E_sol = E(⟨Φ⟩), dowolne E): F_substrat = −E′·∇⟨Φ⟩ = 0 w nasyconym bulku (∇⟨Φ⟩ = 0)
— model-independent. Residual tła Eulera res(ε) = (3ε−2)/9·x/t²: **res(2/3) = 0, Δ_bulk(2/3)
= 0 tożsamościowo**. Postulowana w FCR siła O(1)H²x okazuje się dokładnym zerem — caveat
**rozpuszczony, nie sfinansowany**. Warunkowość: funkcjonalna forma EQ-1 (concept LOCKED,
klasa rygoru §3.6.12) — zadeklarowana.

## §3 — COR-3: most A2 na poziomie ledgera energii (PARTIAL)

**Strona popytu (FP5):** marginalność ⇒ księgi mechaniczne = 0 przy wejściu ((1/2)v_c² −
GM/ct = 0 tożsamościowo) ⇒ substrat musi sfinansować WYŁĄCZNIE energię spoczynkową:
**popyt = Ṁc² = 2c⁵/9G = const EXACT**.

**Strona podaży (FP6):** konwersja offsetu próżni przy przejściu ściany:
**podaż = ΔV·4π(ct)²·c = πλΦ₀⁴c³t² ∝ t²**. Stosunek podaż/popyt = (t/t_*)² z progiem
**t_* = (√2/3√π)·c/(√(Gλ)·Φ₀²)** (jedyny dodatni; symbolicznie — forbidden move #4).

**R-FM-3 RESTRUCTURED (uczciwie, nie ukryte):** naiwne porównanie gęstości (ΔV const vs
ρ̄c² ∝ t⁻²) było mylące — po stronie całkowitych strumieni nie ma przeszkody późnoczasowej:
dla t ≥ t_* podaż ≥ popyt, nadwyżka (1 − (t_*/t)²) → 1 idzie w kinetykę ściany (spójne
z v → c, Phase 1 FP5). Pozostaje **deficyt wczesnoepokowy t < t_*** — uczciwe znalezisko
strukturalne (INFORMATIONAL; możliwy związek z epoką początkową; bez numeryki λ, Φ₀).

**Luki rezydualne (FP7, jawne — werdykt PARTIAL):**
1. Księgowość polowa EQ-5 (S_Φ vs S_matter w jednostkach pola) — schematyczna w concept
   paper (jego §11.2: no sympy) — luka odziedziczona;
2. Bottom-up rate nukleacji (funkcjonał J_source) niewyprowadzony — rate stoi TOP-DOWN
   z LOCKED marginalności;
3. A-iv (wejście monochromatyczne) — wsparte no-caustics (FP2), niewyprowadzone z mikrofizyki
   ściany.

## §4 — Status warunków PR-022 (po Phase 3)

| Warunek | Status |
|---|---|
| (i) tiebreaker derived | **SATISFIED** — warunkowość A-ii z Phase 2 wyprowadzona (§1); rezydualnie: A-iv (declared) + uniqueness caveat |
| (ii) A-ii | **DERIVED_SELF_CONSISTENT** (caveats §1 jawne) |
| (iii) C-2 | **SATISFIED** (dissolved §2) |
| (iv) A2 | **PARTIAL** (§3) — czy „bridge specified" na poziomie ledgera spełnia próg FCR Phase_FINAL §3 — **decyzja użytkownika w Phase FINAL**, nie oceniana pobłażliwie tutaj |

**NO PR-022 w tej fazie** (forbidden move #6: ścisłe odczytanie ALL conditions).

## §5 — Anti-Lakatos (Phase 3): COMPLIANT ✓

0/8 hardcoded; pętla samouzgodnienia COR-1 nazwana wprost (nie udaje liniowej dedukcji);
R-FM-3 zmierzone i zrestrukturyzowane z pełnym ujawnieniem (deficyt wczesnoepokowy zapisany,
nie zatuszowany); COR-3 oceniony PARTIAL wbrew pokusie domknięcia (decyzja progowa oddana
użytkownikowi); λ, Φ₀ symboliczne; G_obs nieobecne (FP8); 0 predecessor verdicts modified;
0 nowych stałych.
