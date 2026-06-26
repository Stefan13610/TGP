---
title: "Phase 0 — pre-registration: op-PR004-SPARC-fit-execution — wykonanie LOCKED falsyfikatora PR-004 (SPARC 175: χ²_red TGP-Newton-baryon vs MOND simple)"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-12
cycle: op-PR004-SPARC-fit-execution-2026-06-12
authorization: "User 2026-06-12: 'do rotacji galaktyk podchodziłem w TGP kilkukrotnie i za każdym razem się nie udawało, ale spróbujmy' → aktywacja wykonania PR-004 (LOCKED-PENDING-FIT od 2026-05-13)"
methodology_binding: "PR-004 decision rule IMMUTABLE (PRE_REGISTERED_FALSIFIERS.md §PR-004); CALIBRATION §3.6; pipeline poniżej LOCKED PRZED pobraniem danych i jakimkolwiek χ²"
anti_lakatos_lock: PRESERVED
---

# Phase 0 — pipeline fitu (LOCKED przed danymi)

## §1 — Co jest testowane (z LOCKED PR-004, bez zmian)

**Model TGP (per retrofit op-L01-N3 2026-05-13, L2 chain LOCKED):** w reżimie galaktycznym
(v²/c² ~ 10⁻⁷) TGP z g_eff[Φ̄ ≈ Φ₀] redukuje się do Newtona z wyłącznie barionowym źródłem:
v_TGP²(R) = v_bar²(R) = V_gas|V_gas| + Υ_d·V_disk|V_disk| + Υ_b·V_bul². **Zero wolnych
parametrów globalnych; zero ρ_DM (S05 enforcement — forbidden direction).**

**Benchmark (PR-004):** MOND simple: g = ν(y)·g_bar, ν(y) = 1/2 + √(1/4 + 1/y),
y = g_bar/a₀, a₀ = 1.2×10⁻¹⁰ m/s² (fixed, Lelli+2017). Ten sam v_bar, te same dane.

**Decision rule (IMMUTABLE, verbatim PR-004):** jeśli χ²_red(TGP) > χ²_red(MOND simple)
z 5σ na próbie 175 galaktyk ⇒ mechanizm rotacyjny TGP z g_eff[Φ̄] **insufficient** ⇒
(a) potrzebna osobna ρ_DM (S05 violated; framework revision) lub (b) dedykowany mechanizm —
„framework needs structural amendment, NOT continued recovery".

## §2 — Pipeline (LOCKED; wszystkie wybory PRZED jakimkolwiek wynikiem)

1. **Dane:** SPARC (Lelli+2016): MassModels_Lelli2016c.mrt (R [kpc], V_obs ± e_V, V_gas,
   V_disk, V_bul [km/s]) + tabela własności (jakość Q, inklinacja). Źródło:
   astroweb.cwru.edu/SPARC (LITERATURE_ANCHORED). V_obs używane jak opublikowane
   (inklinacja już skorygowana przez SPARC).
2. **Υ (mass-to-light, 3.6 μm) FIXED (fiducial Lelli+2016/2017):** Υ_disk = 0.5,
   Υ_bulge = 0.7 [M☉/L☉]. ZERO dopasowywanych parametrów w obu modelach (uczciwe
   porównanie bezparametrowe; standard literaturowy benchmarku ≈ 2.0).
3. **χ²_red:** per galaktyka: χ²_g = Σ_i [(V_obs,i − V_model,i)/e_V,i]² / N_g;
   agregat: (a) GLOBAL = Σχ²_g·N_g / ΣN_g (punktowo), (b) MEDIAN per-galaxy (odporny).
   Raportowane oba; decyzyjny = GLOBAL (deklaracja).
4. **Operacjonalizacja 5σ (deklarowana tu, bo PR-004 jej nie precyzuje):** test sparowany
   per-galaktyka: d_g = χ²_g(TGP) − χ²_g(MOND); kryterium: mean(d)/SEM(d) > 5
   (+ sanity: bootstrap 10⁴ percentyl). Pomocniczo: Δχ²_total vs √(2·dof).
5. **Próba:** primary = pełna próba z kompletem kolumn; secondary (allowed direction
   PR-004) = podpróba Q = 1+2. Punkty z e_V ≤ 0 lub V_obs ≤ 0 odrzucone (deklarowane).
6. **Sanity anchor (nie-werdyktowy):** χ²_red(MOND) musi wyjść rzędu literaturowego ~2
   (jeśli rażąco inne ⇒ błąd pipeline'u, nie fizyka — naprawa plumbingu dozwolona
   PRZED odczytem werdyktu TGP, klasa „better seeds").
7. **Circularity guard:** żadna wielkość TGP nie jest dopasowywana do V_obs; a₀ należy
   wyłącznie do benchmarku MOND.

## §3 — Forbidden moves (inherit PR-004 + wykonawcze)

1. Kolumna ρ_DM / halo (S05) — zakaz. 2. Tuning Υ, a₀, ani żadnego parametru po obejrzeniu
χ² — zakaz (Υ fixed §2.2). 3. Zmiana reguły decyzyjnej / progu 5σ — zakaz (IMMUTABLE).
4. Selekcja galaktyk po wyniku — zakaz (podpróba wyłącznie po Q, zadeklarowana z góry).
5. Hardcoded T_pass — zakaz. 6. OR-clause rozszerzeń — zakaz. 7. Ukrycie wyniku FAIL —
zakaz (wynik idzie do PR-004 status update niezależnie od kierunku).

## §4 — Anticipated outcome (INFORMATIONAL, zapisane PRZED rachunkiem)

Newton+baryony na SPARC historycznie daje χ²_red ≫ MOND (to jest obserwacyjny problem
ciemnej materii); użytkownik: wcześniejsze podejścia TGP do rotacji „za każdym razem się
nie udawało". **Oczekiwany kierunek: TRIGGERED — mechanizm g_eff[Φ̄ ≈ Φ₀] insufficient.**
Zapisujemy to, by werdykt FAIL nie był „niespodzianką" łagodzoną ex post; PASS byłby
anomalią wymagającą audytu pipeline'u. Wartość naukowa wyniku FAIL: twardy constraint
strukturalny na sektor grawitacyjny TGP w reżimie niskich przyspieszeń (analog F8-klasy).

## §5 — Stałe i anchory

a₀ = 1.2×10⁻¹⁰ m/s² (MOND benchmark only); Υ_d = 0.5, Υ_b = 0.7 (LITERATURE_ANCHORED
fiducial); G, c — standard; dane SPARC = LITERATURE_ANCHORED. Nowe stałe fundamentalne: 0.
