---
title: "Phase0_balance (LOCK) — op-oscillon-small-amplitude: właściwy test predykcji P1b (ω₂=−139/24<0): czy gałąź zdrowa ma oscylony MAŁEJ amplitudy?"
date: 2026-09-14
type: phase0-balance
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
status: LOCKED
anti_lakatos_lock: ACTIVE
related:
  - "[[HANDOFF_PROMPT.md]]"
  - "[[README.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
---

# Phase 0 — LOCK (zero obliczeń przed zapisem tego pliku)

**Autoryzacja:** user-gate 2026-09-14 — wybór „N1: cykl małych amplitud"
z [[../op-r3-stationary-states-2026-09-14/NEEDS.md]] (N1, priorytet).

## 0. Kontekst i pytanie

Poprzednik (`op-r3-stationary-states`, CLOSED: Q-E-INCONCLUSIVE) pre-rejestrował
predykcję analityczną **P1b (Lindstedt–Poincaré): ω₂ = −139/24 ≈ −5.792 < 0**
(miękka nieliniowość ⟹ oscylony małej amplitudy oczekiwane, ω(a)<m=1),
ale jego zamrożona rodzina startów (|a|≥0.15) NIE sondowała domeny predykcji.
Ten cykl jest **konfrontacją predykcji P1b w jej naturalnej domenie** (a≤0.10,
t_max=10⁴). Wynik negatywny = falsyfikacja P1b w klasie zbadanej (pełnoprawny wynik).

**Q-G (centralne):** czy dynamika 2. rzędu gałęzi zdrowej (M,𝒦,𝒰 — CYTAT niżej)
posiada długożyciowe oscylony małej amplitudy — stany z E_core utrzymanym
≥100 T₀ i częstością dominującą **ω_peak < 1** (poniżej progu kontinuum;
kryterium „rdzeń związany" z P1a poprzednika: κ²=ω²−1<0 ⟺ profil zlokalizowany)?

**Zakres:** istnienie. ZAKAZ claimów o stosunkach mas leptonów i o dyskretności
rodzin (dyskretność = osobny przyszły cykl, warunkowy na Q-G-PASS).

## 1. Model (formy FROZEN — CYTAT, zero modyfikacji)

Z zamkniętego `op-action-audit-spectrum-insert-2026-09-13/Phase1_output.txt`
(Q-D1-PASS; konwencja |g^tt|, user-gate N2 2026-09-14); K_geo=γ=c₀=1:

- M(ψ)=ψ⁶/(4−3ψ)², M′(ψ)=12ψ⁵(2−ψ)/(4−3ψ)³, 𝒦(ψ)=ψ⁴, 𝒰(ψ)=ψ⁴/4−ψ³/3
- π=Mψ̇; EOM: Mψ̈+½M′ψ̇² = (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′
- Dziedzina ψ∈(0,4/3); próżnia ψ*=1; m²=𝒰″(1)=1; T₀=2π; 100 T₀=628.319
- Gęstość energii potencjalnej w ewaluatorze: tożsamość bez kancelacji
  𝒰(ψ)−𝒰(1)=(ψ−1)²(3ψ²+2ψ+1)/12 (korekta 2 poprzednika — dziedziczona).

Silnik: adaptacja `../op-r3-stationary-states-2026-09-14/engine_core.py`
(uogólniony Störmer–Verlet na (ψ,π), iteracje punktu stałego do stagnacji
maszynowej; siatka przesunięta r_i=(i+½)h; ψ′(0)=0, l'Hôpital w r=0) —
kopia do własnego katalogu z cytatem w method_decisions; zmiany TYLKO
parametryczne (R, t_max, sponge okno) + staging (§3), zero zmian schematu.

## 2. Predykcja pre-rejestrowana (dziedziczona P1b + mapowanie)

- **P1b (CYTAT, nienaruszalna):** ω₂=−139/24<0. Kwalitatywnie: oscylacje
  o skończonej amplitudzie mają ω<1 (spadek częstości pod próg kontinuum).
- **Mapa ilościowa (MIĘKKA, deskryptywna — nie werdykt):** dla modu
  jednorodnego ω(a)≈1+ω₂a²: a=0.02→0.9977; 0.05→0.9855; 0.08→0.9629;
  0.10→0.9421. Dla startów zlokalizowanych współczynnik może się różnić
  (inna geometria) — pre-rejestrowany jest ZNAK (ω_peak<1) i monotonia w a;
  wartości raportowane vs mapa bez progu PASS/FAIL.
- **Phase 1 (analityczna, przed numeryką):** P1a′ — weryfikacja cytatów form
  i tożsamości π-formulacji (sympy, simplify=0, gate 1e−12 w ψ∈{0.9,1,1.1});
  P1b′ — CYTAT ω₂ z poprzednika (bez ponownego wyprowadzania; wpis wartości
  i źródła); tabela przewidywanych ω(a) j.w. zapisana PRZED Phase 3.

## 3. Protokół numeryczny (FROZEN)

**Siatka/parametry:** radialna przesunięta r_i=(i+½)h; **h=0.05 (główna)**,
h=0.025 (potwierdzenia wg reguły niżej); **R=400**; sponge smootherstep
γ₀=1.0 na r∈[320,400]; **dt=0.005** (dt/2=0.0025 w potwierdzeniach);
**t_max=10000**; E_core: r≤80; E_ref=E_core(t=50); zapis ψ(0,t) i E_core(t)
co dt_out=0.1; pełny profil co 100 j.cz.; checkpointy npz co 500 j.cz.

**Rodzina startów (FROZEN, deterministyczna, π₀=0, brak seeda):**
ψ(r,0)=1+a·exp(−r²/(2σ²)) dla **a∈{0.02,0.05,0.08,0.10} × σ∈{3,6,10}**
(12 startów) + próżnia kontrolna ψ≡1.

**Staging (FROZEN — uczciwa oszczędność budżetu, nie zmiana kryteriów):**
- **Etap A:** wszystkie 13 biegów, h=0.05, do t=2000.
- **Triage (reguła zamrożona):** do Etapu B (kontynuacja do t_max=10000)
  przechodzi bieg ŻYWY w t=2000: **E_core(2000)≥0.2·E_ref** ORAZ zero zdarzeń
  brzegowych. Biegi martwe klasyfikowane z Etapu A (ich los jest już
  rozstrzygnięty: energia w tym układzie tylko odpływa przez sponge).
- **Etap B:** kontynuacja żywych do t_max=10000 (z checkpointów).

**Detektor OSCILLON (FROZEN — litera identyczna z poprzednikiem + kryterium ω):**
1. E_core(t) ≥ 0.5·E_ref nieprzerwanie przez ≥100 T₀ (628.319 j.cz.) w oknie
   zaczynającym się ≥ t=50;
2. ≥50 przejść ψ(0,t) przez 1 w tym oknie;
3. **ω_peak ≤ 0.99** (FFT ψ(0,t), okno Hanna, segment stabilny ≥1000 j.cz.;
   0.99 = próg kontinuum 1 minus margines 2Δω przy Δω≈0.006);
4. potwierdzenie: h=0.025 ORAZ dt/2 — czas życia ±10%, ω_peak ±2%.

**Kategorie (FROZEN):**
- **OSCILLON** — detektor j.w. w całości.
- **OSCILLON-WEAK** (deskryptywna, NIE pozytyw Q-G): plateau E_core poniżej
  0.5·E_ref, ale quasi-stacjonarne (|ΔE_core|/E_core ≤ 1e−3 na ostatnich
  1000 j.cz.) ORAZ ω_peak≤0.99; raportowana obowiązkowo, wymaga tych samych
  potwierdzeń co OSCILLON, NIE zmienia werdyktu Q-G.
- **RADIATED** — E_core(t_end) < 0.05·E_ref, spadek bez plateau; τ = czas
  zejścia pod 0.5·E_ref.
- **COLLAPSE** (nadkategoria — realizacja propozycji N4 poprzednika w NOWYM
  LOCKu): wejście w pas graniczny (ψ>4/3−1e−6 lub ψ<1e−6) LUB
  niefinityczność, jeżeli oba zdarzenia typu kolapsowego mieszczą się
  w oknie ≤1 j.cz. między siatkami; podtyp (BOUNDARY-UPPER/-LOWER/NONFINITE)
  raportowany deskryptywnie.
- **INCONCLUSIVE-RUN** — kategoria niezbieżna między siatkami (poza regułą
  nadkategorii COLLAPSE).

**Potwierdzenia h=0.025 (reguła FROZEN):** (i) każdy kandydat
OSCILLON/OSCILLON-WEAK — pełne potwierdzenie; (ii) obowiązkowo start
(a=0.05, σ=6) niezależnie od klasy (kontrola zbieżności negatywu);
(iii) dt/2 dla kandydatów i przy zdarzeniach COLLAPSE.

## 4. Phase 2 — bramka maszynerii (FROZEN; FAIL ⟹ STOP)

- P2a: próżnia ψ≡1, sponge ON, 100 T₀, obie siatki: ‖ψ−1‖∞ ≤ 1e−10.
- P2b: regresja vs poprzednik — start (a=+0.15, σ=3) na R=400/h=0.05 do
  t=300: τ (zejście pod 0.5·E_ref) zgodne z poprzednikiem (206.9) do ±5%
  (różnica dopuszczalna od R=200→400 i innego okna sponge).
- P2c: dryf energii ≤1e−6/100 T₀ (start a=0.05 σ=6, do t=700, obie siatki);
  odbicie sponge ≤1e−3 (test różnicowy vs R=800 do t=350).

## 5. Werdykty (litera; INCONCLUSIVE ≠ pozytyw)

- **Q-G-PASS:** ≥1 start z klasą OSCILLON potwierdzoną (h/2 i dt/2).
- **Q-G-FAIL:** WSZYSTKIE 12 startów RADIATED zbieżnie (kategoria i τ ±10%
  na h i h/2 w zakresie potwierdzeń §3) ORAZ zero OSCILLON-WEAK.
  = **falsyfikacja P1b w klasie zbadanej** (a≤0.10, t≤10⁴).
- **Q-G-INCONCLUSIVE:** każdy inny rozkład kategorii.

## 6. Drzewo decyzyjne (pre-rejestrowane)

- **Q-G-PASS** → hipoteza ratunkowa MA nośnik; następny cykl (osobny LOCK,
  user-gate): dyskretność/rodziny (analog Q-F) + konfrontacja ω(a) z LP;
  kandydat dopisku core (user-gate).
- **Q-G-FAIL** → P1b sfalsyfikowana w domenie naturalnej; łącznie
  z Q-E-INCONCLUSIVE poprzednika = brak nośnika oscylonowego w gałęzi
  zdrowej w OBU klasach amplitud → user-gate: status łańcucha leptonowego
  (dopisek core rem:psi-EOM-R3-branch-status).
- **Q-G-INCONCLUSIVE** → NEEDS metodologiczny (t_max/staging/detektor);
  jeżeli dominują COLLAPSE także przy a≤0.10 — eskalacja do wyniku
  cyklu-bliźniaka `op-collapse-matter-source-2026-09-14` (Q-H2).

## 7. Forbidden moves (egzekwowane)

Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko odczyt;
formy M,𝒦,𝒰 bez modyfikacji (cytaty); ZAKAZ podłóg/barier; rodzina startów/
detektor/progi/triage/sponge niezmienialne po pierwszym biegu produkcyjnym;
INCONCLUSIVE nie reinterpretowane; mapa ilościowa ω(a) pozostaje miękka
(zakaz podnoszenia jej do rangi werdyktu POST HOC); ZAKAZ claimów o masach
leptonów; correction note wyłącznie dla błędu implementacji, PRZED użyciem
wyniku, pierwotne outputy zachowane.

## 8. Deliverables

`Phase_method_decisions.md` (FROZEN przed kodem) · `engine_core.py` (kopia
z cytatem) · `Phase1_analytic.py`+output · `Phase2_gate.py`+output ·
`Phase3_evolve.py`+output + `Phase3_results/` (json+npz per bieg,
verdict.json) · `Phase_FINAL_close.md` · `NEEDS.md` · dopis logu `README.md`.
Cykl bez FINAL+NEEDS+README NIE jest zakończony.
