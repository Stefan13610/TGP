---
title: "ANALIZA N1 — pochodzenie faz ogona z dodatekH lin. 1126–1129: ZNALEZIONE (p131, układ eq:J-ode Formulacji B, nie układ O-L5)"
date: 2026-08-23
type: analiza-post-close
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: ROZSTRZYGNIĘTE
verdict: "Fazy δ_e=−81.4°, δ_μ=+38.6°, Δ(e→μ)=120.01° REPRODUKUJĄ SIĘ CO DO SETNEJ STOPNIA — ale z układu eq:J-ode (F(g)·g″+(2/r)g′=V′, F=1+2α_eff·ln g; forma logarytmiczna typu F-S), zaimplementowanego w _archiwum p131_eta_refinement.py, NIE z układu potęgowego K=g^(2α), który opisuje status_map O-L5 / dodatekH jako 'kanoniczną Form A'. Rozjazd dokumentacyjny, nie fabrykacja."
related:
  - "[[PhaseA_report.md]]"
  - "[[NEEDS.md]]"
  - "[[N1_provenance_check.py]]"
  - "[[../../core/formalizm/dodatekH_lancuch_wyprowadzen.tex]]"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]"
---

# N1 — pochodzenie faz w dodatekH: ROZSTRZYGNIĘTE

**Charakter dokumentu:** POST-CLOSE (cykl CLOSED-GATE-STOP; to śledztwo NEEDS N1
autoryzowane przez użytkownika 2026-08-23, nie kontynuacja zalockowanych faz).

## 1. Znalezisko

Wartości z dodatekH lin. 1126–1129 pochodzą ze skryptu
`_archiwum/scripts_exploratory/advanced/p131_eta_refinement.py` (sesja v42+;
łańcuch p127→p128→p131→p134f: p134f cytuje wprost „From p131:
Delta(e->mu) = 120.01 deg"). Skrypt rozwiązuje układ:

```
(P131)  F(g)·g″ + (2/r)·g′ = V′(g)
        F(g) = 1 + 2·α_eff(g)·ln(g),  α_eff(g) = 2/(1+η_K(g−1)²)
        V′(g) = g²(1−g)
```

To jest **struktura eq:J-ode** (sek08b lin. 38–40: `f(g)g″+(2/r)g′=V′(g)`,
f=1+2α·ln g — równanie solitonu dod. app:ogon-masy, **Formulacja B**,
z runningiem α_eff wstawionym pod logarytm) — a **NIE** układ potęgowy
`g″+2g′/r+(α/g)g′²=g²(1−g)` (K=g^(2α)), który status_map O-L5 i narracja
dodatekH opisują jako „kanoniczną Form A" i który audytowała bramka A.

## 2. Reprodukcja numeryczna ([[N1_provenance_check.py]] → [[N1_provenance_output.txt]])

Układ P131, η_K=181/15, g₀_e=0.90548, g₀_μ=φ·g₀_e, okno fitu [120,260],
konwencja atan2(C,B) — identyczne z p131:

| wielkość | zmierzone (P131) | dodatekH | audytowane M2 (kontrola krzyżowa) |
|---|---|---|---|
| δ_e | **−81.43°** | −81.4° ✓ | −75.34° |
| δ_μ | **+38.58°** | +38.6° ✓ | +97.20° |
| Δ(e→μ) | **120.01°** | 120.01° ✓ (co do setnej) | ~172° |
| A_μ(g₀=1.4550) | 0.3584 | 0.3861 (mapa p127) | 0.6159 |

Baseline p131 (g₀_e=0.89926559): δ_e=−80.85°, δ_μ=+39.67°, Δ=120.52° —
zgodne z mapą fazową p127–p128 (−81.1°, +43.8° przy nieco innych parametrach).
A_μ: 0.358 vs 0.3861 mapy (7%; mapa p127 mogła używać innego okna/η — fazy
się zgadzają, amplituda z dokładnością do konwencji fitu).

**Czułość na η_K:** Δ(e→μ)=120.51–120.52° dla η∈[12,14] — wynik „120°" jest
**niewrażliwy na strojenie η_K**; to własność formy logarytmicznej F, nie fitu.

## 3. Wnioski

1. **Fazy nie są sfabrykowane** — są w pełni reprodukowalne z konkretnego,
   zachowanego w archiwum skryptu.
2. **Rozjazd jest dokumentacyjny:** dodatekH/status_map przypisują te liczby
   narracyjnie układowi potęgowemu („kanoniczna Form A, K=g⁴"), podczas gdy
   faktycznie pochodzą z układu logarytmicznego eq:J-ode (Formulacja B
   z runningiem). Werdykt A3 bramki pozostaje poprawny: **z układu, który
   dokumentacja deklaruje, te fazy nie wychodzą** (Δ≈164–172°).
3. Uwaga metodologiczna z p131 (do rejestru): η_K było tam **bisekowane
   w [12,18] do trafienia r₃₁=3477** (jawny fit do celu), a g₀τ wybierane
   skanem maksymalizującym r₃₁; późniejsza derywacja analityczna η_K=181/15
   (p139–p140, 2 ppm) jest od tego niezależna, ale Δ=120° i tak nie zależy
   od η (patrz §2).
4. **Δ(e→μ)=120.01° to twarda własność formy logarytmicznej** — po korekcie
   dokumentacji (przypisanie do właściwego układu) wynik można przywrócić
   do obiegu, ze statusem zależnym od statusu samej Formulacji B (patrz
   [[ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]] — dwie gałęzie znaku W).

## 4. Proponowane korekty rdzenia (USER-GATED, nie wykonane)

- dodatekH lin. ~1126–1129 + p131/p134: dopisek, że fazy pochodzą z układu
  eq:J-ode (F=1+2α_eff ln g), z odsyłaczem do p131 i tej analizy; usunięcie
  sugestii, że to własność układu potęgowego O-L5.
- status_map O-L5: analogiczna flaga przy „Fazy ogona".
- Mapa p127 (A_μ=0.3861): flaga [WINDOW-DEPENDENT] albo re-derywacja.

## 5. POST-SCRIPTUM 2026-08-31 — potwierdzenie niezależne + realizacja korekt (user-gate)

- **Niezależne potwierdzenie (op-bath-two-sectors, Phase 1, fazy ZMIERZONE):**
  Δ_ML(e→μ)=120.01° zreprodukowane na fazach mierzonych niezależnie
  (okna [50,150] i [120,260], R=300/450, R-kontrola |Δδ|=0.0000°);
  drabina minimów 2π±5% PASS 6/6 (odchył 0.96–3.34%, malejący z d).
  Wniosek tej analizy (fazy = własność układu logarytmicznego eq:J-ode)
  potwierdzony pomiarowo.
- **Korekty rdzenia z §4: WYKONANE 2026-08-31** (autoryzacja użytkownika
  „działaj z krokami 1,2,3"): dopisek o pochodzeniu faz + flaga przy
  statusie O-K1 w dodatekH; flaga [2026-08-31] w wierszu Koide Q=3/2
  w status_map. Mapa p127: pozostaje z uwagą w dopisku dodatekH
  (A_μ=0.3861 sprzeczne z atail_functional).
