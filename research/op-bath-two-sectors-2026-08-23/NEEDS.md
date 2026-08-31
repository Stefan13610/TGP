---
title: "NEEDS — propozycje po op-bath-two-sectors (user-gated; Q1-INCONCLUSIVE, Q2-FAIL)"
date: 2026-08-29
type: needs
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: "EXECUTED 2026-08-31 — user-gate przyznany (autoryzacja „działaj z krokami 1,2,3", krok 2); N1–N4 wykonane, log przy pozycjach"
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# NEEDS — wnioski cyklu (wszystkie addytywne, user-gated)

Zgodnie z drzewem LOCKa §6 i werdyktami: **Q1-INCONCLUSIVE**, **Q2-FAIL**.

> **LOG WYKONANIA (2026-08-31, user-gate „działaj z krokami 1,2,3" krok 2):**
> - **N1 ✅ EXECUTED** — remark `rem:W-sign-axiomatic` w sek08a (przy
>   prop:field-eq-from-action): lokalizacja rozdwojenia (ANALIZA_N2) +
>   wynik Q2 (ω²_min(d) dodatnie, rosnące z gęstością) + nota poziomu 0
>   (spinodala MFT); wniosek „znak W = otwarty problem aksjomatyczny".
> - **N2 ✅ EXECUTED (oba miejsca z propozycji)** — dopisek w Limitations
>   `tgp_lepton_masses.tex` ORAZ nota metodologiczna 2026-08-31 w README
>   cyklu #63 (`op-nonlinear-charge-constraint`): pytanie o stabilizację
>   kąpielą OTWARTE + ustalenie metodologiczne (wariant komórkowy
>   zero-flux NIEZDOLNY rozstrzygnąć; następca: tło samouzgodnione
>   + periodyczne BC / Bloch — realizowany jako
>   `op-bloch-chain-stability-2026-08-31`).
> - **N3 ✅ EXECUTED** — flaga TAU-NEAR-THRESHOLD w dodatekH (dopisek
>   „istnienie τ": 1.9% pod progiem 8/5, wrażliwość ±2%) + flaga
>   w status_map przy Koide Q=3/2.
> - **N4 ✅ EXECUTED** — dopisek w dodatekH przy fazach: Δ_ML(e→μ)=120.01°
>   zmierzone niezależnie potwierdza ANALIZA_N1 co do setnej.
> - Build gate: `build_gate_2026-08-31.log` — main.pdf 557 str., zero
>   NOWYCH undefined refs (porównane z main_build.log/main_build_62.log).

## N1 (HIGH, z Q2-FAIL) — `sek08a` okolice prop:field-eq-from-action / Noty kanonicznej: znak W = problem aksjomatyczny

Dopisek datowany (addytywny) do miejsca, gdzie N2 zlokalizował rozdwojenie
znaku członu źródłowego W:

> Test derywacyjny hipotezy „dwa sektory jednej akcji" (izolowany obiekt
> vs obiekt w kąpieli; research/op-bath-two-sectors-2026-08-23, Phase 3):
> w akcji gałęzi STABILNEJ (eq:field-eq-reproduced, m²=+γ) samouzgodnione
> tła sieci źródeł d∈{8,6,4,2} dają ω²_min(d) = +1.34, +1.57, +1.88, +2.47
> (zbieżne; kontrola d=∞ odtwarza Yukawę m²=γ±0.00%) — wszystkie DODATNIE
> i ROSNĄCE z gęstością; odpowiedź statyczna monotoniczna dla wszystkich d.
> Znak tachionowy NIE emerguje z gęstości. **Wybór znaku W pozostaje
> otwartym problemem aksjomatycznym (decyzja ontologiczna, nie numeryka).**
> Zbadana klasa: komórka 1D periodyczna, źródła gaussowskie σ_s=0.5,
> q∈{0.3,1.0} [INPUT].

## N2 (MED, z Q1-INCONCLUSIVE) — Limitations (paper_lepton_masses lub notatka przy #63): wariant komórkowy nie rozstrzyga stabilności w kąpieli

Pytanie N3 poprzedniego cyklu (czy skończona gęstość sąsiadów stabilizuje
mod runaway) pozostaje OTWARTE, ale z ustaleniem metodologicznym:

- Komórka WS z obciętym profilem izolowanym + δg_bath i zero-flux BC jest
  NIEZDOLNA rozstrzygnąć Q1: ω²_min podąża za rozmiarem komórki (R<π: +,
  R>π: −), 0/15 punktów zbieżnych; komparator BEZ kąpieli daje niemal
  identyczne spektra i breakdowny (t*/t*_izo 0.03–1.81 również przy
  amp=0) — artefakt niestacjonarnego obcięcia dominuje nad efektem kąpieli
  (|c_bath| ≤ 0.45 przy zalockowanych d*).
- Ewentualny następca (do osobnego Phase 0, jeśli user zdecyduje):
  (a) tło SAMOUZGODNIONE w komórce (relaksacja przed linearyzacją),
  (b) periodyczne BC zamiast zero-flux, (c) analiza pasmowa (Bloch)
  zamiast pojedynczej komórki — czyli przeniesienie konstrukcji, która
  w Q2 dała czystą zbieżność, na sektor dynamiczny M0-f_eps.

## N3 (LOW, deskryptywne) — flaga TAU-NEAR-THRESHOLD

τ w M-P (g₀=1.5696 z Q_K=3/2 [INPUT]) leży 1.9% pod progiem kolapsu 8/5:
fit [120,260] PASS (R²=0.99972), fit bliski [50,150] FAIL (R²=0.99878),
wrażliwość +2% → KOLAPS. Jeśli gdziekolwiek w core/papers τ jest używany
jako regularny soliton M-P — dopisać flagę near-threshold z cytatem cyklu.

## N4 (LOW, higiena) — pozytywny produkt uboczny Phase 1

Δ_ML(e→μ)=120.01° zmierzone niezależnie (okna [50,150] i [120,260],
R=300/450) potwierdza N1 co do setnej stopnia; drabina minimów 2π±5%
(odchył 0.96–3.34%, malejący z d) potwierdza model referencyjny locka
oscylacyjnego na PARACH gatunków. Jeśli user zechce, można to dopisać
jako uwagę potwierdzającą przy ANALIZA_N1 poprzedniego cyklu.
