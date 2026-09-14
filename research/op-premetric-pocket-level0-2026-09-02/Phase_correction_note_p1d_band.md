---
title: "Correction note (PRZED użyciem skorygowanego wyniku): P1d — pomiar szerokości pasa przedmetrycznego bez interpolacji podsiatkowej = błąd implementacji TESTU (dyskretyzacja ~dx≈2% > próg 1%); fizyka zbieżna (σ_w zgodne do 6 cyfr, ΔΦ_min/Φ_bar=4.5e−4)"
date: 2026-09-02
type: correction-note
tgp_owner: research/op-premetric-pocket-level0-2026-09-02
status: DOCUMENTED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase1_output_pre_correction.txt]]"
---

# Korekta 1 — implementacja pomiaru szerokości pasa (P1d)

**Pierwotny wynik (zachowany, niezamazany):**
`Phase1_output_pre_correction.txt` — P1d FAIL: Δszer(Φ<ε) = 1.551%
przy progu <1% (dx=0.05 vs 0.025); jednocześnie σ_w zgodne do 6 cyfr
(0.429765 vs 0.429771) i ΔΦ_min/Φ_bar = 4.5e−4.

**Diagnoza (błąd implementacji ZNALEZIONY):** funkcja `band()`
w `Phase1_gate.py` liczyła szerokość pasa jako `(i1−i0)·dx` po węzłach,
z martwym/wadliwym członem interpolacyjnym (`Phi[i0−1+1]` ≡ `Phi[i0]`
— wyrażenie nie interpoluje niczego; brzeg prawy nieinterpolowany
w ogóle). Pomiar po węzłach ma błąd dyskretyzacji O(dx) ≈ 0.05/2.39
≈ 2% — czyli sam TEST nie może spełnić progu 1% niezależnie od fizyki.
To odpowiednik incydentu QF-4c poprzednika (kryterium porównane
z artefaktem estymatora, nie z fizyką).

**Korekta (wyłącznie kod pomiaru, zero zmian progów/kryteriów/modelu):**
szerokość pasa = x_prawy − x_lewy, gdzie oba punkty przecięcia Φ=próg
wyznaczane interpolacją liniową między sąsiednimi węzłami na OBU
brzegach pasa. Kryterium LOCKa (<1% między dx a dx/2) NIETKNIĘTE.

Zapisano PRZED ponownym uruchomieniem bramki i przed jakimkolwiek
biegiem Phase 2.
