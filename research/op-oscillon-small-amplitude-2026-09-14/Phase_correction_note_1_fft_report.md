---
title: "Phase_correction_note_1 — korekta warstwy ANALIZY/RAPORTU: (a) fft_peak zwracał pozorny pik (pierwszy bin ω≥0.1, moc=0) dla sygnału tożsamościowo zerowego (kontrola próżni); (b) addendum deskryptywne: kolumna ω_desc dla obowiązkowej tabeli ω_peak vs mapa LP. Detektor/progi/kategorie/trajektorie NIETKNIĘTE."
date: 2026-09-15
type: phase-correction-note
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
status: APPLIED-BEFORE-USE
related:
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_output_interim_stageA.txt]]"
---

# Correction note 1 (zapisana PRZED użyciem poprawionych wyników)

## (a) Błąd implementacji — fft_peak na sygnale zerowym

W `Phase3_evolve.py::fft_peak` dla sygnału tożsamościowo zerowego
(kontrola próżni: ψ(0,t)≡1 dokładnie ⟹ PSD ≡ 0) `argmax` zwracał
pierwszy bin ω≥0.1 jako „pik" (moc = 0), co w tabeli pośredniej
(`Phase3_output_interim_stageA.txt`, ZACHOWANY) wypisało dla vac
pozorne ω=0.1031±0.0032. Pik o mocy zerowej nie jest pikiem —
błąd wyłącznie warstwy raportującej.

**Korekta:** `fft_peak` zwraca None, gdy moc piku nie jest > 0.

**Zasięg:** ZERO wpływu na trajektorie, energię, detektor,
klasyfikacje i triage (kontrola vac i tak miała kandydat=False —
przejścia=0; żaden warunek progowy nie używał tej wartości).
Pierwotny output pośredni zachowany j.w.

## (b) Addendum deskryptywne (nie-korekta; zero zmian detektora)

LOCK §3 (handoff) nakazuje OBOWIĄZKOWO deskryptywną tabelę
„ω_peak vs przewidywane ω(a)" (miękka, bez progu). Zamrożona
operacjonalizacja pomiaru detektorowego (MD §4: segment
[max(50,T−2000),T], wymagana długość ≥1000 j.cz.) czyni ω_peak
niemierzalnym dla biegów o podtrzymaniu <1000 j.cz. — czyli dla
WSZYSTKICH biegów gasnących Etapu A. Aby tabela obowiązkowa nie
była pusta, raport dostaje dodatkową kolumnę **ω_desc**: pik FFT
(to samo okno Hanna + interpolacja paraboliczna) na CAŁYM oknie
podtrzymania [t₁,t₂] (długość ≥200 j.cz.), z raportowanym
Δω=2π/L — oznaczaną gwiazdką i podpisem „deskryptywnie, segment
krótszy niż wymóg detektora".

**ω_desc NIE wchodzi do żadnego warunku detektora (cond3 nadal
wyłącznie z pomiaru wg MD §4), nie zmienia kategorii ani werdyktu;
mapa ω(a) pozostaje miękka.**
