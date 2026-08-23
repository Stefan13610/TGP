---
title: "NEEDS — proponowane edycje core po op-spectral-analysis-Phi (user-gated)"
date: 2026-07-03
type: needs
tgp_owner: research/op-spectral-analysis-Phi-2026-07-03
status: "N1–N5 EXECUTED 2026-07-03 (user-gate przyznany: 'dopisz NEEDS, zaktualizuj rdzeń'); N6 → cykl op-wall-dynamics-2026-07-03 (Phase 0 rozpisane, realizacja: następna sesja, osobny agent)"
---

# NEEDS — edycje core wynikające z CP-7 (wszystkie addytywne, user-gated)

## N1 (priorytet HIGH) — `sek08b` `ssec:spectral-synthesis-L03` / `thm:spectral-synthesis-L03`

Twierdzenie syntezy głosi σ(L̂) ⊂ [0,∞) „dla wszystkich fizycznych tł".
CP-7 (BVP numeryczny): prawdziwe TYLKO w formulacji grawitacyjnej F-A;
w formulacji solitonowej (funkcjonał, którego EL = ODE korony) próżnia
g=1 ma tachioniczne kontinuum od −γ, a solitony μ/τ mają 2/3 mody
zlokalizowane λ<0 (siodłowość). Proponowana edycja: zawęzić zakres
twierdzenia do F-A + dodać `rem:spectral-CP7` z wynikiem negatywnym
(cytując ten cykl) — analogicznie do wzorca CP-2 (Twierdzenie →
Aproksymacja/zakres).

## N2 (HIGH) — `sek08b` `cor:ghost-artifact`

Zawęzić zakres: „artefakt log-aproksymacji" jest poprawny dla sektora
grawitacyjnego; dla mechanizmu selekcji N=3 korony (ls10/#56) ściana
jest AKTYWNYM składnikiem dynamiki (odbicia μ/τ; min f(g)≈0,04), a
odbicie jest regularyzacją ad-hoc (nie-EL). Dodać remark o napięciu +
odnośnik do tego cyklu.

## N3 (MED) — `papers_external/paper_lepton_masses` Limitations (T-OP4 nowy)

Dodać punkt: „spectral stability of the mu/tau solitons is OPEN:
they are saddle points (index 2/3) of the static energy functional;
stability presumably requires the wall-reflection dynamics (non-
variational) or a constrained (fixed-charge) formulation — analysis:
op-spectral-analysis-Phi-2026-07-03". Elektron: 0 modów
zlokalizowanych (czysty).

## N4 (MED) — `audyt/L03_K_phi_stability/` POST_ACTION_UPDATE_2026-07-03

Dyspozycja rozdzielona: F-A = CLOSED-RESOLVED (numerycznie), F-S =
OPEN-RECLASSIFIED (zmierzony wynik negatywny, łączy się z L04).

## N5 (LOW) — `dodatekA_notacja.tex` l.404–412

Sync z dokładną formą Q (Phase 1): Q = W″ − ½F″(u₀′)² − F′[u₀″+(2/r)u₀′]
(dotychczasowy zarys pomija człony F′/F″ na tłach niejednorodnych).

## N6 (RESEARCH, nie-core) — interpretacja dynamiczna ściany
### [ROZSZERZONE 2026-07-03 o interpretację autora; → cykl `op-wall-dynamics-2026-07-03`]

**Interpretacja fizyczna autora (sesja #61, do sformalizowania jako
hipoteza robocza cyklu):** ściana NIE jest artefaktem — wynika z
**wewnętrznej energii solitonu**, tj. z **ilości tworzonej przestrzeni**,
która w rdzeniu obiektu przekracza próg stabilności. Konsekwencja
formalna: ilość tworzonej przestrzeni jest wielkością **budżetowaną**
⇒ fizycznie dopuszczalne fluktuacje muszą zachowywać budżet
(więz typu ∫δΦ d³x = 0 lub odpowiedni ładunek), a stabilność μ/τ
należy liczyć **na podprzestrzeni więzu** (analogia: Q-ball — funkcjonał
bez więzu ma kierunki ujemne, przy ustalonym ładunku konfiguracja jest
minimum).

Ścieżki cyklu: (a) stabilność z więzem budżetu przestrzeni (projekcja
operatora CP-7 na podprzestrzeń ∫δΦ=0 i warianty ważone); (b) ściana
jako warunek jednostronny dla fluktuacji (g ≥ g*, problem komplementarny);
(c) regularyzacje ściany (soft wall) i zależność spektrum μ/τ od
regularyzacji; (d) relacja ściany dolnej g*≈0,78 do górnego ogranicznika
rdzenia g_crit=8/5 (H7/H8) — dwa progi jednego mechanizmu budżetowego?;
(e) czy indeks siodłowy 0/2/3 koreluje z generacją (SPECULATIVE).
LOCK: [[../op-wall-dynamics-2026-07-03/Phase0_balance.md]].

---

## Log wykonania (2026-07-03, sesja #61)

- N1 ✅ `sek08b` — zakres `thm:spectral-synthesis-L03` zawężony do F-A + `rem:spectral-CP7`.
- N2 ✅ `sek08b` — `rem:ghost-artifact-scope-CP7` po `cor:ghost-artifact`.
- N3 ✅ lepton paper — Limitations T-OP4 (spectral stability μ/τ OPEN).
- N4 ✅ `audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-03.md`.
- N5 ✅ `dodatekA_notacja.tex` — sync pełnej formy Q.
- N6 → Phase 0 zalockowane w `op-wall-dynamics-2026-07-03` (realizacja: następna sesja).
