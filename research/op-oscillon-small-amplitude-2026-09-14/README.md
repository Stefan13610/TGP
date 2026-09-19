---
title: "README — op-oscillon-small-amplitude: właściwy test predykcji P1b (ω₂=−139/24<0) — oscylony małej amplitudy w gałęzi zdrowej?"
date: 2026-09-14
type: cycle-readme
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
folder_status: closed
status: CLOSED
verdict: "Q-G: INCONCLUSIVE wg litery (0×OSCILLON, 0×OSCILLON-WEAK, 9×RADIATED [τ=174–699; kontrola negatywu a=0.05 σ=6: τ=418.2 identyczne na h i h/2], 2×COLLAPSE zbieżnie dt/2 [σ=10, a≥0.08, BOUNDARY-UPPER t≈40–55], 1×INCONCLUSIVE-RUN [a=0.05 σ=10: podtrzymanie 880>628 j.cz., ale E_end=7.8%·E_ref w szczelinie kategorii, bez plateau]). Deskryptywnie: wszystkie mierzalne częstości rdzenia 1.0013–1.0016 = próg kontinuum, zero śladu zmiękczenia LP (mapa 0.9977→0.9421) i zero monotonii w a — pole małej amplitudy nie samopułapkuje się w klasie zbadanej (gauss a≤0.10, σ≤10, t≤10⁴). P1a′ PASS 24/24; P2 PASS 6/6 (regresja τ=206.9 dokładnie). INCONCLUSIVE ≠ pozytyw; falsyfikacja P1b NIE orzeczona (litera FAIL niespełniona)."
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
---

# op-oscillon-small-amplitude (2026-09-14)

Realizacja NEEDS **N1** poprzednika (`op-r3-stationary-states`,
Q-E-INCONCLUSIVE): konfrontacja pre-rejestrowanej predykcji
**P1b: ω₂=−139/24<0** (miękka nieliniowość ⟹ oscylony małej amplitudy)
w jej naturalnej domenie — **a≤0.10, t_max=10⁴, R=400** (poprzednik
sondował tylko |a|≥0.15, t≤10³).

- **Q-G:** czy istnieją długożyciowe oscylony małej amplitudy
  (E_core≥0.5·E_ref przez ≥100 T₀ + ≥50 przejść + **ω_peak≤0.99**
  — rdzeń związany, poniżej progu kontinuum; zbieżnie h/2, dt/2)?
- **Q-G-FAIL = falsyfikacja P1b w klasie zbadanej** (pełnoprawny wynik):
  wtedy brak nośnika oscylonowego w OBU klasach amplitud → user-gate
  o statusie łańcucha leptonowego.

12 startów gauss (a∈{0.02,0.05,0.08,0.10}×σ∈{3,6,10}) + próżnia;
staging FROZEN (Etap A t=2000 → triage → Etap B t=10⁴). Zakres:
istnienie; ZAKAZ claimów o masach i dyskretności (osobne cykle).
Cykl-bliźniak: [[../op-collapse-matter-source-2026-09-14/README.md]].

Status: **CLOSED — Q-G-INCONCLUSIVE** (werdykt w frontmatterze;
pełne zamknięcie: [[Phase_FINAL_close.md]], otwarte pozycje:
[[NEEDS.md]]).

## Log
- 2026-09-14 — LOCK zapisany (sesja główna; autoryzacja: wybór usera
  „N1: cykl małych amplitud" po zamknięciu poprzednika
  Q-E-INCONCLUSIVE). Obliczeń zero.
- 2026-09-15 — realizacja w jednej sesji (agent-implementator):
  MD FROZEN (silnik = kopia poprzednika, zmiana tylko parametryczna
  sponge [320,400], R=400, t_max=10⁴; tabela ω(a) zapisana PRZED
  Phase 3) → Phase 1 PASS (P1a′ 24/24; P1b′ cytat ω₂=−139/24) →
  Phase 2 PASS 6/6 (próżnia 0.0; regresja τ=206.9 odchyłka 0.000;
  dryf 9.65e−8; odbicie 1.13e−7) → Phase 3: Etap A 13 biegów h=0.05
  do t=2000 (równolegle), triage FROZEN: 12/12 startów martwych
  (E_core(2000)/E_ref=0.001–0.078 <0.2 lub kolaps), Etap B: tylko vac
  (do t=10⁴, ψ≡1 dokładnie); potwierdzenia: h/2 dla a=0.05 σ=6
  (τ=418.2 identyczne), dt/2 dla obu COLLAPSE (Δt≤0.06 j.cz.) →
  **WERDYKT Q-G-INCONCLUSIVE wg litery** (9×RADIATED, 2×COLLAPSE,
  1×INCONCLUSIVE-RUN, 0×OSCILLON/WEAK). Korekta 1 (tylko warstwa
  raportu: pozorny pik FFT sygnału zerowego + addendum ω_desc) —
  nota PRZED użyciem, output pośredni zachowany. Integralność:
  LOCK/MD/silnik UNCHANGED (SHA256). Deskryptywnie: ω_desc≈1.0013–
  1.0016 na progu kontinuum, zero śladu mapy LP — negatyw merytorycznie
  mocny; status łańcucha leptonowego → NEEDS N2 (user-gate).
