---
title: "Phase_FINAL_close — zamknięcie: bramka A (audyt maszynerii 2) FAIL na A3 — STOP całego cyklu przed Phase 1"
date: 2026-08-23
type: phase-final-close
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: CLOSED-GATE-STOP
verdict: "A1 PASS (g0_crit=8/5 z |Δ|=3.7e−11; formuła 2(α+2)/[2(α+2)−d] potwierdzona na 4 parach (α,d)≠(2,3) z |Δ|≤1.1e−10); A2 PASS (ω_tail=1.00000 dla α∈{1,2,3}, |ω−1|<1e−4, oscylacje w oknie ważności, kontrola negatywna procedury czysta); A3 FAIL MERYTORYCZNY (δ_e=−75.50° vs −81.4±2°, δ_μ=+88.48° vs +38.6±2°, Δ(e→μ)=163.98° vs 120±1°; niereprodukowalne w ŻADNYM z 6 wariantów układu, 7 oknach fitu, 2 konwencjach fazy, R-kontroli; podtest amplitudy A≈|g0−1| PASS 0.993–1.011); A5 ROZSTRZYGNIĘTE (różnica = znak członu źródłowego W→−W; dwa różne układy, oba wewnętrznie spójne); A6 wykonane. Reguła bramki LOCKa: A3 FAIL → STOP CAŁEGO CYKLU. Phase 1–3 NIEURUCHOMIONE (zero obliczeń poza bramką)."
anti_lakatos_lock: PRESERVED
tags: [machinery2-audit, ODE-form-A, oscillating-tail, phase-map, gate-stop, negative-result, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[PhaseA_report.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/AUDYT_KRYTYCZNY_2026-07-28.md]]"
---

# Phase FINAL — zamknięcie cyklu (STOP na bramce A)

## 0. Werdykt

**Bramka A ZAMKNIĘTA: A3 FAIL merytoryczny.** Zgodnie z regułą bramki
(LOCK §3) i drzewem decyzyjnym (LOCK §6): *„A1–A3, A5: FAIL
któregokolwiek → STOP; raport «maszyneria 2 pada»"* — **cały cykl
zatrzymany PRZED Phase 1.** Rachunek centralny (runaway w kąpieli
sąsiadów, ω²(n)) **NIE ZOSTAŁ WYKONANY** — pytanie binarne cyklu
pozostaje otwarte, ale nie może być liczone na maszynerii, której
warstwa fazowa jest niereprodukowalna (fazy per gatunek to dokładnie
wejście Phase 1 → Phase 2 tego cyklu).

Wynik negatywny bramki jest pełnoprawnym rozstrzygnięciem: maszyneria 2
pada na audycie tak, jak maszyneria 1 padła 2026-08-10 — ale w INNYM
miejscu i w węższym zakresie (patrz §2: co przeżyło audyt).

## 1. Werdykty per test (liczby)

| Test | Werdykt | Kluczowe liczby |
|---|---|---|
| **A1** próg kolapsu | **PASS** | (α=2,d=3): num=1.6000000000, \|Δ\|=3.66e−11 (gate ≤1e−6); pary (1,3): 2.0, (3,3): 10/7, (2,4): 2.0, (1.5,2): 1.4 — wszystkie \|Δ\|≤1.1e−10; niezależny stałokrokowy RK4 w f=g^(α+1), batch-refinement (nie bisekcja rdzenia) |
| **A2** ogon ω=1 | **PASS** | ω_fit=1.00000 dla α∈{1,2,3} (g₀=0.9, okno [40,150], R²≥0.99998, 35 zer (g−1) w oknie, \|g−1\|≤0.002≤0.2); kontrola negatywna: syntetyczny e^(−r) — 0 zer, R²=0.003 → procedura poprawnie odrzuca (osiągalny FAIL potwierdzony) |
| **A3** fazy i amplitudy | **FAIL** | δ_e=−75.50° (target −81.4±2°, odchył 5.90°); δ_μ=+88.48° (target +38.6±2°, odchył 49.88°); Δ(e→μ)=163.98° (target 120±1°, odchył 43.98°); **podtest amplitudy PASS**: A/\|g₀−1\| = 0.993/0.997/1.004/1.011 dla g₀∈{0.95,0.98,1.02,1.05} (gate 1±10%); **kontrola R PASS**: dA/A=0.0, dδ=0.0000° (R₀=200 vs 300); okno ≤0.75·R ✓ |
| **A4** audyt 8 skryptów | **WYKONANE** (deskryptywne) | tabela per skrypt w [[PhaseA_report.md]]; outputy w `PhaseA_A4_output_*.txt` |
| **A5** niespójność istnienia | **ROZSTRZYGNIĘTE** | różnica zlokalizowana: ZNAK członu źródłowego (maszyneria 2: +g²(1−g) ⟺ W̃=g⁷/7−g⁸/8, W̃″(1)=−1, ogon oscylacyjny; AUDYT Dod. A: u²(u−1) ⟺ W=u⁸/8−u⁷/7, W″(1)=+1, runaway profilu). Numerycznie: (M2) g₀=1.2491/1.5 ograniczone+oscylacje, g₀=2.0212/3.1891 kolaps (próg 1.6 ✓); (AUD) wszystkie 4 g₀ runaway ✓. Obiekt ogona dla ew. Phase 1 = układ (M2) |
| **A6** korekta p134e−g | **WYKONANE** | powtórzone jawnie: A_τ monotonicznie rosnąca; r₃₁=3477.5 = kres skanu (g₀≤4), nie fizyczne maksimum; Δ(μ→τ) nie osiąga ±120°; Q_K=3/2 = INPUT |

## 2. Diagnoza A3 — dlaczego FAIL jest merytoryczny, nie implementacyjny

Trzy iteracje diagnostyczne (`PhaseA_A3_diag*.py` + outputy), wykonane
i udokumentowane PRZED werdyktem (wzorzec T3 z RETRO):

1. **Konwencja fazy wykluczona:** Δ(e→μ) jest niezależna od konwencji
   (offsety znoszą się w różnicy); testowane atan2(C,B) i atan2(−C,B);
   żadna nie zbliża obu faz jednocześnie (najlepsza: odchyły 5.9°/49.9°).
2. **Okno fitu wykluczone:** 7 okien od [20,35] (konwencja a3d) po
   [150,300] — dryf δ ≤ 2.8°, A stabilne do 3 cyfr. R-kontrola czysta.
3. **Wariant układu wykluczony:** testowane (i) plain α=2 kanoniczna,
   (ii) α=1 (hipoteza „siodła" mapy przy g₀≈1.99 = próg kolapsu α=1),
   (iii) Form B (α=1, źródło 1−g), (iv) F-S log bez odbić,
   (v) biegnące α_eff(g)=2/(1+η_K(g−1)²), η_K=181/15, g₀_e=0.90548
   (podstawienie bezpośrednie), (vi) wariant EL z K=g^(2α_eff(g)).
   δ_μ we wszystkich wariantach: +68…+93° — nigdy +38.6° ani +43.8°.
4. **Walidacja mojego integratora o WŁASNE skrypty rdzenia** (nie
   o .tex): A_tail(0.95)=0.049650 (mój) vs 0.0496530 (atail_functional,
   zgodność 5 cyfr); A_tail(1.4550)=0.5939 (mój) vs ≈0.594 (interpolacja
   siatki rdzenia: A(1.44474)=0.5726, A(1.47105)=0.6291). **Skrypty
   rdzenia same przeczą wartości dodatekH A_μ=0.3861 przy g₀=1.4550.**
5. **τ przy g₀=4 (biegnące α_eff): KOLAPS** w niezależnej integracji —
   spójne z tym, że lim(α→0) progu 2(α+2)/[2(α+2)−3] = 4, tj. g₀_τ=4
   leży NA granicy istnienia (dodatkowa rysa maszynerii, raportowana).

**Wniosek:** wartości dodatekH lin. 1126–1129 (δ_e=−81.4°, δ_μ=+38.6°,
Δ=120.01°≈2π/3) — w tym flagowy wynik „Δ(e→μ) potwierdzone 2π/3" —
**nie są reprodukowalne z udokumentowanego ODE maszynerii 2 w żadnym
zbadanym wariancie**, a warstwa amplitudowa tych samych tabel jest
niezgodna z bieżącymi skryptami rdzenia. Pochodzenie tych liczb jest
niezidentyfikowane (→ NEEDS N1).

## 3. Co przeżyło audyt (uczciwy bilans — nie wszystko pada)

- **Próg kolapsu g₀_crit = 2(α+2)/[2(α+2)−d] = 8/5**: potwierdzony
  niezależnie do 1e−10 (A1) — solidny.
- **ω_tail = 1, uniwersalne względem α; ogon oscylacyjny** (nie zanik
  wykładniczy): potwierdzone niezależnie (A2) — solidne. Mechanizm
  drabiny odległości z RETRO (locku oscylacyjnego) opiera się na ω=1
  i NIE pada wraz z fazami.
- **A_g ≈ |g₀−1| przy próżni** (1±1.1%): potwierdzone (A3-podtest).
- **Q_K=3/2 przy g₀_τ=1.5696 (1.9% pod progiem)**: skrypt rdzenia
  reprodukuje (r₃₁=3477.44, +0.006% vs PDG) — NIE audytowane niezależnie
  w tym cyklu (poza zakresem bramki), odnotowane deskryptywnie w A4.
- **Pada:** warstwa FAZOWA ogona per gatunek (δ_e, δ_μ, Δ=2π/3) —
  dokładnie ta, którą Phase 1 tego cyklu miała ekstrahować jako wejście
  do rachunku centralnego Phase 2.

## 4. Czego NIE policzono (jawnie)

- Phase 1 (κ, φ, A per gatunek; d* dla par) — **NIEURUCHOMIONA**.
- Phase 2 (baseline #63 V3, skan n, spektra, ewolucje, ω²_min(n)) —
  **NIEURUCHOMIONA**. Pytanie binarne cyklu (czy mod runaway dostaje
  ω²>0 przy skończonym n) — **NIEROZSTRZYGNIĘTE**.
- Phase 3 — nie dotyczy (wymaga V-PASS).
- Zero obliczeń poza bramką A i jej diagnostyką (zgodnie z LOCKiem).

## 5. Zgodność z LOCKiem

- Kryteria, progi i okna — bez zmian po starcie obliczeń.
- Diagnostyka A3 (diag 1–3) to poszukiwanie błędu IMPLEMENTACJI
  dopuszczone przez LOCK §4 pkt 1; udokumentowana w skryptach i tutaj
  PRZED użyciem werdyktu; błędu implementacji NIE znaleziono → FAIL
  utrzymany (żadna korekta kryteriów nie nastąpiła).
- Kontrole negatywne wykonane i czyste (A2-kontrola). P1c/P2c —
  nie dotyczy (fazy nieosiągnięte).
- Rdzeń `.tex` nietknięty; STATE.md nietknięty; brak commitów.
- Q_K=3/2 flagowane jako INPUT wszędzie, gdzie się pojawia (A6).

## 6. Następne kroki

User-gated — patrz [[NEEDS.md]]. Ten cykl pozostaje CLOSED-GATE-STOP;
ewentualna realizacja Phase 1–2 na fazach ZMIERZONYCH (zamiast
dokumentowanych) wymagałaby NOWEGO cyklu z nowym lockiem.

## 7. Pliki cyklu

- `PhaseA_audit_machinery2.py` + `PhaseA_output.txt` — bramka A1–A3, A5, A6.
- `PhaseA_A3_diag.py/.txt`, `PhaseA_A3_diag2.py/.txt`,
  `PhaseA_A3_diag3.py/.txt` — diagnoza A3 (wykluczenie błędu implementacji).
- `PhaseA_A4_runner.py`, `PhaseA_A4_run_a3d.py` (przechwycenie bajtowe —
  a3d wymaga UTF-8, capture tekstowy cp1250 pada), `PhaseA_A4_runner_log.txt`,
  `PhaseA_A4_output_*.txt` (8 szt.) — faktyczne outputy skryptów rdzenia.
- `PhaseA_report.md` — pełny raport z tabelą A4.
- `Phase_FINAL_close.md` (ten plik), `NEEDS.md`, `README.md` (status).
