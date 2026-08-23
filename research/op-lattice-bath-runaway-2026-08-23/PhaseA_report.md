---
title: "PhaseA_report — audyt maszynerii 2 (bramka): A1 PASS, A2 PASS, A3 FAIL, A5 rozstrzygnięte, A4/A6 wykonane → STOP"
date: 2026-08-23
type: phase-report
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: GATE-FAIL-A3
anti_lakatos_lock: PRESERVED
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
---

# Phase A — raport audytu maszynerii 2 (BRAMKA)

**Werdykt bramki: A1 PASS · A2 PASS · A3 FAIL (merytoryczny) ·
A5 ROZSTRZYGNIĘTE · A4/A6 wykonane → BRAMKA ZAMKNIĘTA, STOP CAŁEGO
CYKLU przed Phase 1** (reguła bramki LOCK §3; drzewo §6).

Metoda: niezależna re-implementacja ODE — własny stałokrokowy klasyczny
RK4 (bez `scipy.integrate`), wektoryzowany po warunkach początkowych;
próg kolapsu przez batch-refinement przedziału; całkowanie Form A
w zmiennej kanonicznej f=g^(α+1) (kolaps regularny). Rdzeń używa
adaptacyjnego `solve_ivp` (RK45/DOP853) i bisekcji — metody rozłączne.
Pełny output: `PhaseA_output.txt`.

## 1. A1 — próg kolapsu (PASS)

Gate: |Δ|≤1e−6 dla (α,d)=(2,3); formuła 2(α+2)/[2(α+2)−d] na ≥3 parach
≠(2,3); FAIL przy >1e−3.

| α | d | formuła | numerycznie | \|Δ\| |
|---|---|---|---|---|
| 2.0 | 3 | 1.6000000000 | 1.6000000000 | 3.66e−11 |
| 1.0 | 3 | 2.0000000000 | 2.0000000000 | 3.66e−11 |
| 3.0 | 3 | 1.4285714286 | 1.4285714286 | 3.66e−11 |
| 2.0 | 4 | 2.0000000000 | 2.0000000001 | 1.10e−10 |
| 1.5 | 2 | 1.4000000000 | 1.4000000000 | 3.66e−11 |

Granica przeżycie/kolaps monotoniczna we wszystkich rundach (brak
anomalii). **A1 PASS.**

## 2. A2 — ogon oscylacyjny ω=1 (PASS)

Gate: ω_tail=1.000±0.01 dla ≥3 α; oscylacje w oknie |g−1|≤0.2.

- α=1: ω_fit=1.00000, R²=0.999992, 35 zer (g−1) w oknie [40,150], max|g−1|=0.002
- α=2: ω_fit=1.00000, R²=0.999987, 35 zer
- α=3: ω_fit=1.00000, R²=0.999981, 35 zer
- Kontrola negatywna procedury (syntetyczny zanik e^(−r), bez oscylacji):
  0 zer, R²(fit oscylacyjny)=0.003 → procedura poprawnie ODRZUCA;
  osiągalny FAIL potwierdzony. **A2 PASS.**

Uniwersalność ω=1 względem α potwierdzona niezależnie — ta warstwa
maszynerii 2 (fundament locku oscylacyjnego z retrospektywy 2026-08-16)
audyt PRZEŻYWA.

## 3. A3 — fazy i amplitudy (FAIL merytoryczny)

Gate (LOCK): δ_e=−81.4°±2°, δ_μ=+38.6°±2°, Δ(e→μ)=120°±1°
(dodatekH lin. 1126–1129; układ: biegnące α_eff(g)=2/(1+η_K(g−1)²),
η_K=181/15, g₀_e=0.90548, g₀_μ=φ·g₀_e); A_g≈|g₀−1| (1±10%);
okno ≤0.75·R; R-kontrola ≤1%.

**Zmierzone (okno [40,150], R₀=200, konwencja δ=atan2(C,B) jak
dodatekH lin. 1059):**

| wielkość | zmierzone | target | odchył | gate |
|---|---|---|---|---|
| δ_e | **−75.50°** | −81.4±2° | 5.90° | FAIL |
| δ_μ | **+88.48°** | +38.6±2° | 49.88° | FAIL |
| Δ(e→μ) | **163.98°** | 120±1° | 43.98° | FAIL |
| A/\|g₀−1\| (g₀=0.95/0.98/1.02/1.05) | 0.993/0.997/1.004/1.011 | 1±10% | — | PASS |
| R-kontrola (200 vs 300) | dA/A=0.0; dδ=0.0000° | ≤1% | — | PASS |
| δ_τ(g₀=4) | +43.87° (fit przed KOLAPSEM) | −27.3° (info) | — | poza gate |

**Diagnoza (3 iteracje, `PhaseA_A3_diag*.py`, wykonane i
udokumentowane PRZED werdyktem — wzorzec T3 z RETRO; błędu
implementacji NIE znaleziono, korekta kryteriów NIE nastąpiła):**

1. **Konwencja fazy — wykluczona.** Δ(e→μ) jest niezależna od
   konwencji (offset znosi się w różnicy); testowane atan2(C,B)
   i atan2(−C,B): najlepsza daje odchyły 5.9°/49.9° — żadna nie mieści
   obu faz w ±2°. Zmierzone Δ≈164° vs wymagane 120±1°.
2. **Okno fitu — wykluczone.** 7 okien od [20,35] (konwencja a3d) do
   [150,300]: dryf δ_μ ≤2.8°, δ_e ≤1.3°, A stabilne do 3 cyfr.
   R-kontrola czysta (IVP — brak ściany pudła; odnotowane).
3. **Wariant układu — wykluczony.** Testowane: (i) plain α=2
   (kanoniczna Form A); (ii) α=1, źródło g²(1−g) — hipoteza z „siodła"
   mapy fazowej przy g₀≈1.99 = próg kolapsu α=1 (2.0); (iii) Form B
   (α=1, źródło 1−g); (iv) F-S log bez odbić; (v) biegnące α_eff
   (podstawienie bezpośrednie); (vi) biegnące α_eff jako EL
   z K=g^(2α_eff(g)) (człon α_eff′·ln g). δ_μ we wszystkich:
   +68…+93° — nigdy +38.6° (target biegnący) ani +43.8° (mapa plain).
4. **Walidacja mojego integratora o WŁASNE skrypty rdzenia** (nie
   o .tex): A_tail(0.95)=0.049650 (mój RK4) vs 0.0496530
   (atail_functional, solve_ivp) — zgodność 5 cyfr; A_tail(0.90):
   0.0997 (mój, g₀=0.8993) vs 0.098982 (rdzeń, g₀=0.90); ω=1
   i A≈|g₀−1| zgodne z atail_asymptotic. **Implementacja wykluczona
   jako źródło rozbieżności.**
5. **Sprzeczność dodatekH z własnymi skryptami rdzenia:** mapa fazowa
   p127–p128 podaje A_μ=0.3861 przy g₀=1.4550; siatka
   atail_functional daje A(1.44474)=0.5726 i A(1.47105)=0.6291 ⟹
   A(1.4550)≈0.59 (mój RK4: 0.5939). Wartość 0.3861 nie pochodzi
   z tego równania. Analogicznie „siodło" mapy przy g₀≈1.99 leży NAD
   progiem kolapsu α=2 (1.6) — nie istnieje w kanonicznej Form A.
6. **τ przy g₀=4 (biegnące α_eff): KOLAPS** w niezależnej integracji
   (h=0.001 i 0.0005; też g₀=3.999). Spójne z lim(α→0) progu
   2(α+2)/[2(α+2)−3]=4 — g₀_τ=4 leży na granicy istnienia (→ NEEDS N4).
7. **Żaden z 8 audytowanych skryptów rdzenia nie liczy faz δ_e/δ_μ
   z dodatekH** — wartości −81.4°/+38.6°/Δ=120.01° nie mają
   zidentyfikowanego skryptu-źródła w audytowanym zestawie (→ NEEDS N1).

**Wniosek A3: FAIL merytoryczny.** Warstwa fazowa ogona per gatunek —
w tym flagowe „Δ(e→μ)=120.01°≈2π/3 potwierdzone" — jest
niereprodukowalna z udokumentowanego ODE maszynerii 2. To dokładnie ta
warstwa, którą Phase 1→2 tego cyklu miały brać jako wejście (fazy par
w drabinie d*), więc STOP jest merytorycznie konieczny, nie tylko
formalny.

## 4. A4 — audyt 8 skryptów rdzenia (deskryptywny, obowiązek kompletności)

Wszystkie 8 skryptów uruchomiono i przechwycono faktyczne outputy
(`PhaseA_A4_output_*.txt`); wyniki skryptów NIE są traktowane jako
prawda (LOCK). Kolumny wg LOCK A4: (i) czy istnieje test zdolny dać
FAIL; (ii) zgodność SUMMARY z faktycznym outputem; (iii) rejestr
INPUT vs DERIVED.

| skrypt (czas) | (i) osiągalny FAIL? | (ii) SUMMARY vs output | (iii) INPUT → DERIVED |
|---|---|---|---|
| **ngen_collapse_proof_v47b** (rc=0, 1309 s) | TAK (próg err<1e−6 per przypadek) — i FAIL WYSTĄPIŁ: (α,d)=(0.5,3) err=4.63e−6 | **CZĘŚCIOWO NIEZGODNE:** SUMMARY „VERIFIED (all d) …to 1e−10" — własna tabela: błędy 5.1e−9…4.6e−6, w tym 1 FAIL; ostatnia linia „All tests: SOME FAIL" uczciwa, ale narracja 1e−10 (powielona w status_map O-L5) nie jest poparta tym przebiegiem. „First integral verification: max rel err 7.97e−01" bez flagi (średnia 6.2e−5) | INPUT: g₀_e=0.8677 (kalibracja r₂₁), φ-FP, g₀_τ=1.5696 (z Q_K=3/2 = INPUT). DERIVED: formuła progu (potwierdzona), marginesy 1.9%/58.7% |
| **gcrit_energy_proof_v47b** (**NIE UKOŃCZONO W OKNIE AUDYTU**: 2 niezależne przebiegi przerwane po ~1665 s CPU każdy bez SUMMARY; stdout blokowo buforowany → zero częściowego outputu; przyczyna: ~3000+ całkowań rtol=1e−12/max_step=0.01/r_max=500 w bisekcjach ×80 iteracji; brak trybu szybkiej weryfikacji — ustalenie audytowe, szczegóły w pliku outputu) | NIE — zero kryteriów PASS/FAIL (analiza statyczna; skrypt eksploracyjny: „czy całka tarcia ma zamkniętą formę?") | SUMMARY (statycznie) uczciwie otwarte: całka tarcia = funkcjonał trajektorii, NIE domyka dowodu g₀_crit dla d>1; zgodności z outputem nie dało się zweryfikować w oknie (odnotowane) | INPUT: formuła g₀_crit — **zakładana w procedurze bracketingu, czyli dowodzona rzecz jest wejściem własnej weryfikacji**. DERIVED (deklarowane): numeryczny bilans −W(g₀c)=całka tarcia |
| **gcrit_pohozaev_v47b** (rc=0, 640 s) | NIE — zero kryteriów; same „Approaches" | **NIEZGODNE:** SUMMARY pkt 3 „DERRICK/POHOZAEV: Gives T/V = d/(d−2)" — własny output: T/V=0.0192 vs 3.0 (d=3), 0.0053 vs 2.0 (d=4), 0.0017 vs 1.67 (d=5). Relacja wirialna NIE potwierdzona przebiegiem (całki zdominowane niecałkowalnym ogonem oscylacyjnym r²·A²), a SUMMARY podaje ją jako wynik | INPUT: formuła g₀_crit. DERIVED: algebraiczna tożsamość W(g₀c) (zweryfikowana), zapis (g₀c−1)/g₀c=d/N; dowód analityczny d>1 — OPEN (uczciwie) |
| **atail_asymptotic_v47b** (rc=0, 165 s) | NIE — brak kryteriów; narracja badawcza korygowana w locie („Wait…", „EUREKA") | **CZĘŚCIOWO:** wnioski końcowe (ω=1, A≈\|g₀−1\|, decay rate −0.024≈0 ⟹ oscylacje) zgodne z liczbami; ale sekcje 1–3 zawierają błędną linearyzację w g (twierdzą „h″+2h′/r−h=0, exponential" — znak sprzeczny z własnym kodem, który daje +h); γ≈0.234 podawana jako wykładnik, choć własny fit ma resid 0.09 (por. collapse_exponent: γ dryfuje do 0.14) | INPUT: g₀ leptonów (kalibracja+φ+Koide). DERIVED: ω=1, A≈\|g₀−1\| (0.989–1.011), γ fit; model interpolacyjny NIE odtwarza Q_K=3/2 (uczciwie) |
| **atail_functional_v47b** (rc=0, 104 s) | NIE — czysta tabela, zero kryteriów | ZGODNE (deskryptywne; Q_K=1.5003 z x=A² wydrukowane) | INPUT: g₀ leptonów. DERIVED: siatka A_tail(g₀) — **kluczowa dla A3: A(1.455)≈0.59 ≠ 0.3861 z dodatekH p127** |
| **ode_koide_formA_exact_v47b** (rc=0, 269 s) | CZĘŚCIOWO — jawna ścieżka negatywna („No Koide bracket") istnieje; brak formalnych PASS/FAIL | ZGODNE: Q_K=1.499999999993, r₃₁=3477.44 (+0.006% PDG), margines 1.90% — liczby = narracja | INPUT: r₂₁=206.768 (kalibruje g₀_e), φ-FP (g₀_μ), **Q_K=3/2 (definiuje g₀_τ — INPUT)**. DERIVED: r₃₁=3477.44 (jedyny nietrywialny output), g₀_τ=1.5696 (1.9% pod progiem) |
| **collapse_exponent_v47b** (rc=0, 458 s) | NIE — eksploracyjny | ZGODNE i uczciwe: SUMMARY („γ NIE jest czystą liczbą wymierną; dryfuje; brak zamkniętej formy") wprost z outputu (lokalna γ: 0.37→0.12 przy δ→0) | DERIVED: tabela γ(α,d); ustala, że γ≈0.20–0.23 z atail_* to artefakt zakresu δ |
| **a3d_soliton_brannen_r** (rc=0, 2 s; UWAGA: wymaga konsoli UTF-8 — przechwycenie tekstowe cp1250 pada UnicodeDecodeError, odnotowane) | TAK — 6 testów z progami; FAIL WYSTĄPIŁ (T5: brak μ* dającego r_B=√2; zmierzone μ=1.37 vs deklarowane 4.12; r_B(μ_num)=1.572) → 5/6 | **CZĘŚCIOWO NIEZGODNE:** docstring „Wynik oczekiwany: 6/6 PASS"; finalny box „POTWIERDZONE… NUMERYCZNIE ZAMKNIĘTY" pomija T5 FAIL; sekcja „POTWIERDZONE NUMERYCZNIE: A_tail~(g₀−g*)^μ z μ~4.1" sprzeczna z własnym pomiarem μ=1.37 w tym samym przebiegu | INPUT: g₀_e=1.24915 (gałąź **F-S log z odbiciami ad hoc** od ściany-ducha — obiekt z AUDYT Dod. A!), φ-FP, g₀_τ=3.18912 z Q_K=3/2 (INPUT). DERIVED: r₃₁≈3477.4; r_B=√2 wynika z Q_K=3/2 przez tożsamość algebraiczną (thm T3-sqrtN1) — nie jest niezależnym potwierdzeniem. Fazy tej gałęzi: −89.2°/+78.1° — też ≠ dodatekH |

**Wnioski przekrojowe A4:**
- 5 z 8 skryptów nie ma żadnego testu zdolnego dać FAIL (anty-wzorzec
  z audytu maszynerii 1); 2 mają i w obu FAIL faktycznie wystąpił
  w tym przebiegu (ngen: (0.5,3); a3d: T5) — bez odnotowania
  w narracjach SUMMARY.
- Dwie jawne niezgodności SUMMARY↔output: gcrit_pohozaev (relacja
  wirialna „dana", faktycznie sfalsyfikowana przez własne liczby)
  i ngen (precyzja „1e−10 dla wszystkich" vs 4.6e−6 + FAIL).
- Rejestr wejść: **Q_K=3/2 jest INPUT** w każdym miejscu, gdzie
  pojawia się g₀_τ (ngen, ode_koide, a3d, atail_*); r₂₁ i φ-FP —
  również INPUT (kalibracja). Jedyny nietrywialny DERIVED łańcucha
  Koide: r₃₁=3477.44 (+0.006% PDG) — przy czym per p134e−g NIE wolno
  cytować r₃₁ jako „wyprowadzonego" w wariancie biegnącym (kres skanu).
- Żaden audytowany skrypt nie generuje faz δ_e/δ_μ z dodatekH
  lin. 1126–1129 (→ NEEDS N1).

## 5. A5 — rozstrzygnięcie niespójności istnienia (ROZSTRZYGNIĘTE)

Równanie-po-równaniu:

- **(M2) maszyneria 2** (dodatekH/O-L5, wszystkie skrypty v47b):
  `g″ + 2g′/r + (α/g)g′² = +g²(1−g)`, α=2 — EL funkcjonału
  z F=g^(2α)=g⁴ i **W̃(g)=g⁷/7−g⁸/8** (W̃″(1)=−1 <0: próżnia
  statycznie „tachionowa" → ogon oscylacyjny ω=1; solitony istnieją
  dla g₀<8/5).
- **(AUD) AUDYT_KRYTYCZNY Dodatek A** („F-A kanoniczna"):
  `u″ + 2u′/r + 2u′²/u = u³−u² = −u²(1−u)` — EL z **W(u)=u⁸/8−u⁷/7**
  (W″(1)=+1 >0: próżnia stabilna → zanik wykładniczy; profil z g₀>1
  nie zawraca → runaway → brak solitonu).

**Różnica: wyłącznie ZNAK członu źródłowego (W→−W).** Zmienna, człon
kinetyczny (2u′²/u), wymiar i BC — identyczne. Weryfikacja numeryczna
(własny RK4): (M2) g₀=1.2491/1.5 → ograniczone+oscylacje, 2.0212/3.1891
→ kolaps (próg 1.6 ✓); (AUD) wszystkie 4 → runaway ✓ (zgodnie
z tabelą AUDYT-u). Oba układy wewnętrznie spójne; opisują RÓŻNE teorie.
**Obiekt, którego ogon miał iść do Phase 1 = układ (M2).**

Uwaga trzecia (raportowana, nie gate): EL akcji F-S (W_FS=g³/3−g⁴/4)
z K=g⁴ daje trzecie równanie `u″+2u′/r+2u′²/u=(1−u)/u²` — INNE od obu;
„kanoniczna Form A" rdzenia definiuje więc źródło na poziomie równania,
a nie wyprowadza go z akcji. Który znak W wynika z AKCJI TGP —
nierozstrzygnięte w rdzeniu (→ NEEDS N2).

## 6. A6 — korekta p134e−g (powtórzona jawnie)

1. A_τ(g₀^τ) jest **monotonicznie rosnąca** — nie istnieje fizyczne
   maksimum amplitudy.
2. **r₃₁=3477.5 to kres skanu (g₀^τ≤4), nie fizyczne maksimum.**
   Zakaz cytowania r₃₁ jako „wyprowadzonego".
3. Δ(μ→τ) **nie osiąga ±120°** dla żadnego g₀^τ.
4. **Q_K=3/2 = trzeci parametr WEJŚCIOWY TGP** (12 ścieżek wyprowadzenia
   zawiodło — status_map O-L5); wszystko, co od niego zależy
   (g₀_τ w A3/A4, łańcuch Koide), flagowane jako INPUT.

## 7. Werdykt i skutek

**A3 FAIL → STOP całego cyklu przed Phase 1** (Phase 1–3
NIEURUCHOMIONE; zero obliczeń poza bramką i jej diagnostyką).
Szczegóły zamknięcia: [[Phase_FINAL_close.md]]; dalsze kroki
(user-gated): [[NEEDS.md]].
