---
title: "Phase_FINAL_close — zamknięcie: Q-H1-PULL (źródło ρ̂ przy λ̃≥0.2 INDUKUJE ucieczkę z dziedziny — COLLAPSE zbieżnie przy λ̃=0.5 [t=0.375 obie siatki], centralny runaway w rdzeniu źródła; λ̃_crit∈(0.05,0.2]; przy λ̃≤0.05 DEFORMATION Yukawa-podobna z saturacją nieliniową — gate liniowy 10.7%>5% FAIL już przy λ̃=0.01) + Q-H2-INCONCLUSIVE (STABILIZED=0; 8/12 COLLAPSE; 1/12 RADIATED potwierdzone h/2+dt/2 [qR3−0.20@λ̃=0.05: kolaps baseline'u UCHYLONY, pole osiada w zdeformowanej próżni]; 3/12 INCONCLUSIVE-RUN — w tym 2 przetrwańcy t=1000 z E_ref≤0 i 1 rozjazd siatek na krawędzi stabilności)"
date: 2026-09-15
type: phase-final-close
tgp_owner: research/op-collapse-matter-source-2026-09-14
status: CLOSED
verdict: "Q-H1: PULL wg litery (LOCK §4: ≥1 λ̃ daje THRESHOLD-PULL lub COLLAPSE zbieżnie — λ̃=0.5: COLLAPSE zbieżnie h05/h025 z t_end 0.375/0.375; λ̃=0.2: COLLAPSE na h05 [BREAKDOWN t=0.600]; λ̃∈{0.01,0.05}: DEFORMATION [0.01 zbieżnie obie siatki]; λ̃_crit∈(0.05,0.2] opisowo). Gate liniowy (λ̃=0.01): FAIL 10.7%>5% — odpowiedź jest Yukawa-KSZTAŁTNA, ale spłycona nieliniowo (δψ(0) meas/lin = 0.893 przy λ̃=0.01, 0.691 przy λ̃=0.05), więc Q-H1-DEFORMATION nie byłoby orzeczone nawet bez PULL. Q-H2: INCONCLUSIVE wg litery (LOCK §4: PASS wymaga ≥1 STABILIZED potwierdzonego — jest 0; FAIL wymaga 12/12 COLLAPSE zbieżnie — jest 8 COLLAPSE [kontrola negatywu qR3+0.20@0.05 h/2 zgodna: t_end 1.415/1.415] + 1 RADIATED potwierdzone h/2 i dt/2 [qR3−0.20@0.05: τ=44.7/44.6/44.7, kolaps baseline t=4.74 NIE zachodzi; pole osiada w zdeformowanej próżni E_core→−2.837=E_static(0.05)] + 3 INCONCLUSIVE-RUN [g−0.30σ3@0.05 i @0.2: przeżywają t=1000 bez zdarzenia brzegowego, ale E_ref=E_core(50)≤0 wyłącza zamrożone progi energetyczne — osiadłe stany ψ(0)=0.866 i 0.611; qR3−0.20@0.2: rozjazd potwierdzeń h025=COLLAPSE vs h05/dt2=RADIATED]). Deskryptywnie: kolaps indukowany źródłem = centralny runaway w rdzeniu ρ̂ (r<1), przejście przez pas w JEDNYM kroku (|Δψ|≫1; podtyp sufit/podłoga = kierunek overshootu, nieinformatywny); przetrwańcy relaksują DO zdeformowanej próżni (E_core→E_static: −2.837@0.05, −27.49@0.2), nie zachowują energii startu — stąd zero STABILIZED wg litery. P1-H1 PASS (𝒰_mat=λ̃ρ̂ψ²/(4−3ψ) WYPROWADZONE, gate'y simplify=0); P1-H2 PASS (δψ=−5λ̃(G_Yuk∗ρ̂), znak δψ<0 POTWIERDZONY w Phase 3); P1-H3 PASS (𝒰_mat→+∞ przy ψ→4/3⁻, →0 przy ψ→0⁺). P2 PASS 3/3 po korekcie 1 (estymator sekularny; pierwotny FAIL 2.28e−6 zdiagnozowany jako dt²-owy offset ΔC2 stanu początkowego — skalowanie 4.00/4.00, plateau okien 2–10; po korekcie 8.8e−9/9.0e−9 ≪ 1e−6, przewidywanie noty trafione). INCONCLUSIVE ≠ pozytyw; zakazy claimów dotrzymane."
anti_lakatos_lock: PRESERVED
tags: [collapse-matter-source, matter-coupling, L-mat-unified, yukawa-deformation, threshold-pull, source-induced-collapse, breakdown-boundary, radiated, inconclusive-qh2, second-order-dynamics, healthy-branch, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_1_p2c_estimator.md]]"
  - "[[NEEDS.md]]"
  - "[[README.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# Phase FINAL — zamknięcie cyklu op-collapse-matter-source

**Status: CLOSED-EXECUTED (2026-09-15, jedna sesja: MD FROZEN →
Phase 1 (sympy) → uzupełnienie 𝒰_mat w MD §2 wynikiem P1-H1 →
Phase 2 (1 korekta estymatora, progi/konfiguracja nietknięte) →
Phase 3 Q-H1 (6 biegów) → Phase 3 Q-H2 (12 głównych + 8 potwierdzeń
+ 1 kontrola negatywu) → zamknięcie).** Kryteria LOCKa stosowane
DOSŁOWNIE; zero zmian kryteriów/progów/detektorów/startów/listy λ̃/
ρ̂/sponge po pierwszym biegu produkcyjnym.

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q-H1** (binarne M911-N2: deformacja czy indukcja) | **Q-H1-PULL** (litera LOCK §4) | λ̃=0.5: COLLAPSE zbieżnie (t_end=0.375 na h05 i h025); λ̃=0.2: COLLAPSE (h05); λ̃≤0.05: DEFORMATION — źródło POTRAFI wciągnąć pole za progi; λ̃_crit∈(0.05,0.2] |
| **Q-H2** (stabilizacja kolapsu — centralne NEEDS-N2) | **Q-H2-INCONCLUSIVE** (litera LOCK §4) | STABILIZED=0 (nie PASS); COLLAPSE=8/12 (nie FAIL); RADIATED=1 potwierdzone (kolaps uchylony!), INCONCLUSIVE-RUN=3 |
| P1-H1 (wyprowadzenie 𝒰_mat) | **PASS** | 𝒰_mat=λ̃ρ̂ψ²/(4−3ψ) z literalnego √−g·(q/Φ₀)ψρ; ∂𝒰_mat/∂ψ=λ̃ρ̂ψ(8−3ψ)/(4−3ψ)²; tożsamość próżniowa λ̃ρ̂(ψ−1)(ψ+4)/(4−3ψ); simplify=0, num ≤3.6e−14 |
| P1-H2 (odpowiedź zlinearyzowana, PRE-REJESTRACJA) | **PASS** | (−∇²+1)δψ=−5λ̃ρ̂ ⟹ δψ=−5λ̃(G_Yuk∗ρ̂); kwadratura I₀=0.776061935; **znak δψ<0 POTWIERDZONY** w Phase 3 |
| P1-H3 (fakty brzegowe, PRE-REJESTRACJA) | **PASS** | lim 𝒰_mat = +∞ (ψ→4/3⁻), = 0 (ψ→0⁺); przewidziana asymetria — konfrontacja: §6 |
| P2 (bramka) | **PASS 3/3** (po korekcie 1) | próżnia 0.0 (≤1e−10); regresja qR3−0.20: BB t=4.7400 (4.74±2%); dryf sekularny 8.8e−9/9.0e−9 (≤1e−6/100T₀) |
| Gate liniowy Q-H1 (λ̃=0.01) | **FAIL** (10.7% > 5%) | Saturacja nieliniowa mierzalna już przy λ̃=0.01 — WYNIK deskryptywny, nie błąd (LOCK: „nieliniowe odchylenie = wynik") |

**Uwagi zalockowane, stosowane dosłownie:** INCONCLUSIVE ≠ pozytyw;
ZAKAZ claimów o masach leptonów/oscylonach dotrzymany (cykl-bliźniak);
S05 dotrzymane (ρ statyczne, zero pól dynamicznych).

## 1. Wejścia i formy (rejestr [INPUT]; MD §1–§3, §9)

- Sektor pola CYTAT (poprzednik `op-r3-stationary-states`, Q-D1-PASS,
  |g^tt|): M=ψ⁶/(4−3ψ)², 𝒦=ψ⁴, 𝒰=ψ⁴/4−ψ³/3; K_geo=γ=c₀=1; dziedzina
  (0,4/3), pas 4/3−1e−6 / 1e−6, ZERO podłóg/barier.
- Człon materii Z RDZENIA (tylko odczyt): eq:L-mat-unified
  L_mat=−(q/Φ₀)ψρ + L01 ρ≡−T^μ_μ/c₀² + eq:vol-element-M911
  √−g=c₀ψ/(4−3ψ) ⟹ (P1-H1) 𝒰_mat=λ̃ρ̂ψ²/(4−3ψ), λ̃=(qc₀/Φ₀)ρ₀>0;
  ρ̂=exp(−r²/18) FROZEN. EOM: Mψ̈+½M′ψ̇²=(1/r²)(r²𝒦ψ′)′−½𝒦′ψ′²−𝒰′−∂𝒰_mat/∂ψ.
- Silnik: kopia engine_core poprzednika + człon materii w F (punktowo,
  wariacyjnie zgodny z energią: d[λ̃ρ̂(ψ−1)(ψ+4)/(4−3ψ)]/dψ =
  λ̃ρ̂ψ(8−3ψ)/(4−3ψ)² dokładnie) i w ewaluatorze energii (tożsamość
  sfaktoryzowana bez kancelacji — dziedzictwo correction note 2
  poprzednika). Dla λ̃=0 bitowo równoważny (P2b: t_end 4.7400 = 4.74
  baseline'u co do kroku).
- Störmer–Verlet uogólniony na (ψ,π), dt=0.005 (potwierdzenia
  dt=0.0025 na h=0.05); h∈{0.05,0.025}, R=200, sponge smootherstep
  γ₀=1 na [160,200] (ON produkcyjnie, OFF w P2c); E_core r≤80
  gęstość PEŁNA (z 𝒰_mat−𝒰_mat(1), próżniowo odjęta), E_ref=E_core(50).
- Starty Q-H2 IDENTYCZNE jak u poprzednika (definicje MD §6 tamtego
  cyklu): qR3 a=±0.20 (sinc(r/π)·exp(−r²/450)), gauss a=−0.30 σ=3,
  gauss a=+0.15 σ=6; π₀=0; źródło od t=0. Baseline λ̃=0: 4× COLLAPSE
  (BREAKDOWN-BOUNDARY górny; t=4.74 / 2.275 / 6.64 / 18.185 na h05).

## 2. Phase 1 (sympy; `Phase1_output.txt`)

- P1-H1: 𝒰_mat literalnie = −√−g·L_mat = λ̃ρ̂ψ²/(4−3ψ) — gate
  simplify=0; pochodna i tożsamość próżniowa simplify=0; kontrola
  numeryczna ψ∈{0.5,1,1.2}: max |Δ|=3.55e−14 (próg 1e−12·max(1,|v|)).
- P1-H2: O(ε) EOM statycznego = f″+2f′/r−f−5λ̃ρ̂ (gate simplify=0);
  redukcja sferyczna kernela Yukawy (jakobian i granice zweryfikowane
  sympy): δψ(r)=−(5λ̃/r)∫r′ρ̂(r′)[e^{−|r−r′|}−e^{−(r+r′)}]/2 dr′;
  δψ(0)=−5λ̃I₀, I₀=0.776061935; kontrola rezydualna ODE kwadratury:
  5.7e−6 względnie. PREDYKCJA: δψ<0 wszędzie (True na siatce), ogon
  e^{−r}/r.
- P1-H3: granice sympy: +∞ / 0; siła materii znika na podłodze
  (lim ∂𝒰_mat/∂ψ = 0 przy ψ→0⁺).

## 3. Phase 2 — bramka (PASS 3/3 po korekcie 1; `Phase2_output.txt`)

| Gate | Wynik | Próg |
|---|---|---|
| P2a próżnia (λ̃=0, sponge ON, 100T₀, obie siatki) | ‖ψ−1‖∞ = 0.0 / 0.0 | ≤1e−10 |
| P2b regresja (λ̃=0, qR3−0.20, h05) | BREAKDOWN-BOUNDARY t_end=4.7400 | 4.74±2% |
| P2c energia ze źródłem (λ̃=0.05, a=+0.05 σ3, sponge OFF, t=700) | dryf sekularny 8.801e−9 (h05) / 8.993e−9 (h025) | ≤1e−6/100T₀ |

**Korekta 1** (`Phase_correction_note_1_p2c_estimator.md`, PRZED
użyciem): pierwotny estymator okienny [0,10T₀]vs[90T₀,100T₀] dał
2.28e−6/2.30e−6 (FAIL) — diagnostyka (`Phase2_diag_output.txt`)
wykazała czyste skalowanie dt² (ilorazy 4.00/4.00) i plateau okien
2–10 (znaki „++-+---++"): mierzony był dt²-owy offset hamiltonianu-
cienia stanu początkowego (ΔC2, znany z diagnostyki poprzednika),
nie zmiana sekularna. Korekta: okna plateau [10T₀,20T₀]vs[90T₀,100T₀]
(×100/80 na 100T₀); próg i konfiguracja BEZ ZMIAN; przewidywanie noty
(~1e−8) trafione. Pierwotny output zachowany
(`Phase2_output_pre_correction1.txt`). Deskryptywnie: offset ΔC2
2.28e−6, fit LSQ [10,100]T₀ 1.6e−7, max|E−E0|=3.4e−4 (dt²).

## 4. Phase 3 — Q-H1 (źródło na próżni; `Phase3_qh1_output.txt`)

| λ̃ | siatka | klasa | δψ(0) zmierzone | δψ(0) liniowe (−5λ̃I₀) | meas/lin |
|---|---|---|---|---|---|
| 0.01 | h05 | DEFORMATION | −0.034647 | −0.038803 | 0.893 |
| 0.01 | h025 | DEFORMATION | −0.034653 | −0.038803 | 0.893 |
| 0.05 | h05 | DEFORMATION | −0.134018 | −0.194015 | 0.691 |
| 0.2 | h05 | COLLAPSE (BREAKDOWN t=0.600) | n/a | −0.776062 | — |
| 0.5 | h05 | COLLAPSE (BB górny t=0.375) | n/a | −1.940155 | — |
| 0.5 | h025 | COLLAPSE (BB górny t=0.375) | n/a | −1.940155 | — |

- Klasy końcowe: 0.01 DEFORMATION (zbieżnie), 0.05 DEFORMATION,
  0.2 COLLAPSE, 0.5 COLLAPSE (zbieżnie, t_end identyczne) ⟹
  **Q-H1-PULL**; λ̃_crit ∈ (0.05, 0.2] (opisowo).
- Osiadłość deformacji: V=8.8e−5 (0.01) / 1.15e−4 (0.05) ≪ 0.01·D;
  min ψ̄ = 0.9654 (0.01) / 0.86598 (0.05) — oba NAD progiem 5/6
  (0.8333); zero THRESHOLD-PULL w zbadanej liście.
- Gate liniowy (λ̃=0.01, r≤40): err=10.7% > 5% FAIL — kształt
  Yukawa-podobny, ale spłycony; saturacja rośnie z λ̃ (0.893→0.691).
  Kierunek zgodny z P1-H3: siła materii słabnie przy ψ<1
  (∂𝒰_mat/∂ψ ∝ ψ(8−3ψ)/(4−3ψ)²).
- Deskryptywnie (probe lokusa, `Phase3_desc_locus_output.txt`):
  kolaps indukowany = centralny runaway w rdzeniu źródła (r<1,
  ρ̂≈1); pas przekraczany w JEDNYM kroku z ogromnym overshootem
  (λ̃=0.5: max ψ=3339 @ r=0.175, min ψ=−0.29 obok; ψ(0)=10.5) —
  podtyp „górny/dolny" jest kierunkiem overshootu, nie fizyczną
  preferencją granicy.

## 5. Phase 3 — Q-H2 (stabilizacja; `Phase3_qh2_output.txt`, `Phase3_results/verdict.json`)

Kategorie końcowe par (baseline λ̃=0: COLLAPSE, wszystkie 4 starty):

| start | λ̃=0.05 | λ̃=0.2 | λ̃=0.5 | baseline t_end (λ̃=0) |
|---|---|---|---|---|
| qR3 a=−0.20 | **RADIATED** (potw. h/2+dt/2; τ=44.7/44.6/44.7) | INCONCLUSIVE-RUN (rozjazd: h05/dt2 RADIATED τ=30.5 vs h025 COLLAPSE t=0.570) | COLLAPSE t=0.230 | 4.74 |
| qR3 a=+0.20 | COLLAPSE t=1.415 (kontrola negatywu h/2: 1.415 — zgodna) | COLLAPSE t=0.790 | COLLAPSE t=0.495 | 2.275 |
| gauss a=−0.30 σ3 | INCONCLUSIVE-RUN (przeżywa t=1000; E_ref=−1.37≤0) | INCONCLUSIVE-RUN (przeżywa t=1000; E_ref≤0) | COLLAPSE t=0.160 | 6.64 |
| gauss a=+0.15 σ6 | COLLAPSE t=1.600 (LOWER) | COLLAPSE t=0.825 (BREAKDOWN) | COLLAPSE t=0.535 | 18.185 |

KLASY: STABILIZED=0, RADIATED=1, COLLAPSE=8, INCONCLUSIVE-RUN=3 ⟹
**Q-H2-INCONCLUSIVE** (litera: PASS wymaga ≥1 STABILIZED
potwierdzonego; FAIL wymaga 12/12 COLLAPSE zbieżnie).

Deskryptywnie (obowiązkowe konteksty; `Phase3_desc_survivors_output.txt`):
- **Kolaps UCHYLONY (3 pary bez zdarzenia brzegowego do t=1000):**
  qR3−0.20@0.05 (RADIATED potwierdzone) oraz g−0.30σ3@0.05 i @0.2
  (INCONCLUSIVE-RUN wyłącznie przez E_ref≤0). Wszyscy przetrwańcy
  relaksują DO zdeformowanej próżni: ψ(0,1000)=0.86598@0.05
  (= deformacja Q-H1 co do 1e−5; osiadły statycznie V≤0.01D) i
  0.6105–0.6077@0.2 (stan PODPROGOWY ψ<5/6, osiadły dla g−0.30σ3);
  E_core→E_static: −2.837@0.05, −27.49@0.2.
- **Zero STABILIZED wg litery:** kategoria wymaga zachowania energii
  startu (E_core(t_max)≥0.05·E_ref>0); przetrwańcy oddają energię
  startu do promieniowania i osiadają w studni materii (E_core<0) —
  litera klasyfikuje to jako RADIATED (E_ref>0) albo INCONCLUSIVE-RUN
  (E_ref≤0), nie jako stabilizację.
- **Krawędź stabilności:** qR3−0.20@0.2 — kategoria niezbieżna
  między siatkami (h05 przeżywa, h025 kolabuje w t=0.570);
  uczciwie INCONCLUSIVE-RUN.
- Wszystkie 8 COLLAPSE: t_end=0.16–1.6 ≪ baseline'y (4.74–18.185) —
  przy silnym źródle kolaps jest SZYBSZY niż bez materii (indukcja,
  spójnie z Q-H1-PULL).

## 6. Konfrontacja z P1-H2/P1-H3 (obowiązkowa wg handoffu)

- **Znak δψ<0 (P1-H2): POTWIERDZONY.** Deformacje ujemne przy
  0.01/0.05; głębokość spłycona względem liniowej (0.893/0.691) —
  saturacja nieliniowa zgodna kierunkowo z formą siły (∂𝒰_mat/∂ψ
  maleje przy ψ<1). Ogon Yukawy: kształt zgodny, gate ilościowy 5%
  nie przechodzi (10.7%) już przy λ̃=0.01.
- **Asymetria sufit/podłoga (P1-H3): BEZ rozstrzygnięcia w formie
  przewidzianej.** Przewidywano: stabilizacja łatwiejsza dla kolapsów
  górnych. Obserwacja: (i) jedyny w pełni potwierdzony przypadek
  uchylenia kolapsu (qR3−0.20@0.05) dotyczy baseline'u SUFITOWEGO —
  kierunkowo zgodne; (ii) ale drugi sufitowy start (qR3+0.20) kolabuje
  przy każdym λ̃, szybciej niż baseline; (iii) mechanizm rzeczywisty
  kolapsu z materią to centralny runaway w rdzeniu ρ̂ z przejściem
  przez pas w jednym kroku — kierunek granicy (sufit/podłoga) jest
  artefaktem overshootu, nie sygnaturą asymetrii 𝒰_mat. Asymetria
  potencjału jest realna analitycznie (P1-H3 PASS), lecz dynamika
  silnego sprzężenia jest zdominowana przez indukcję, nie przez
  barierę sufitową.
- **Niespodzianka programowa (poza pre-rejestracją, deskryptywnie):**
  materia w formie korpusowej przy λ̃≥0.2 nie stabilizuje — sama
  DESTABILIZUJE próżnię (Q-H1-PULL). Odczyt (a) NEEDS-N2 („kolaps =
  artefakt braku materii") uzyskuje TYLKO częściowe, nie-literowe
  wsparcie: 3/12 par unika kolapsu baseline'owego przez relaksację do
  zdeformowanej próżni, ale żadna nie spełnia litery STABILIZED.

## 7. Korekty i incydenty (pełna lista; anti-Lakatos PRESERVED)

1. **Correction note 1** (P2c, estymator dryfu sekularnego):
   udokumentowana PRZED użyciem wyniku; pierwotny output zachowany;
   diagnoza dt² (4.00/4.00) + plateau; próg/konfiguracja nietknięte;
   przewidywanie noty (~1e−8) trafione (8.8e−9/9.0e−9). Uwaga
   kosmetyczna: nagłówkowa linia REJESTR w `Phase2_output.txt`
   nadal wymienia stare okna — wiążący jest opis w sekcji P2c
   i correction note.
2. Incydent harnessu (bez wpływu na wyniki): TypeError w księgowaniu
   rekordu COLLAPSE w `Phase3_qh1_response.py` (duplikat klucza
   t_end) — crash PO deterministycznej klasyfikacji biegu λ̃=0.2;
   poprawka księgowania, biegi powtórzone deterministycznie
   (cache json; wyniki ukończonych biegów niezmienione).
3. Incydent harnessu (bez wpływu): znak nie-ASCII w komentarzu
   `Phase3_qh2_stabilize.py` (SyntaxError przed startem obliczeń).
4. Zero zmian form/progów/detektorów/startów/listy λ̃/ρ̂/sponge po
   pierwszym biegu produkcyjnym; MD §2 uzupełnione wynikiem P1-H1
   PRZED pierwszym biegiem numerycznym (procedura z LOCKa/handoffu).

## 8. Higiena

- Ścieżki relatywne od rootu vaulta; ZERO `cd`; artefakt zagnieżdżenia
  `TGP/TGP_v1/TGP/...` — sprawdzany po zapisach, NIE wystąpił.
- Rdzeń `.tex`, STATE.md, git — NIETKNIĘTE (sek08a tylko odczyt).
- Katalogi innych cykli — tylko odczyt.
- `integrity_snapshot.txt`: hashe LOCK+MD+silnik po FROZEN oraz stan
  po korekcie 1 i zamknięciu (zweryfikowane przy zamknięciu).
- Batch w tle + checkpointy npz co 100 j.cz. (`Phase3_results/
  checkpoint_*.npz`); aktywne czekanie między fazami (`wait_for.py`).

## 9. Pliki cyklu

`Phase0_balance.md` (LOCK) · `Phase_method_decisions.md` (FROZEN;
§2 = wynik P1-H1) · `engine_core.py` (kopia+cytat+materia) ·
`Phase1_matter_analytic.py` + `Phase1_output.txt` · `Phase2_gate.py`
+ `Phase2_output.txt` (+ `Phase2_output_pre_correction1.txt`,
`Phase2_diag_energy.py`, `Phase2_diag_output.txt`,
`Phase_correction_note_1_p2c_estimator.md`) · `Phase3_qh1_response.py`
+ `Phase3_qh1_output.txt` · `Phase3_qh2_stabilize.py` +
`Phase3_qh2_output.txt` · `Phase3_desc_locus.py` +
`Phase3_desc_locus_output.txt` · `Phase3_desc_survivors.py` +
`Phase3_desc_survivors_output.txt` · `Phase3_results/` (json+npz
per bieg, `qh1_summary.json`, `verdict.json`, checkpointy, logi) ·
`integrity_snapshot.txt` · `wait_for.py` · `NEEDS.md` · `README.md`.
