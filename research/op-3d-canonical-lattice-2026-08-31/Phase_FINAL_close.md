---
title: "Phase_FINAL_close — zamknięcie: Q-FAIL (PRIMARY po tłach istniejących) — sieć sc solitonów kanonicznych μ istnieje wyłącznie przy d=2π i jest tam TACHIONOWA: ω²_min(2π)=−1.674350 (zbieżnie, po odjęciu translacji), ucieczka nieliniowa t_esc≈4.06–4.20 ≤ 2·t*_izo=9.42 przy obu znakach; most P1c-kan PASS (kotwica kanoniczna −1.646589 reprodukowana w 3D z odchyłem 0.33% przy h=0.3); trend vs 1D: 3D NIE podnosi ω²_min (−1.674 < −1.222) — hipoteza stabilizacji gęstością zamknięta także w 3D sc (model kanoniczny)"
date: 2026-08-31
type: phase-final-close
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: CLOSED
verdict: "GATE Phase 1: PASS (wszystkie składowe; P1a/P1d dziedziczone od poprzednika z cytatem: maxerr 6.022e−4/3.441e−4, ratio 3.999; |ΔE|/E ≤7.42e−15). P1b-kan PASS: kotwica kanoniczna λ_min(w1)=−1.646589 (h=0.0125; zbieżność |λ(0.025)−λ(0.0125)|=1.12e−4 ≤ 1.65e−3; residuum profilu 2.89e−15 ≤ 1e−10); diagnostyka kieszeni: V(r)=4g−5g² ma minimum −4.872 W RDZENIU (r≈0), głębokość 3.87, FWHM=1.20 — kieszeń szeroka i rozdzielcza (FWHM/h_mostu=4.0), pułapka f_ε poprzednika NIE występuje; t*_ref=4.336 (dt-konsystencja 0.0000 ≤1%; kontrola K_ε=0.1: 4.260). P1c-kan PASS: λ_min(3D)=−1.656051 (N=76) / −1.651965 (N=100) vs kotwica — odchył 0.57%/0.33% (gate ±5% + trend poprawy); t*_izo(3D)=4.710 (biegi 4.720/4.710/4.900, wszystkie w oknie t*_ref±15% [3.686,4.986]). Phase 2: istnienie tła sieci sc WYŁĄCZNIE d=2π (oba starty NIESTAŁE-ZBIEŻNE: ‖R‖∞≤6.7e−10, ‖g−1‖∞=0.469≥0.05, dgrid 4.58e−3/4.59e−3 ≤5e−3); d∈{π, 3.0790}: UCIECZKA-g→0; d=3π: NIEZBIEŻNE-SIATKOWO (A1.0; N=48 sam daje kandydata ‖R‖=6.7e−11)/KOLAPS (A0.7); d=4π: KOLAPS-DO-PRÓŻNI. Phase 3 (rachunek centralny): ω²_min(2π)=−1.674350 (N=48; ZBIEŻNE |Δ(32→48)|=1.47e−2 ≤ 8.37e−2; argmin Γ; per k: Γ −1.6743/X −1.6639/M −1.6541/R −1.6448) po odjęciu 3 modów translacyjnych zidentyfikowanych PRZED interpretacją (overlap=1.0, coverage=1.0, λ_trans=−1.84e−3=O(h²), Rayleigh zgodny −1.83e−3); cross-check u-formy |Δ|≤1.85e−3. P3b: PASS 10/10 po korekcie 1 (d=4π N=48: eigsh k=10 nie mieści 8-krotnej degeneracji kończącej się na pozycji 10 — dowód przeciw FD-exact 1.3e−13, korekta k=32/ncv=240, maxerr 1.849e−3 ≤1e−2; correction note PRZED użyciem). P3c PASS: ω²_min=+1.06518/+1.06520 >0. Phase 4 (oba tła d=2π, pełna macierz): UCIECZKA (mechanizm dev>‖tło‖) t_esc = 4.060/4.060/4.120/4.200 (A0.7) i 4.120/4.120/4.060/3.980 (A1.0) — WSZYSTKIE ≤ 2·t*_izo=9.42 przy obu znakach; dt-konsystencja t_esc dokładna (4.060→4.060; 4.120→4.120); gate energii ≤2.7e−8 ≤1e−6. WERDYKT: Q-FAIL wg rulingu PRIMARY (Phase3_verdict_ruling.md, zapisany przed Phase 4): ω²_min<0 zbieżnie dla wszystkich istniejących teł ORAZ ucieczka ≤2·t*_izo. Strict obok: dla d∈{π, d*₁, 3π, 4π} sieć sc kanoniczna NIE ISTNIEJE w zbadanej klasie (negatyw istnienia, nie stabilności). Deskryptywnie (obowiązkowy punkt LOCKa): ω²_min(3D)=−1.674 vs 1D −1.222…−1.230 — trzeci wymiar NIE podnosi ω²_min, wręcz pogłębia niestabilność; hipoteza stabilizacji gęstością zamknięta także w 3D sc w modelu kanonicznym; decyzja aksjomatyczna o znaku W pozostaje bez podpory rachunkowej."
anti_lakatos_lock: PRESERVED
tags: [3d-lattice, canonical-action, tachyonic-sector, bloch-3d, negative-result, q-fail, density-stabilization-closed, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_verdict_ruling.md]]"
  - "[[Phase_correction_note_p3b_multiplicity.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-3d-canonical-lattice

Obliczenia i zamknięcie 2026-08-31 (jedna sesja implementatora). Kryteria
LOCKa (Phase0_balance.md §3, §6) stosowane DOSŁOWNIE; zero zmian
progów/punktów/siatek po starcie. Jedyna korekta = udokumentowany błąd
implementacji solvera (correction note PRZED użyciem wyniku, §5).

**Rejestr WEJŚĆ (flagowane):** g₀_e=0.90548 [r₂₁/φ-FP], φ (kalibracja μ:
g₀_μ=φ·g₀_e=1.4650974 — formuła zamrożona w MD; wartość „1.46507" w LOCKu
to zaokrąglenie transkrypcji, odchył 3e−5), d*₁=3.0790 [użyte w skanie d —
punkt bez istnienia], ε=0.2 (kontrola 0.1) [ewolucje K_ε], β=γ=1.

---

## 0. Werdykt

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| Gate Phase 1 (P1b-kan/P1c-kan) | **PASS** | most radialny→kartezjański w modelu kanonicznym reprodukuje kotwicę −1.646589 z odchyłem 0.33% przy h=0.3 i t*_izo w oknie ±15% — bariera rozdzielczości poprzednika była własnością f_ε, nie 3D |
| Phase 2 (istnienie) | **CZĘŚCIOWE** | sieć sc istnieje wyłącznie przy d=2π; π/d*₁ → ucieczka g→0, 3π/4π → kolaps do próżni (3π/A1.0: niezbieżność siatkowa) |
| **Q** (∃d: ω²_min(d)>0?) | **Q-FAIL** (PRIMARY po tłach istniejących; ruling przed Phase 4) | jedyne istniejące tło (d=2π) jest tachionowe ω²_min=−1.674350 (zbieżnie) i ucieka nieliniowo w t_esc≈4.1 ≤ 2·t*_izo=9.42 przy obu znakach |

**Sens (deskryptywnie):** to pełnoprawny negatyw RACHUNKOWY w 3D — po raz
pierwszy w tej linii istnieje liczba ω²_min(3D). Argument nieprzenośności
negatywu 1D (Hill, ogony 1/r) został skonsumowany: trzeci wymiar nie
uratował hipotezy stabilizacji gęstością (ω² spadło z −1.22 do −1.67).

## 1. Phase 1 — kotwica kanoniczna i most (Phase1_output.txt)

- **Dziedziczone (LOCK §1, cytaty w output):** P1a PASS (maxerr najniższej
  gałęzi 6.022e−4 tach / 3.441e−4 stab przy N=32; ratio N32/N64 =
  3.999/3.999), P1d PASS (|ΔE|/E ≤ 7.42e−15, dt-konsystencja ≤1.29e−13)
  — ten sam model kanoniczny, ten sam integrator, to samo środowisko
  (CPython 3.14.2, numpy 2.4.3, scipy 1.17.1); nie liczone ponownie.
- **P1b-kan (kotwica):** profil radialny kanoniczny μ (g₀=1.4650974;
  EL: g″+(2/r)g′+(2/g)g′²=g²(1−g)); residuum (operacjonalizacja FROZEN)
  2.89e−15 ≤ 1e−10; profil g∈[0.7773, 1.4651], oscylacyjny ogon 1/r.

  | h | λ_min(w1, u-forma) | waga-K (kontrola) |
  |---|---|---|
  | 0.05 | −1.647151 | −1.646608 |
  | 0.025 | −1.646701 | −1.646565 |
  | 0.0125 | **−1.646589 (KOTWICA)** | −1.646555 |

  Gate wewnętrzny: 1.12e−4 ≤ 1.65e−3 PASS. Formy w1/waga-K zgodne do
  3.4e−5 przy h=0.0125.
- **Diagnostyka kieszeni Q(r) (obowiązkowa, anty-pułapka po f_ε):**
  V(r)=4g−5g²: V_vac=−1; min = −4.872 przy r≈0 (rdzeń solitonu, nie
  powłoka!); głębokość D=3.872; **FWHM = 1.20** (r∈[0.006, 1.194]) ⟹
  FWHM/h_mostu(0.3) = 4.0 — kieszeń szeroka, gładka, w pełni rozdzielcza
  na siatkach mostu. Kontrast z poprzednikiem: f_ε miało powłokę −15.5
  o skali ~0.02–0.1 przy r=3.38. Forma ważona Q/K: min −3.509 przy r≈0.
- **t*_ref (ewolucja radialna, K_ε=0.2):** t*=4.336 (a=+0.01, dt=0.004
  i 0.002 — identyczne, dt-konsystencja 0.0000 ≤ 1%), 4.392 (a=−0.01);
  t*_ref := 4.336; kontrola K_ε=0.1: 4.260 (−1.8%). Detektor FROZEN:
  floor 0.01 / ceil 2g₀ / non-finite.
- **P1c-kan (MOST — bramka główna): PASS.**

  | siatka | h | λ_min(3D, waga-K) | odchył od kotwicy | u-forma 3D (deskr.) |
  |---|---|---|---|---|
  | N=76³ | 0.3947 | −1.656051 | 0.57% | −1.669432 |
  | N=100³ | 0.3000 | **−1.651965** | **0.33% (gate ≤5%)** | −1.659655 |

  Trend poprawy 0.0095 → 0.0054 PASS. t*_izo(3D, N=100): biegi
  (a,dt) = (+,0.02): 4.720; (+,0.01): 4.710; (−,0.02): 4.900 — wszystkie
  w oknie t*_ref±15% [3.686, 4.986] PASS; **t*_izo(3D) := 4.710**
  (Phase 4: t_end = 14.130).

## 2. Phase 2 — istnienie sieci sc (Phase2_output.txt, Phase2_backgrounds3d.npz)

| d | start A1.0 | start A0.7 | istnienie |
|---|---|---|---|
| π = 3.1416 | UCIECZKA-g→0 | UCIECZKA-g→0 | NIE |
| 2π = 6.2832 | **NIESTAŁE-ZBIEŻNE** | **NIESTAŁE-ZBIEŻNE** | **TAK** |
| d*₁ = 3.0790 [INPUT] | UCIECZKA-g→0 | UCIECZKA-g→0 | NIE |
| 3π = 9.4248 | NIEZBIEŻNE-SIATKOWO¹ | KOLAPS-DO-PRÓŻNI | NIE |
| 4π = 12.5664 | KOLAPS-DO-PRÓŻNI | KOLAPS-DO-PRÓŻNI | NIE |

Tło d=2π (N=48): g∈[0.6100, 1.4688], ‖g−1‖∞=0.4688 ≥0.05,
‖R‖∞ ≤ 6.7e−10 ≤ 1e−8, dgrid(16³) = 4.58e−3/4.59e−3 ≤ 5e−3.
Dedup: oba starty dają tła różniące się >1e−8 numerycznie (zapisane
osobno jako 2pi/A1.0 i 2pi/A0.7), ale spektralnie i dynamicznie
identyczne (¶4) — efektywnie JEDNO tło fizyczne.
¹ d=3π/A1.0: N=48 samodzielnie daje kandydata (‖R‖=6.7e−11, amp 0.586),
N=32 NO-CONV ⟹ dgrid=0.486 — punkt niezbieżny siatkowo per LOCK
(odnotowany w NEEDS jako niepewność jednostronna).

## 3. Phase 3 — RACHUNEK CENTRALNY (Phase3_output.txt, Phase3_results3d.json)

**Identyfikacja translacji (nieusuwalna, PRZED interpretacją):** w Γ
dokładnie 3 mody z overlap=1.000 (indeksy 1–3), λ_trans = −4.24e−3 (N=32)
/ −1.84e−3 (N=48) — skaluje jak O(h²) ✓; coverage c_i=1.0 (bez eskalacji);
iloraz Rayleigha ∂_ig niezależnie: −4.21e−3/−1.83e−3 — zgodny z ARPACK.

**Tabela ω²_min(d, k) po odjęciu translacji (N=48; oba tła identycznie):**

| k | ω²_min(2π) N=48 | N=32 |
|---|---|---|
| Γ | **−1.674350 (argmin)** | −1.689007 |
| X | −1.663878 | −1.678994 |
| M | −1.654096 | −1.669529 |
| R | −1.644804 | −1.660493 |

Zbieżność: |Δ(32→48)| = 1.47e−2 ≤ 0.05·max(1.674, 0.1) = 8.37e−2 —
**ZBIEŻNE**; wartość główna ω²_min(2π) = **−1.674350**. Pełne minimum
(z translacjami) identyczne (−1.674350 — mod tachionowy leży poniżej
translacji). Cross-check u-formy (N=32): Γ |Δ|=1.76e−3, R |Δ|=1.85e−3 ✓.

**P3b (próżnia, nieusuwalna):** 9/10 punktów PASS wprost (2.5e−3–9.5e−3);
d=4π N=48 FAIL 2.482e−1 → zdiagnozowany artefakt ARPACK (krotność 8
kończąca się na pozycji k=10; operator dokładny do 1.3e−13 vs FD-exact),
korekta 1 (k=32, ncv=240; PRZED użyciem): maxerr = 1.849e−3 ≤ 1e−2 —
**P3b po korekcie: PASS 10/10** (Phase_correction_note_p3b_multiplicity.md).
**P3c (sektor stabilny, nieusuwalna): PASS** — ω²_min = +1.06518 (N=32) /
+1.06520 (N=48) > 0 (‖R‖ ≤ 3.3e−10; źródło q=1, kontynuacja frakcjami).

## 4. Phase 4 — nieliniowa weryfikacja (Phase4_output.txt + Phase4_run0–3.txt)

Ruling kwantyfikatora zapisany PRZED Phase 4 (Phase3_verdict_ruling.md).
Superkomórka 2×2×2 (96³), mod minimalny NIE-translacyjny z Γ (λ=−1.674350,
zgodność z Phase 3 |Δ|=0), a = 0.01·‖tło‖ = 0.00469, t_end = 3·t*_izo
= 14.130; ucieczka wg detektora FROZEN (dominujący mechanizm: odchył
przekracza 100% tła).

| tło | bieg | t_esc | vs 2t*=9.42 | gate |ΔE|/E |
|---|---|---|---|---|
| 2π/A0.7 | +a ε=0.2 dt=0.01 | 4.060 | ≤ ✓ | 8.9e−11 |
| 2π/A0.7 | +a ε=0.2 dt=0.005 | 4.060 (dt-kons. dokładna) | ≤ ✓ | 5.6e−12 |
| 2π/A0.7 | −a ε=0.2 dt=0.01 | 4.120 | ≤ ✓ | 2.7e−08 |
| 2π/A0.7 | +a ε=0.1 dt=0.01 | 4.200 | ≤ ✓ | 1.0e−10 |
| 2π/A1.0 | +a ε=0.2 dt=0.01 | 4.120 | ≤ ✓ | 2.7e−08 |
| 2π/A1.0 | +a ε=0.2 dt=0.005 | 4.120 (dt-kons. dokładna) | ≤ ✓ | 1.7e−09 |
| 2π/A1.0 | −a ε=0.2 dt=0.01 | 4.060 | ≤ ✓ | 8.9e−11 |
| 2π/A1.0 | +a ε=0.1 dt=0.01 | 3.980 | ≤ ✓ | 2.7e−08 |

Wszystkie 8 biegów: UCIECZKA ≤ 2·t*_izo przy obu znakach; gate energii
≤ 2.7e−8 ≤ 1e−6 ✓. Pary {4.060, 4.120} między tłami zamienione znakami —
zgodne ze znakiem modu ARPACK (konwencja swobodna, macierz ± pokrywa oba)
i z identycznością fizyczną obu teł (diagnostyka dedup).

**Q-FAIL (litera LOCKa §3 Phase 4, PRIMARY):** ω²_min < 0 dla wszystkich
istniejących teł (zbieżnie) ✓ ORAZ ucieczka ≤ 2·t*_izo ✓.
Strict obok (ruling §2): dla d∈{π, d*₁, 3π, 4π} negatyw ISTNIENIA
(nie stabilności) w klasie zbadanej.

## 5. Obowiązkowy punkt LOCKa: trend ω²_min(3D) vs 1D (−1.22)

**ω²_min(3D, d=2π) = −1.674350** vs 1D (bloch-chain, d=2π…):
−1.222…−1.230. Trzeci wymiar **NIE podnosi** ω²_min — pogłębia
niestabilność o ~37%. Kierunek przewidywany argumentem nieprzenośności
(Hill nie obowiązuje w 3D — mogło pójść w obie strony) rozstrzygnięty
rachunkiem: w dół. Hipoteza „stabilność = własność konfiguracji
o skończonej gęstości" jest w klasie sc/kanonicznej zamknięta w 1D i 3D.

## 6. Korekty, incydenty, INCOMPLETE (pełna lista)

1. **Doprecyzowanie transkrypcyjne przed startem (MD, FROZEN):**
   g₀_μ zamrożone FORMUŁĄ φ·0.90548 = 1.4650974 (LOCK ma zaokrąglenie
   „1.46507"; źródło bath-two-sectors: 1.46510; odchył 3e−5 bez wpływu
   na bramki ±5%/±15%).
2. **Korekta 1 (Phase_correction_note_p3b_multiplicity.md, PRZED
   użyciem):** P3b d=4π N=48 — eigsh(k=10) konstrukcyjnie nie mieści
   8-krotnej degeneracji na pozycjach 3–10; dowód przeciw analitycznemu
   FD-symbolowi (|eigsh−FD|=2.5e−1 przy k=10; 1.3e−13 przy k=32);
   korekta parametru solvera (k=32, ncv=240) wyłącznie dla tego punktu
   kontroli. Pierwotny Phase3_output.txt zachowany bez zmian.
   Rodzina błędu = korekta 1 poprzednika (tol=0 przyjęte od startu
   w całym cyklu).
3. **Incydent środowiskowy (bez wpływu na liczby):** proces Phase 4
   ubity przez limit czasu tła (~55 min) PO ukończeniu 4/4 biegów tła
   A0.7 i przed biegami A1.0; kontynuacja `Phase4_continue_A10.py`
   (delta WYŁĄCZNIE implementacyjna: jeden bieg na proces; fizyka,
   parametry i detektory verbatim), 4 biegi A1.0 w osobnych procesach.
   Sandbox: blokada `/dev/null` i ścieżek poza vaultem — wszystkie
   rachunki przez pliki skryptów, pełne ścieżki, zero `cd`.
4. **Dedup teł:** progu 1e−8 nie osiągnięto (tła A1.0/A0.7 różnią się
   ~1e−4 na N=48), więc oba tła liczone w Phase 3 i Phase 4 w całości —
   wyniki identyczne do 6 miejsc (spektra) i lustrzane w znakach (t_esc):
   potwierdzenie, że to jedno tło fizyczne.
5. **INCOMPLETE: BRAK** — wszystkie zalockowane punkty policzone.

## 7. Pliki cyklu

Obliczeniowe: `Phase1_canonical_gate.py` + `Phase1_output.txt`;
`Phase2_relax_lattice3d.py` + `Phase2_output.txt` +
`Phase2_backgrounds3d.npz`; `Phase3_bloch3d.py` + `Phase3_output.txt` +
`Phase3_results3d.json`; `Phase3_addendum_p3b_diag.py` +
`Phase3_output_addendum_p3b.txt`; `Phase4_nonlinear3d.py` +
`Phase4_output.txt`; `Phase4_continue_A10.py` + `Phase4_run0.txt` …
`Phase4_run3.txt`. Metodologiczne: `Phase_method_decisions.md`,
`Phase3_verdict_ruling.md`, `Phase_correction_note_p3b_multiplicity.md`.
Zamykające: ten plik, `NEEDS.md` (user-gated), `README.md`
(zaktualizowany). Rdzeń `.tex`, STATE.md, git, katalogi innych cykli —
NIETKNIĘTE (odczyt wzorców).

## 8. Mapowanie na drzewo decyzyjne LOCKa §6

Gałąź: **„Q-FAIL → NEEDS: hipoteza stabilizacji gęstością zamknięta
także w 3D sc (model kanoniczny); decyzja aksjomatyczna o znaku W bez
podpory"** — wykonana dosłownie. Do NEEDS: konsekwencje ontologiczne
(user-gate), niepewność jednostronna d=3π, opcje bcc/fcc.
