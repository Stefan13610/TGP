---
title: "Phase0_balance — LOCK: re-lock cyklu 3D w modelu KANONICZNYM (opcja b z NEEDS poprzednika) — most radialny→kartezjański bez ściany f_ε, potem sieć sc + Bloch"
date: 2026-08-31
type: phase0-lock
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-31: wybór opcji 'b' z NEEDS op-3d-lattice-bath-stability (»Most w modelu KANONICZNYM (K=g⁴)«) — JAWNA DECYZJA USER-GATE: kotwice #63 (λ_min=−1.3896, t*=3.62) PRZESTAJĄ być bramką 3D; nowa kotwica kanoniczna mierzona wewnątrz cyklu. Autoryzacja obejmuje realizację ('rozpiszę LOCK i odpalę' → 'b')."
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/NEEDS.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
---

# Phase 0 — LOCK cyklu `op-3d-canonical-lattice`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria,
punkty skanu, kontrole, progi i forbidden moves zamknięte PRZED kodem.**

---

## 0. Pytanie binarne (Q) — identyczne co u poprzednika, w modelu kanonicznym

**Czy istnieje separacja d sieci prostej kubicznej solitonów kanonicznych μ
w sektorze tachionowym, przy której ω²_min(d) = min_k λ_min(L̂_d(k)) > 0
(po odjęciu modów translacyjnych)?**

Zmiana względem poprzednika (GATE-FAIL-STOP na P1c) — wyłącznie MODEL
mostu i kotwic: wszystko w akcji kanonicznej (K=g⁴, U=g⁷/7−g⁸/8),
spektra bez regularyzacji. Uzasadnienie (close poprzednika §2): kieszeń
−15.5 na powłoce g≈0.75 jest własnością REGULARYZACJI f_ε modelu #63
(f_log=0 przy g=e^{−1/4}), nie fizyki sektora; operator kanoniczny jest
wielomianowy — bez wąskiej powłoki. Pre-rejestracja nieprzenośności
negatywu 1D (Hill, ogony 1/r) dziedziczona bez zmian.

## 1. Dziedziczenie i reużycie (jawne)

- **Skrypty poprzednika gotowe i NIEURUCHOMIONE** (zero wyników):
  `../op-3d-lattice-bath-stability-2026-08-31/Phase2_relax_lattice3d.py`,
  `Phase3_bloch3d.py` — napisane od początku w akcji kanonicznej
  (method_decisions poprzednika); wolno skopiować do tego katalogu
  i użyć (z korektą eigsh tol=0 + pokrycie translacji z
  `Phase_correction_note_eigsh.md`).
- **P1d poprzednika (gate energii ewolucji kanonicznej K_ε=0.2/0.1)
  DZIEDZICZONY bez powtarzania:** PASS |ΔE|/E ≤ 7.42e−15, dt-konsystencja
  1.29e−13 (`../op-3d-lattice-bath-stability-2026-08-31/Phase1_output.txt`)
  — ten sam model, ten sam integrator.
- **P1a poprzednika (dyspersja próżni 3D kanoniczna) DZIEDZICZONY:**
  PASS maxerr 6.0e−4/3.4e−4, rząd 2 (tamże).
- Kalibracja gatunku: g₀_μ = φ·g₀_e = 1.46507 (g₀_e=0.90548) [INPUT,
  M-P z bath-two-sectors Phase 1]; poniżej progu 8/5 — soliton regularny.
- d*₁(μμ, M-P) = 3.0790 [INPUT — zmierzone w układzie KANONICZNYM M-P,
  spójne z tym cyklem].

## 2. Model ZAMKNIĘTY

- **M-3D-kan:** akcja kanoniczna; statyka radialna: g″+(2/r)g′+(2/g)g′²
  = g²(1−g); spektra weight-1 BEZ regularyzacji; ewolucje nieliniowe
  z K_ε (ε=0.2, kontrola 0.1 [INPUT]).
- **Obiekt PRIMARY:** soliton μ kanoniczny (g₀=1.46507 [INPUT]);
  deskryptywnie wolno dodać e (g₀=0.90548), NIE zastępować.
- **Skan d (zalockowany):** d ∈ {π, 2π, 3.0790 [INPUT], 3π, 4π};
  wolno dodać, NIE usuwać.
- **Siatki:** komórka sieci: N ∈ {32, 48}/wymiar (wolno dodać 64);
  most P1c: pudło L=30, h ∈ {0.4, 0.3} (N=76³, 100³).
  Spójność zalockowana TERAZ: dla d=4π i N=48 h=0.26 ≤ 0.3 —
  jeśli P1c przejdzie przy h=0.3, siatki sieci są adekwatne z zapasem.
- **Zbiór k:** Γ, X, M, R (wolno dodać, NIE usuwać).
- **Rejestr WEJŚĆ:** g₀_e=0.90548, φ (kalibracja μ), d*₁=3.0790,
  ε=0.2/0.1, β=γ=1 — flagować w każdym zależnym wyniku.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — nowa kotwica kanoniczna + most (bramka)
- **P1b-kan (kotwica radialna kanoniczna, tanio w 1D):** profil μ
  z radialnego EL (strzelanie/kolokacja; residuum ≤1e−10); λ_min(w1)
  operatora kanonicznego na h ∈ {0.05, 0.025, 0.0125}, R=60:
  **kotwica = wartość przy h=0.0125; gate wewnętrzny: |λ(0.025)−λ(0.0125)|
  ≤ 1e−3·|λ|**. OBOWIĄZKOWA diagnostyka pre-rejestrowana: profil
  pointwise Q(r) linearyzacji — min i szerokość FWHM ewentualnej kieszeni
  (żeby pułapka f_ε nie powtórzyła się po cichu). Dodatkowo t*_ref:
  radialna ewolucja izolowanego μ z K_ε (a=±0.01, dt×2) → t*_ref
  zmierzone (dt-konsystencja ≤1%).
- **P1c-kan (MOST, bramka główna):** profil interpolowany na 3D
  (L=30, padding do próżni): **λ_min(3D, h=0.3) = kotwica ±5%** ORAZ
  poprawa |λ(h=0.4)−kotwica| > |λ(h=0.3)−kotwica| (trend). t*_izo(3D,
  h=0.3) = t*_ref ± 15%. **FAIL ⟹ STOP** (raport bez Phase 2–4;
  NEEDS → opcja c poprzednika).
- P1a/P1d: DZIEDZICZONE (patrz §1) — odnotować w output, nie liczyć.

### Phase 2 — istnienie tła sieci sc (kanonicznie)
Jak u poprzednika (LOCK §3 Phase 2, przejęte dosłownie): relaksacja
‖residuum‖∞ ≤ 1e−8, 5 d × 2 starty (superpozycja profili PRIMARY
+ 0.7× amplitudy) × N ∈ {32,48}; niestałość ‖g−1‖∞ ≥ 0.05; zbieżność
siatkowa ≤ 5e−3. Brak istnienia WSZĘDZIE ⟹ CLOSED-GATE-STOP
(pełnoprawny negatyw istnienia sieci sc).

### Phase 3 — RACHUNEK CENTRALNY
Jak u poprzednika, przejęte dosłownie: λ_min(L̂_d(k)) w Γ/X/M/R,
N ∈ {32,48}; identyfikacja 3 modów translacyjnych (λ~O(h²), overlap
z ∂g/∂x_i) PRZED interpretacją — nieusuwalna; ω²_min po ich odjęciu;
zbieżność |Δω²_min| ≤ 0.05·max(|ω²_min|, 0.1); **P3b** próżnia
w komórce (gałęzie |k+G|²−γ ≤ 1e−2); **P3c** sektor stabilny m²=+γ
ze źródłem gaussowskim 3D (ω²_min>0 wymagane) — obie nieusuwalne.

### Phase 4 (warunkowa — przy zbieżnym Phase 3)
Superkomórka 2×2×2, perturbacja ±0.01·‖tło‖ wzdłuż modu minimalnego,
ewolucja K_ε do t = 3·t*_izo(3D); gate energii ≤ 1e−6.
- **Q-PASS:** ∃d: ω²_min(d)>0 (zbieżnie) ORAZ brak ucieczki do
  3·t*_izo przy obu znakach.
- **Q-FAIL:** ω²_min(d)<0 dla wszystkich ISTNIEJĄCYCH teł (zbieżnie)
  ORAZ ucieczka ≤ 2·t*_izo. (Ruling kwantyfikatora przy istnieniu
  częściowym: PRIMARY po tłach istniejących, strict obok — zapisać
  PRZED Phase 4.)
- **Q-INCONCLUSIVE:** pozostałe, wprost.
- Deskryptywnie obowiązkowo: trend ω²_min(3D) vs 1D (−1.22).

## 4. Forbidden moves
Jak u poprzednika (LOCK §4, przejęte dosłownie), plus:
8. Kotwica kanoniczna (P1b-kan) jest MIERZONA — zakaz jej strojenia;
   wartość i tabela zbieżności w output PRZED jakimkolwiek rachunkiem 3D.
9. Dziedziczenia P1a/P1d wolno użyć TYLKO bez zmian modelu/integratora;
   jakakolwiek zmiana ⟹ bramka liczona od nowa.

## 5. Deliverables
`Phase_method_decisions.md` (delta względem poprzednika), `Phase1_canonical_gate.py`
+ output, `Phase2_relax_lattice3d.py` (kopia+ewentualna delta) + output
+ `Phase2_backgrounds3d.npz`, `Phase3_bloch3d.py` (kopia+delta) + output
+ json, warunkowo `Phase4_nonlinear3d.py`, `Phase_FINAL_close.md`,
`NEEDS.md`, `README.md`.

## 6. Drzewo decyzyjne
```text
P1b-kan niespójny (brak zbieżności kotwicy) → STOP (raport)
P1c-kan FAIL → STOP; NEEDS: opcja c poprzednika (bazy sferyczne/O_h/AMR)
P2 brak istnienia wszędzie → CLOSED-GATE-STOP (negatyw istnienia sieci sc
   w modelu kanonicznym; NEEDS: bcc/fcc / inna geometria)
Q-PASS → NEEDS: „stabilność = własność konfiguracji 3D o skończonej
   gęstości" (user-gate: retrospektywa + sek08b) + PILNE konsekwencje
   dla drabiny mas
Q-FAIL → NEEDS: hipoteza stabilizacji gęstością zamknięta także w 3D sc
   (model kanoniczny); decyzja aksjomatyczna o znaku W bez podpory
Q-INCONCLUSIVE → NEEDS: wniosek metodologiczny, bez ontologii
```

---

**LOCK ZAMKNIĘTY 2026-08-31. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
