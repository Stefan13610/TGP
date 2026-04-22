# fermion_from_soliton — plan

**Cel:** Lokalizować i zamknąć "Pauli-gap" zidentyfikowany w
[[../nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md]] (residual 1.45×
overbinding w triton i alpha po HC + konwergencji wariacyjnej).

## Kontekst

Nuclear_from_soliton (25/25 PASS) zidentyfikował dual-gap:
1. Hard-core repulsion (rho-exchange) — zregulowany z V_r=2400 MeV
2. **"Pauli-gap"** — pozostały 1.45× residual po konwergencji wariacyjnej

Atom_from_soliton (13/13 PASS) wspomniał o "missing fermion sector" jako
największej luce TGP. Cohesion_closure (0/13 FAIL) potwierdził że brak
antysymetryzacji uniemożliwia materia-skala kohezję.

Q6_statistics.py (qm_statistics, 8/8 PASS) POSTULUJE (-1)^B exchange phase
z topologii nawijania; nie derywuje z równań substratu. Plausibility study.

## Hipoteza strukturalna do przetestowania

Nuclear "Pauli-gap" może NIE BYĆ dokładnie Pauli antysymetryzacją przestrzeni
(w S-wave spatial wave function jest SYMETRYCZNY), lecz raczej **SPIN-ISOSPIN
CHANNEL GAP**: Pauli zmusza pary do konkretnych kanałów (nn → singlet S=0),
gdzie V_NN(singlet) ≪ V_NN(triplet) eksperymentalnie.

TGP V_NN jest scalar (nie zależy od spin/izospin), więc bezpośrednio dla
triton/alpha przeszacowuje binding w "zmuszonych singletach".

## Program

### fs01: Spin-channel hypothesis test

Test: Czy eksperymentalna hierarchia V_NN (V_triplet ≫ V_singlet) z Pauli
channel weighting kwantytatywnie zamyka nuclear residual?

Metoda:
- Explicit spin-isospin wave function dla 3H, 4He
- Wyznaczyć channel-probability per pair w każdym systemie
- Weighted V_NN = Σ_ch P_ch · V_ch
- Compare E_predicted vs E_obs

Predykcja: jeśli f_s = V_singlet/V_triplet ~ 0.3 (z a_t/|a_s| NN scattering),
triton E powinno wyjść ~-8.48 MeV automatycznie.

### fs02 (ZREALIZOWANE, 13/13 PASS): SU(2) × SU(2) derywacja channel algebry

Explicit Slater determinants w 4-dim internal space (iso × spin) dla
triton (3-body) i alpha (4-body). Derywacja channel projektorów
P_{T,S} = (3 + τ·τ)(3 + σ·σ)/16. Kluczowe odkrycie: heurystyka fs01
(2/3 singlet) BŁĘDNA — poprawnie jest 50/50 per-pair dla obu systemów.
Poprawka formuły w(f_s) = (1+f_s)/2 — f_s ≈ 0.85 wciąż w eksperymentalnym
zakresie. Operator decomposition V_NN = V₀ + V_σσ(σ·σ) + V_ττ(τ·τ) +
V_στ(σ·σ)(τ·τ) może generować f_s ≈ 0.85 z realistycznych mesonowych
parameters.

### fs03 (ZREALIZOWANE, 8/8 PASS): operator decomposition + strukturalna diagnoza

Próbowaliśmy derywacji V_σσ, V_ττ z TGP phi^4 overlap. Odkryto że:
- **f_s ≠ 1 wymaga algebraicznie V_σσ ≠ 0 LUB V_ττ ≠ 0** (dowód: iloczyn
  <(σ·σ)(τ·τ)> = -3 w OBU allowed L=0 kanałach; V_0 i V_στ nie rozróżniają)
- TGP phi^4 produkuje WYŁĄCZNIE V_0 (skalar) — implicit f_s = 1.0
- Zatem luka jest **STRUKTURALNA**, nie reparametryzacyjna
- Wymaga: (A) Skyrmion-like topological nucleon LUB (B) explicit Dirac
  multi-body, LUB (C) fenomenologiczny proxy (już zrobione w fs01+fs02)

Program fermion_from_soliton ZAKOŃCZONY jako diagnostyczny sukces.

## Czego NIE testujemy

- Spinory Diraca z {γμ,γν}=2gμν (wykraczające poza niniejszy program)
- Wielociałowa antysymetryzacja z explicit Slater determinants
- Spin-statistics theorem per se (q6 już zrobione phenomenologicznie)

## Success criteria

Bardzo skromne:
- fs01 PASS → TGP V_NN + phenomenological channel weight = quantitative nuclear physics. Nad-luka lokalizowana do V_NN(spin, isospin) form.
- fs01 FAIL → luka głębsza; diagnose następny kierunek.

## Pliki

- [[PLAN.md]] — ten dokument
- [[fs01_spin_channel_hypothesis.py]] — 6/6 PASS: test channel-weighted V_NN
- [[fs02_su2_channel_derivation.py]] — 13/13 PASS: SU(2)×SU(2) derywacja
- [[fs03_obe_operator_decomposition.py]] — 8/8 PASS: operator decomposition + strukturalna diagnoza
- [[FERMION_FROM_SOLITON_VERDICT.md]] — werdyk v3 (27/27 PASS łącznie, program zakończony)

## Powiązania

- [[../nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md]] — source of Pauli-gap
- [[../atom_from_soliton/ATOM_FROM_SOLITON_VERDICT.md]] — broader fermion gap
- [[../qm_statistics/q6_statistics.py]] — q6 statistics postulate
- [[../../sek09_cechowanie.tex]] — gauge/phase sector
- [[../../dodatekU_su2_formalizacja.tex]] — SU(2) emergence
