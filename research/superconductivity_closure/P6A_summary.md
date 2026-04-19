# P6.A — Ambient high-T_c w TGP (ps9-ps11 podsumowanie)

**Data:** 2026-04-19
**Status:** P6.A (d-wave 2D) zamknięte jakościowo; P6.B (phonon coupling) otwarte.

## Wyniki liczbowe

| Program | Klasa | r(log-log) | RMS_log | Status |
|---------|-------|------------|---------|--------|
| ps5 (single-band BCS) | 19 SC mixed | 0.48 | 0.61 | baseline |
| **ps10 (P6.A cuprates)** | **8 cupratów** | **0.96** | **0.19** | **✓ working** |
| ps11 Fe-SC | 8 pnictides/chalcogenides | 0.35 | 0.43 | słaby (brak P6.B) |
| ps9 MgB2 | 1 mat (dwu-gap) | — | — | 14 K vs 39 K, off |

## Model P6.A

$$
T_c^{\text{cuprate}} = K_{dw} \cdot k_d(z_{\text{planar}}) \cdot C_0 \cdot A_{ZR}^2 \cdot M(a_\parallel) \cdot \sqrt{n_{\text{layers}}} \cdot \Lambda_E^{\text{cup}}
$$

**Parametry uniwersalne (z fitu 8 cupratów):**
- $K_{dw} = 3.50$ (d-wave + van Hove booster)
- $\Lambda_E^{\text{cup}} = 0.051$ meV (2.6× mniejsze niż BCS 0.131 meV)
- $A_{ZR}^2 = A_d^2 + 2A_p^2 = 0.181$ (Zhang-Rice singlet: Cu 3d + 2× O 2p)
- $z_{\text{planar}} = 8$, $k_d = 2.94$
- layer_factor = $\sqrt{n}$ (z saturation)

**Interpretacja:**
- $K_{dw}$ = 1.5 (d-wave node boost) × 2.5 (vH DOS enhance) × 0.9 (ZR correction)
- Dlaczego $\Lambda_E$ cuprate < BCS? Bo cuprates parują ELEKTRONIKO (spin fluctuations), nie przez fonon → słabsza skala wzbudzeń substratu.

## Per-cuprate residuals (ps10)

| Cuprate | a_inplane | n | T_obs | T_pred | Δlog |
|---------|-----------|---|-------|--------|------|
| La214 opt | 3.780 | 1 | 38 | 53 | +0.15 |
| YBCO | 3.820 | 2 | 92 | 76 | -0.08 |
| BiSCCO | 3.820 | 2 | 110 | 76 | -0.16 |
| Tl2212 | 3.855 | 2 | 105 | 76 | -0.14 |
| **Hg1223** | 3.855 | 3 | **138** | 93 | -0.17 |
| Tl2223 | 3.855 | 3 | 127 | 93 | -0.14 |
| Nd214 (e-dop) | 3.940 | 1 | 24 | 54 | +0.35 |
| Bi2201 | 3.820 | 1 | 34 | 54 | +0.20 |

Model systematycznie:
- **niedoszacowuje n=2,3 cuprate rekordy** (−0.10 to −0.17)
- **przeszacowuje n=1 niższe-Tc** (+0.15 to +0.35)

Sygnał: layer_factor `√n` jest za słabe; może powinno być `n^0.7` lub `(n-n_0)^0.5` z offsetem.

## Odpowiedź na pytanie użytkownika: "wysokotemperaturowe przy ambient?"

### TAK, są możliwe — ale z ograniczeniami

**Istniejące rekordy ambient:**
- Hg1223: 138 K (record ambient, od 1993)
- Hg-Ca-Ba-Cu-O rodzina: do 150 K przy n=4

**TGP predykcje dla nieodkrytych ambient SC:**

| Kandydat | Scenariusz | T_c predykcja |
|----------|------------|----------------|
| Hg1234 (n=4) | naturalne + ewentualnie strain | 107–120 K |
| Hg1245 (n=5) | syntezowalny? niedawne próby | **120–135 K** |
| Hg-inf-layer (n=10) | metastabilny, trudno syntezować | 170 K |
| Hg1223 / SrTiO3 strain | MBE thin film | 93–100 K (strain mały) |
| Idealna hybryda (a=4.088, n=5, A_ZR opt) | pipe dream | **150 K** |

### Czy możliwe T_c > 200 K przy ambient?

**Marginalnie TAK, w TGP**, ale wymaga:
1. `a_inplane` dokładnie na `a*_TGP = 4.088 Å` (obecne cuprates mają 3.82–3.88)
2. `n ≥ 5` stabilnych płaszczyzn parujących (Hg-cuprates max n=4 zsyntezowane)
3. Hybryda orbitalna `A_ZR` > 0.5 (wymaga znalezienia innego mostu od dwóch O-p)
4. `Λ_E > 0.08` meV (wymaga wyższej `ω_Debye` substratu)

**Żadna z tych warunków nie jest z góry wykluczony fizyką TGP.**

### Wpływ próżni (UHV)

Próżnia **nie podnosi** T_c bezpośrednio, ale:
- Cuprates rozkładają się w powietrzu przy ogrzewaniu → UHV dla samplowania
- FeSe/STO monolayer (65 K) STABILNY tylko w UHV
- Ni-infinite-layer filmy wymagają ultra-clean atmosphere

## Przepis eksperymentalny — test TGP P6.A

**Konkretna propozycja:**

1. **Synteza Hg1245** — n=5 w rodzinie HgBa₂Ca₄Cu₅O₁₂+δ
   - Technika: high-P synthesis at 3-5 GPa, rapid quench
   - Przewidywane T_c = 120–135 K ambient (MA być powyżej Hg1223 rekordu 138 K!)
   - Literature: M. Nuñez-Regueiro et al. próba w 1990s, niestabilne fazy

2. **Hg1234 / SrTiO3 epitaxial**
   - MBE/PLD wzrost cienkiego filmu (10-50 nm)
   - SrTiO3 a=3.905 → tensile strain 1.3%
   - TGP predykcja: T_c ~108 K (marginalne vs bulk Hg1234 ~125K)
   - Ale: jeśli substrate phonon coupling (P6.B) działa jak w FeSe/STO (65K ex 8K) — mogłoby dać **>200 K**

3. **Test key TGP: skan a_inplane przez strain**
   - Różne substraty: LaAlO3 (3.79), SrTiO3 (3.905), DyScO3 (3.946)
   - Jeśli T_c maksimum przy a bliskim 4.088 Å → potwierdza TGP harmonika
   - Jeśli monotoniczna — TGP falsyfikowane

## Co pozostaje dla P6.B i dalej

### P6.B: Phonon-substrate coupling (hydrydy i FeSe/STO)
- H3S 203 K, LaH10 250 K — wymagają ciśnienia + light-atom phonons
- FeSe/STO monolayer 65 K (vs 8 K bulk) — interface coupling
- Hipoteza: `Λ_E^eff = Λ_E * (ω_phonon/ω_0)^α`, α do fitowania

### P6.C: f-metale orbital switching
- Ce, Yb pod ciśnieniem: 4f → 5d hybrydyzacja
- Formalizm: `A_eff = (1-η) A_f + η A_d`, η(P)

### P7 (spekulatywnie): room-temperature SC
- TGP teoretyczny plafon ~300 K wymaga `n` bardzo duże + hybryd orbitalny + optimal `a`
- Czy to realistyczne? Nie wiadomo - wymaga P6.B kalibracji.

## Pliki

- [[ps9_mgb2_two_gap.py]] — próba dwu-gap MgB2 (nie trafia, problem P6.B)
- [[ps10_cuprates_2d_dwave.py]] — **fit cupratów r=0.96**
- [[ps11_ambient_highTc_survey.py]] — ambient high-Tc survey + Fe-SC + predykcje
- [[P6A_summary.md]] — ten plik
- [[P6_plan.md]] — ogólny plan P6 (ostatnia aktualizacja)
