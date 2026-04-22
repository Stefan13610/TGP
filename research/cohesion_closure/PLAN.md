# cohesion_closure — program badawczy

**Session:** 2026-04-21
**Cel:** Sprawdzić czy TGP ma predykcyjny aparat dla energii kohezji metali
(E_coh ∈ [0.67, 8.9] eV/atom). Zamknąć lukę L3 z ATOMIC_SHELLS_VERDICT:
most między izolowanym atomem (atomic_shells FAIL) a metalem (Tc PASS w SC).

---

## Motywacja

Atomic_shells_closure (as01–as03) pokazał że TGP **nie ma** aparatu atomowego
(A_s, A_sp, A_d, A_f to empiryczne markery Fermiego, nie derywowane).
EM_VERDICT (em01, em02) pokazał że TGP **ma** aparat EM z substratu.

Pytanie ukryte: gdzie TGP kończy (fundamentalne) a gdzie zaczyna (fenomeno-
logiczne)? Energie kohezji to OSTATECZNY test — są to obserwable wielocząstkowe
dla ciał stałych, z jednoznacznie zmierzonymi wartościami.

- Li: 1.63, Na: 1.11, K: 0.93, Rb: 0.85, Cs: 0.80 eV (alkalie)
- Cu: 3.49, Ag: 2.95, Au: 3.81 (coinage)
- Fe: 4.28, Ni: 4.44, W: 8.90, Mo: 6.82, Nb: 7.57 (przejściowe)
- Al: 3.39, Hg: 0.67, Pb: 2.03 (ciekawe outliery)

Hipotezy robocze:

**H1. E_coh ∝ A_orb²** — A_s (alkalie), A_d (coinage+TMs), A_f (lantanoids)
z SC mogłyby skalować amplitudę Fermiego.

**H2. Koide dla E_coh** — K_coh = Σ E_coh / (Σ √E_coh)² = const dla rodzin
(alkalie, coinage, 3d TMs). Jak Koide dla leptonów (2/3) albo dla barionów.

**H3. Fermi-sea prosta formuła** — E_coh ≈ (3/5)·E_F·Z_val — (1/2)·e²/(ε_0·r_s·4π)
gdzie Z_val = walencja, r_s = promień Wignera-Seitza, E_F = energia Fermiego.

**H4. TGP soliton-like tail model** — E_coh = binding energy of electron
soliton w substracie kompresyjnym (metalu), skalujące jak m_e·c²·α²·n_coord.

## Dane bazowe (z Kittel + CRC)

Tabela E_coh (eV/atom) dla 25+ metali w coh00.

## Skrypty

| # | Skrypt | Cel |
|---|---|---|
| coh00 | baseline diagnostic | tabela E_coh, wzorce, kategoryzacja |
| coh01 | test H1: A_orb → E_coh | regresja E_coh = f(A_s, A_d, A_f) |
| coh02 | test H2: Koide-type | K_coh dla rodzin; residua |
| coh03 | test H3: Fermi-sea baseline | prosty model jeli-Ashcroft |
| coh04 | test H4: TGP soliton ansatz | czy m_e c²α² skala pracuje? |
| VERDICT | podsumowanie | które hipotezy przetrwały, luki |

## Kryteria falsyfikacji

- **H1 PASS** jeśli korelacja E_coh vs |A_orb|² ma r²>0.7 dla rodziny
- **H2 PASS** jeśli K_coh ma CV<5% w rodzinie (vs 10⁻¹⁵ dla leptonów)
- **H3 PASS** jeśli Fermi-sea model daje E_coh z błędem <30% dla alkali
- **H4 PASS** jeśli TGP-native formuła przewiduje trend (monotoniczność) dla
  alkali Li→Cs

**Scenariusz „zimny finał"**: wszystkie 4 hipotezy fail → TGP nie ma aparatu
do chemii stałej, cohesion_closure ma VERDICT taki jak atomic_shells
(uczciwy negatywny).

**Scenariusz „gorący finał"**: H1 albo H2 pass → TGP ma dodatkowy uchwyt
wielociałowy, można iść dalej (przewodnictwo termiczne, moduły sprężystości).
