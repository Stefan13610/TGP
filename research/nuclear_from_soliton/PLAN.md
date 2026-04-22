# nuclear_from_soliton — program badawczy

**Session:** 2026-04-21
**Cel:** Uczciwy test czy TGP ma aparat dla fizyki jądrowej. Rozszerzenie łańcucha
atom_from_soliton (wodór PASS) na jądra.

---

## Motywacja

Użytkownik wskazał pierwotnie:
> „Nawet niestabilność jąder atomowych powinna wynikać z samych równań solitonowych."

Po wykazaniu że wodór emerguje z TGP (atom_from_soliton, 13/13 PASS), naturalne
pytanie: czy łańcuch **substrat → solitony → wieloczastkowe stany związane**
rozciąga się na skale MeV (jądra)?

## Co TGP MA

### Infrastruktura wielociałowa (nbody/)

**V₂ (2-body)** — [[nbody/pairwise.py]]:
```
V₂(d) = -4πC₁C₂/d + 8πβC₁C₂/d² - 24πγC₁C₂(C₁+C₂)/(2d³)
```
Derywowane z Yukawa overlap integrals. Zweryfikowane analitycznie + numerycznie.

**V₃ (3-body irreducible)** — [[nbody/three_body_force_exact.py]]:
```
V₃(i,j,k) = (2β - 6γ)·C_i·C_j·C_k·I_Y(d_{ij}, d_{ik}, d_{jk})
```
gdzie I_Y = triple Yukawa overlap integral (Feynman parametryzacja, exact
numerycznie via 2D simplex quadrature).

**Mass scale** — m_sp² = γ (z vacuum condition β=γ). Wartość *dimensionless*
w jednostkach TGP. **Wartość w MeV nie jest przewidywana** — wymagałaby
fiksacji Φ₀.

### Istniejące testy

- [[nbody/tgp_bound_state.tex]]: trzy-ciałowe stany związane z V₂+V₃
- [[nbody/tgp_efimov_analog.tex]]: analog Efimova
- Three-body binding w reżimie kwantowym przy C ≲ C_Pl

## Czego TGP NIE MA (stan obecny)

1. **Nucleon jako topologia TGP** — brak identyfikacji: jaka winding/defekt
   to proton? neutron? Topologia leptonów (sek05-07) jest znana, ale
   baryonów nie.

2. **Masa nukleonu m_p** — TGP przewiduje stosunki masy leptonów (e:μ:τ =
   1:3477:...). m_p = 938 MeV wymagałoby innej topologii — otwarte.

3. **Masa pionu m_π** — TGP nie daje pionu jako Goldstone-bozonu łamanej
   symetrii chiralnej. To osobne zagadnienie sektor silny.

4. **SU(3) confinement** — TGP ma SU(3) defekty (sek10), ale nie ma
   numerycznej implementacji QCD-like confinement.

## Hipotezy robocze

**H1. Struktura V₂(TGP) pasuje do NN potencjału**
Czy można sfitować parametry (C, β, γ, m_sp) tak żeby V₂(r) dało deuteron
z E_B = 2.22 MeV? Jeśli tak, forma V₂ jest strukturalnie dopuszczalna dla
fizyki jądrowej.

**H2. V₃ (TGP) ma właściwy znak i skalę dla 3-ciałowych sił jądrowych**
W fizyce jądrowej trzy-ciałowa siła jest ZNANA potrzeba (Urbana IX, Illinois).
Bez niej ³H binding (8.48 MeV) jest źle przewidywany z samych V₂. Jeśli
TGP V₃ (ze sfitowanymi parametrami V₂) daje POPRAWNY trójciałowy wkład,
to jest nietrywialna weryfikacja strukturalna.

**H3. Alfa cząstka ⁴He jako ciasno związany 4-solitonowy stan**
⁴He ma B.E. = 28.3 MeV, znacząco więcej niż ³H per-nucleon (2.83 vs 2.83 MeV).
To jest klasyczny „alpha clustering" fenomen. Czy V₂+V₃ go odtwarza?

**H4. Niestabilność jąder z wielociałowych równań**
Docelowo: rozpad β, tunelowanie α. Wymaga sektora słabego (SU(2) dla W boson).
**Odłożone — nie w tej sesji**.

## Skrypty

| # | Skrypt | Cel | Zakres sesji |
|---|---|---|---|
| nfs00 | PLAN.md | ten plik | ✓ |
| nfs01 | deuteron_two_body.py | H1: fit V₂ do E_d=2.22 MeV | ✓ |
| nfs02 | triton_three_body.py | H2: predykcja E(³H) z V₂ + V₃ | ✓ |
| nfs03 | alpha_four_body.py | H3: ⁴He variational | opcjonalnie |
| nfs04 | beta_decay_sketch.py | H4: rozpad β z SU(2) | odłożone |

## Kryteria falsyfikacji

- **H1 PASS**: istnieje fizycznie sensowny set parametrów (m_sp ~ m_π range)
  dający E_d = 2.22 MeV ± 0.5 MeV
- **H2 PASS**: z tych samych parametrów V₂ + TGP V₃ daje E(³H) w zakresie
  6-10 MeV (observed 8.48)
- **H3 PASS**: ⁴He B.E. w zakresie 20-40 MeV
- **H1 FAIL, H2 FAIL**: struktura V₂+V₃ nie pasuje do nuclear scales —
  fundamentalna luka

## Scenariusze

**„Konstruktywny"**: H1+H2 PASS → TGP multi-body machinery jest strukturalnie
**kompatybilna** z fizyką jądrową. Brakuje tylko „mostu" (identyfikacja
nucleon↔topologia, m_π z chiral symmetry breaking). TGP jest *bliskie*
jądrom.

**„Pokorny"**: H1 PASS (można zawsze sfitować 1 punkt) ale H2 FAIL (trójciałowe
nie pasuje) → V₃ (TGP) ma zły znak/skalę dla nuclear. Luka głębsza.

**„Zimny"**: H1 FAIL → nawet deuteron niedostępny z V₂(TGP). Fundamentalna
niedoskonałość — forma V₂ nie jest NN potential.

## Ważne ograniczenia

**NIE** pretendujemy do:
- Derywacji m_p z TGP (niemożliwe obecnie)
- Derywacji m_π (wymaga chiral condensate)
- Pełnej fizyki jądrowej (wymaga SU(3) i SU(2) sektorów, nie zrobione)

**TYLKO** testujemy:
- Czy STRUKTURA V₂ (forma funkcji) jest dopuszczalna dla NN?
- Czy STRUKTURA V₃ (znak, skala) jest spójna z potrzebami 3-body jądrowego?

To jest **test strukturalnej gotowości** (structural readiness check),
nie test predykcji absolutnych wartości.

---

## Pliki sesji (do utworzenia)

- `PLAN.md` — ten plik
- `nfs01_deuteron_two_body.py` — test H1
- `nfs02_triton_three_body.py` — test H2
- `NUCLEAR_FROM_SOLITON_VERDICT.md` — werdyk

## Powiązania

- [[research/atom_from_soliton/ATOM_FROM_SOLITON_VERDICT.md]] — łańcuch atom
- [[nbody/tgp_bound_state.tex]] — teoria 3-body w TGP
- [[nbody/pairwise.py]] — implementacja V₂
- [[nbody/three_body_force_exact.py]] — implementacja V₃ (Feynman exact)
- [[sekcje_tex/sek10_N0_wyprowadzenie.tex]] — m_sp definicja
