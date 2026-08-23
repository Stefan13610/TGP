---
title: "POST_CLOSE_BRAINSTORM_UPDATE — korekta zakresu: od golego substratu, nie od solitonow"
date: 2026-07-04
type: post-close-brainstorm
status: BRAINSTORM-ADDENDUM
anti_lakatos_lock: "Phase0 not rewritten; this file records an off-protocol conceptual correction"
related:
  - "[[README.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
  - "[[../../axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex]]"
---

# Post-close brainstorm update

## 1. Korekta autora

Po zamknieciu `op-blocked-soliton-bang` autor doprecyzowal, ze poprzedni
toy modelowal zbyt pozny etap. Wlasciwy obraz nie zaczyna sie od listy
gotowych solitonow.

Poprawny lancuch:

```text
1. Goly substrat relacyjny Gamma.
2. Lokalne zaburzenie s_i != 0, jeszcze nie soliton.
3. Phi_i = s_i^2 > 0 wlacza lokalny zalazek metrycznosci.
4. W tej lokalnej metrycznosci rozpad/kreacja nowych zaburzen jest bezkosztowa.
5. Zaburzenia kaskadowo narastaja.
6. Przy duzej gestosci/overlapie dochodzi do locku.
7. Po locku zaburzenia mozna nazywac solitonami; kreacja i relaksacja zamykaja sie.
```

Ten op testowal tylko etap:

```text
gotowe solitony -> klaster -> blokowana relaksacja
```

czyli etap 6/7, nie geneze 1-5.

## 2. Zgodnosc z sekcja 1

To rozroznienie jest zgodne z `core/sek01_ontologia/sek01_ontologia.tex`:

- substrat `Gamma=(V,E)` jest poziomem 0 i nie ma metryki;
- `Phi=0` oznacza faze niemetryczna substratu, nie pusta przestrzen;
- `Phi>0` oznacza faze metryczna, w ktorej propagacja/geometria maja sens;
- `Phi0` jest pozniejszym tlem referencyjnym wygenerowanym przez globalna
  zawartosc materii, a nie stanem bez materii.

Z `axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex` wynika natomiast
minimalny lancuch:

```text
N0 -> Z2 -> s_i -> Gamma -> Phi -> g_mu_nu
```

oraz fakt, ze niestabilnosc `N0/Phi=0` jest strukturalna/logiczna/miarowa,
nie procesem w czasie fizycznym. Dlatego symulacja genezy musi uzywac
parametru relaksacji/selekcji `tau`, a nie zakladac gotowego czasu
metrycznego.

## 3. Wniosek metodologiczny

Nastepny model nie powinien startowac od:

```text
positions = [soliton_1, ..., soliton_N]
amplitudes = [...]
```

Powinien startowac od:

```text
graph Gamma
site variables s_i
Phi_i = s_i^2
no external metric
```

Minimalna klasa rownan:

```text
ds_i/dtau = - dH_Gamma/ds_i + seed_i + R_i[s]
Phi_i = s_i^2
metric_enabled_i = Phi_i > epsilon
```

gdzie `R_i[s]` jest kandydackim sprzezeniem relacyjnym/metrycznym, ktore
moze obnizac lokalny koszt kreacji kolejnych zaburzen, ale samo musi byc
funkcja konfiguracji substratu, nie arbitralnym spawn-rule.

## 4. Kryteria przyszlego op-a

Proponowany osobny cykl:

```text
research/op-bare-substrate-genesis-2026-07-04/
```

Pytania zalockowane przed kodem:

1. Czy goły substrat relacyjny bez seedu pozostaje niemetryczny?
2. Czy pojedyncze male zaburzenie `s_i != 0` zanika, rozlewa sie, czy
   odpala kaskade?
3. Czy lokalny koszt `DeltaE_create(x)` spada do zera w regionie
   `Phi>0`, zanim pojawi sie lock?
4. Czy kaskada konczy sie samorzutnym lockiem, czyli stanem, w ktorym
   dalsza kreacja i relaksacja sa niemozliwe?
5. Czy struktury po locku maja profile solitonopodobne bez wrzucania
   solitonow jako obiektow wejscia?

Minimalne falsyfikatory:

- kazde zaburzenie zawsze zanika do `Phi=0`;
- kazde zaburzenie zawsze runaway bez locku;
- lock pojawia sie tylko przez jawny clamp/spawn-rule;
- struktury koncowe nie sa lokalne ani solitonopodobne.

## 5. Status obecnego op-a po korekcie

`op-blocked-soliton-bang` pozostaje uzyteczny, ale jako test poznej fazy:

```text
czy gotowy klaster solitonow moze blokowac relaksacje i generowac energie interakcyjna?
```

Nie powinien byc cytowany jako test:

```text
czy solitony powstaja z golego substratu?
```

Tego jeszcze nie policzono.
