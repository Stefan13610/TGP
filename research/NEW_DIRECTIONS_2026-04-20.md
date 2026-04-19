# Nowe kierunki badawcze TGP — 2026-04-20

Pięć nowych folderów badawczych otwartych po zamknięciu P7.2 (nadprzewodnictwo)
i wypuszczeniu papers rdzeniowego + SC na Zenodo.

## Spis

| # | Folder | Obszar | Horyzont | Flagship test |
|---|--------|--------|----------|---------------|
| 1 | [[thermal_transport_molecular/PLAN.md]] | kondensowana materia → transport ciepła | średni | plateau $\lambda_\text{th}$ w rubrenach/perowskitach |
| 2 | [[casimir_mof/PLAN.md]] | próżnia substratowa → nanogeometria | średni | ciśnienie Casimir w MOF-5/177 |
| 3 | [[neutrino_msw/PLAN.md]] | leptony → oscylacje neutrin | długi | masa bezwzględna $\sum m_\nu$ |
| 4 | [[liquid_viscosity/PLAN.md]] | miękka materia → transport pędu | średni | VFT/Adam-Gibbs unifikacja |
| 5 | [[muon_g_minus_2/PLAN.md]] | leptony → poprawki QED | **krótki (2 lata)** | $\Delta a_\mu, \Delta a_e, \Delta a_\tau$ z jednego parametru |

## Wspólna filozofia

Każdy folder ma identyczną strukturę `PLAN.md`:

1. Problem fizyczny + state of the art
2. Dlaczego TGP (hipoteza mechanizmu)
3. Cele badawcze (T1-T4 / N1-N4 / G1-G5 itp.)
4. Plan numeryczny (ps01, ps02, ...)
5. Literatura startowa
6. Relacje z innymi sektorami TGP
7. Falsyfikowalność (konkretna predykcja)
8. Otwarte pytania
9. Link do rdzenia TGP

## Priorytetyzacja (subiektywna)

**Krótki horyzont (0-2 lata), najwyższa falsyfikowalność:**
- **(5) muon g-2** — FNAL Run-6, J-PARC E34: wynik do 2028 z $10^{-11}$
- **(2) Casimir MOF** — pomiary bulk modulus MOF już istnieją

**Średni horyzont (2-5 lat):**
- **(1) thermal transport** — dane dla perowskitów hybrydowych rosną szybko
- **(4) viscosity** — dostępne dane VFT dla setek cieczy, dopasowanie od razu

**Długi horyzont (5+ lat):**
- **(3) neutrino MSW** — następna galaktyczna SN + CMB-S4

## Zalecany następny krok

Najlepszy stosunek "koszt/efekt" — czyli coś co za niewielką pracę da duży
sygnał do uwiarygodnienia TGP:

- **liquid_viscosity** — bazy danych VFT są ogromne, fit o 2-3 parametrach do
  50+ cieczy daje natychmiastową walidację.
- **muon_g-2** — wymaga rozumienia QED-pętli w TGP-metryce (techniczne), ale
  publikowalny test niezależny.

## Relacje z papers

- **Paper 1 TGP Core** [[papers/core/tgp_core.pdf]] — fundament dla wszystkich.
- **Paper 2 TGP-SC** [[papers/sc/tgp_sc.pdf]] — wzór struktury aplikacyjnej.
- Planowane kolejne preprints Zenodo (w kolejności gotowości):
  - `papers/viscosity/tgp_viscosity.tex` — jeśli VFT unifikacja wyjdzie
  - `papers/leptons/tgp_leptons.tex` — mass + Koide + Cabibbo + g-2 unified
  - `papers/thermal/tgp_thermal.tex` — anomalous heat transport

Każdy preprint cytuje Paper 1 jako `[TGP-core]` i inne wcześniejsze papers
(e.g. `[TGP-SC]`).
