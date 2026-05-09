---
title: "PAPER_LAYOUT — TGP_v1 published documents map"
date: 2026-05-09
type: meta-documentation
status: ACTIVE
purpose: "Wyjaśnienie roli każdego z 3 PDF-ów w roocie + .tex źródeł"
related:
  - "[[STATE.md]]"
  - "[[README.md]]"
  - "[[main.tex]]"
  - "[[tgp_letter.tex]]"
  - "[[tgp_companion.tex]]"
---

# PAPER_LAYOUT — co to za 3 PDF-y w roocie

Częsty problem dla nowego kontrybutora / reviewera: w roocie są **3 różne PDF-y**.
Który jest "the paper"? Odpowiedź: **wszystkie trzy są kanoniczne**, każdy w innej
roli.

## Trójdzielny layout

| Plik źródłowy | PDF | Rozmiar | Język | Format LaTeX | Rola |
|---|---|---|---|---|---|
| `main.tex` | `main.pdf` | **5.4 MB** | 🇵🇱 Polski | `article` 11pt A4 | **Full thesis / monograph** — kompletny dokument referencyjny TGP, włącza wszystkie sekcje `core/sek*/` przez `\input{}`. To jest "the book" — autorska forma rozbudowana, ze wszystkimi rozdziałami, dowodami, dyskusjami filozoficznymi |
| `tgp_letter.tex` | `tgp_letter.pdf` | 312 KB | 🇬🇧 Angielski | `revtex4-2` PRL twocolumn | **Letter** (krótki paper) — "Three inputs, forty predictions: The Tensor-Gravitational-Potential as a flavor-cosmological meta-theory of the Standard Model". Format do złożenia w **Physical Review Letters** |
| `tgp_companion.tex` | `tgp_companion.pdf` | 455 KB | 🇬🇧 Angielski | `revtex4-2` PRD twocolumn | **Companion** (długi technical) — pełny techniczny artykuł towarzyszący letterowi. Format do złożenia w **Physical Review D** |

## Strategia publikacji (standard fizyka teoretyczna)

`tgp_letter.pdf` + `tgp_companion.pdf` to **typowa para journalowa**:

- **Letter** (PRL): krótka prezentacja głównych wyników (4 strony + figury), pisana jako "headline result" — czytelna dla szerokiej publiki
- **Companion** (PRD): pełna techniczna ekspozycja — wszystkie wyprowadzenia, sympy verifications, numerical work; czytelnik przychodzi po szczegóły

Ten model jest standardem dla dużych przełomowych wyników w fizyce teoretycznej.

## Pytanie "który jest 'the' paper?"

Odpowiedź zależy od kontekstu pytającego:

| Kto pyta? | Odpowiedź |
|---|---|
| Reviewer journalowy | **`tgp_letter.pdf`** + **`tgp_companion.pdf`** (para). Letter to entry point, companion to dowody. |
| Naukowiec czytający dla zrozumienia | **`tgp_companion.pdf`** (pełna techniczna ekspozycja, EN) lub **`main.pdf`** (PL, jeszcze pełniej) |
| Autor przy edycji aksjomatów / dodawaniu sekcji | **`main.tex`** (Polish, monograph format — to gdzie żyje treść `core/sek*/`) |
| Cytowanie w innym paperze | Letter (PRL DOI gdy ukaże się) lub companion (PRD DOI). Main jest internal/Zenodo. |

## Aktualizacja PDF-ów

Wszystkie trzy są kompilowane lokalnie z `.tex`. **Nie commituj** zaktualizowanych
PDF-ów po małych zmianach źródła — generują się na żądanie. **Commituj PDF** tylko
przy:

- Zenodo deposit (immutable timestamp)
- Pre-submission journal version
- Major version (np. v2 po feedback)

Build artifacts (`*.aux`, `*.log`, `*.bbl`, etc.) są w `.gitignore` — patrz linie 8-19.

## Status post-2026-05-09 falsification

⚠ **Wszystkie 3 dokumenty mają sekcje grawitacyjne (M9.1'') wymagające update'u
post-S07.** Konkretnie:

- `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`: CRITICAL UPDATE banner + dual-V annotation linie 95-126 (gravity-only deprecation)
- `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex`: CRITICAL UPDATE banner (M9.1'' (4-3ψ)/ψ falsified)

Żadna z PDF-ów nie powinna być deponowana / submitted **przed zamknięciem
[[research/op-S07-alternative-f-psi-derivation-2026-05-09/]]**.

## Bibliografie

- `tgp_main.bib` — used by `main.tex`
- `tgp_companionNotes.bib` (~104 bytes) — placeholder for companion bibliography (companion uses `tgp_main.bib` per `\bibliography{tgp_main}` w .tex)
- `tgp_letterNotes.bib` (~104 bytes) — placeholder analogiczne dla lettera

## Decyzja redakcyjna do podjęcia (osobna sesja)

Czy `tgp_companionNotes.bib` i `tgp_letterNotes.bib` są potrzebne? Jeśli letter
i companion używają `tgp_main.bib`, te pliki to noise — można usunąć.
Niskoprodukowy cleanup, dla porządku.
