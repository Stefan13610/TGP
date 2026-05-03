<!--
NEEDS.template.md
Co JEST do zrobienia w tym folderze. Otwarte luki, blokery, pytania.

Reguły:
  1. Każdy item musi mieć WERYFIKOWALNE źródło (`source:`) cytujące plik
     wewnątrz folderu LUB cytat z audytu / sąsiedniego folderu. Jeśli
     `source:` brzmi "general intuition", item NIE jest legalny.
  2. NEEDS.md NIE jest miejscem na nową fizykę / wyprowadzenia. To indeks
     tego, co już zostało gdzieś wskazane jako brakujące.
  3. Po domknięciu luki — przenosimy ją do FINDINGS.md (jako odpowiedź)
     lub usuwamy z listy z notatką w `## Closed needs` poniżej.
  4. Synchronizacja: `tgp_status.open_bridges` w README ma zawierać ID
     najważniejszych otwartych potrzeb (np. "BR-N1", "M9.1''-volume").
-->

---
title: "NEEDS — <nazwa folderu>"
date: <YYYY-MM-DD>
parent: "[[README.md]]"
type: needs
tgp_owner: research/<nazwa-folderu>
tags:
  - needs
  - <topic-tags>
---

# NEEDS — <nazwa folderu>

## Otwarte luki

| ID | Luka (1 zdanie) | Typ | Source | Kandydat dostawcy |
|----|------------------|-----|--------|-------------------|
| N1 | <co dokładnie brakuje> | analytical-bridge / numerical / data / lemma / decision | <plik:linia, np. `program.md §2.3` lub `[[../op-X/Phase2_results.md]] §UV2.5`> | <`research/<dostawca>` lub "open"> |
| N2 | … | … | … | … |

## Pytania otwarte

- **Q1**: <pytanie wymagające decyzji człowieka>
  - Source: `<plik>`
  - Wpływa na: <lista folderów / sekcji core>
- **Q2**: …

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | <co blokuje postęp>  | <ID innej luki / decyzja>  | <plik> |

## Closed needs (historia)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| (puste przy pierwszym zapisie) |  |  |  |
