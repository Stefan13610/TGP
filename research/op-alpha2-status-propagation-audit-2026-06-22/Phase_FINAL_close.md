---
title: "Phase FINAL — op-alpha2-status-propagation-audit. Werdykt value-blind: DO-POPRAWY (propagacja #31/#32 niepełna). 8 trafień twardych w warstwie submission + meta-rdzeń. Rekomendacje ZGŁOSZONE, edycje WSTRZYMANE do autoryzacji."
date: 2026-06-22
type: phase_close
status: 🟢 CLOSED — DO-POPRAWY naprawione (P1+P2+P3+P4 zastosowane; user autoryzacja 2026-06-22; main.tex build exit 0, 553 str.)
phase: FINAL
cycle: op-alpha2-status-propagation-audit-2026-06-22
anti_lakatos_lock: PRESERVED
verdict: "DO-POPRAWY — α=2 prezentowane jako 'derived/algebraic theorem from substrate' w README + 3 dokumentach submission + status_map + dodatekH, sprzecznie z zasądzonym (#31/#32) statusem 'selekcja na gęstości'. Rdzeń sek08/sek10 + FOUNDATIONS już spójne."
---

# Phase FINAL — CLOSE (op-alpha2-status-propagation-audit)

> **Werdykt (value-blind, reguła §3 LOCK):** **DO-POPRAWY** — propagacja statusu
> „α = 2 = selekcja aksjomatyczna na gęstości (NIE derywacja z substratu)" z cykli
> #31 (op-A3, sympy 5/5) i #32 (audit gęstości) była **NIEPEŁNA**.

## §1 — Werdykt
- **8 trafień twardych (DO-POPRAWY):** H1, H2 (`README.md`), H3 (`tgp_letter.tex`),
  H4 (`tgp_companion.tex`), H6, H7 (`tgp_core.tex`), H10 (`status_map.tex`),
  H11 (`dodatekH`). Wszystkie głoszą α = 2 / K = φ⁴ jako „derived / proved / algebraic
  theorem / wyprowadzone **z substratu**" bez kwalifikatora selekcji.
- **2 trafienia niejasne (NIEJASNE):** H5 (`companion:219` „fixes"), H9 (`core:139–141`
  „rigorous theorem" + hedge MC/FRG) — poprawne jako jednoznaczność w klasie, ale brak
  cross-ref do `rem:alpha2-pivot-status`.
- **4 spójne (SPÓJNE, kotwice):** H8 (thm:alpha2 jako jednoznaczność w klasie C1–C3),
  S1 (`sek08 rem:alpha2-pivot-status-pl`), S2 (`sek10:208–216`), S3 (`FOUNDATIONS`).

## §2 — Anti-Lakatos (checklist)
- ✓ **R1 chronione:** relacja EL α = p/2 NIE podważona — formalizm OK.
- ✓ **R2 chronione:** α = 2 NIE odrzucone fenomenologicznie (PPN/masy/Koide stoją);
  audyt dotyczy WYŁĄCZNIE statusu epistemicznego (selekcja vs derywacja).
- ✓ **R3 chronione:** thm:alpha2 jako jednoznaczność *w klasie* = SPÓJNE (H8);
  DO-POPRAWY tylko tam, gdzie dopisano „from substrate".
- ✓ **R4 chronione:** Φ = ⟨ŝ²⟩ kanoniczne — Opcja B NIE reaktywowana.
- ✓ Werdykt **wyliczony** z reguły (8 twardych grep-trafień), nie wybrany.
- ✓ Rdzeń/papers **NIE edytowane** (forbidden move #4) — rekomendacje §3 zgłoszone osobno.
- ✓ Budżet nowych stałych = 0.

## §3 — Rekomendacje (ZGŁOSZONE; NIE wykonane — czekają na autoryzację)

### P1 — README.md (2 edycje, addytywne)
- H1: „α = 2: algebraic theorem, not a fit" → „α = 2: **selekcja aksjomatyczna na gęstości**
  (jednoznaczność w klasie C1–C3; NIE wyprowadzona z substratu — substrat dałby α = ½,
  por. op-A3 / `rem:alpha2-pivot-status`)".
- H2: „established through prop:substrate-action…" → przeramować na „**selected within class
  C1–C3** (consistency, not substrate derivation)".

### P2 — Dokumenty submission (letter / companion / core paper) — NAJWYŻSZY priorytet
- H3 (`tgp_letter.tex:44`), H4 (`tgp_companion.tex:55`), H6 (`tgp_core.tex:64`),
  H7 (`tgp_core.tex:336`): zamienić „derived / proved / fixed … **from the substrate**"
  → „**fixed by the conformal selection conditions (C1)–(C3)** [axiomatic selection on the
  density Φ = ⟨ŝ²⟩, not derived from the microscopic substrate]".
- H8/thm:alpha2: zostawić jako jednoznaczność w klasie, ale dodać zdanie-caveat + cytat op-A3.
- H5, H9: lekki cross-ref do `rem:alpha2-pivot-status`.

### P3 — Meta-rdzeń
- H10 (`status_map.tex:72`): „wyprowadzenie z substratu" → „**selekcja C1–C3 (na gęstości)**".
- H11 (`dodatekH:113`): status „[AX]→[AN] (wyprowadzony z A2)" → doprecyzować, że N0-4
  (K = φ⁴ / α = 2) to **selekcja**, nie derywacja micro→macro (cross-ref op-A3 §2).

### P4 — Deliverable wtórny (parameter-counting; OSOBNA decyzja rdzenia, §7 LOCK)
Zgłoszenie wprost: jeśli α = 2 jest selekcją aksjomatyczną, to nagłówek „**40 predykcji
z 3 inputów**" zależy od decyzji, czy warunki C1–C3 liczą się jako wejście. **Audyt tego
NIE rozstrzyga** — rekomenduje jawną notę w README/abstrakcie: „inputs (g₀ᵉ, Ω_Λ, N)
**plus structural selection axioms** (Z₂, conformal class C1–C3, β=γ)". Analog do uczciwego
+1 C_σ w #35.

## §4 — Edycje ZASTOSOWANE (user autoryzacja 2026-06-22: „Pełny P1+P2+P3" + „Tak, dodaj notę" P4)

| # | Plik | Zmiana |
|---|---|---|
| P1a | `README.md` (Highlights) | „algebraic theorem, not a fit" → „axiomatic selection on the density; uniqueness in class C1–C3; NIE substrate derivation (substrat: α=½, op-A3, rem:alpha2-pivot-status)" |
| P1b | `README.md` (v2 pivot) | „established through prop:substrate-action…" → „**selected** … unique element of conformal class C1–C3 on the density; not derived from substrate" |
| P4 | `README.md` (Abstract) | dodana **Parameter-counting note**: 3 inputy numeryczne **+ structural selection axioms** (Z₂, klasa C1–C3 fixująca α=2, β=γ); honest headline |
| P2a | `tgp_letter.tex:44` | „derived as an algebraic theorem from the substrate" → „fixed as unique element of conformal class C1–C3 on density; axiomatic selection, not substrate derivation (α=½)" |
| P2a' | `tgp_letter.tex:38` | dodany clause „+ structural selection axioms (Z₂, conformal class fixing α=2, β=γ)" |
| P2b | `tgp_companion.tex:55` | „proved as an algebraic theorem from the substrate" → „fixed as unique element of conformal class C1–C3 on density; axiomatic selection, not derivation; RG confirms stability of selected fixed point" |
| P2b' | `tgp_companion.tex:219` | dodany caveat „axiomatic selection on the density, not microscopic-substrate derivation (substrate: α=½)" |
| P2c | `tgp_core.tex:64` | „α=2 as an algebraic theorem" → „fixed by conformal class C1–C3 as axiomatic selection on density (uniqueness within class; not substrate derivation, α=½)" |
| P2c' | `tgp_core.tex` (po thm:alpha2) | dodany `\begin{remark}[Epistemic status of α=2]` (selekcja w klasie, nie derywacja; α=2 fenomenologicznie wymagane) |
| P3a | `status_map.tex:72,77` | „wyprowadzenie z substratu" → „selekcja w klasie konforemnej C1–C3 (na gęstości), NIE derywacja (substrat: α=½)" |
| P3b | `dodatekH:113` | status N0-4 przeramowany: K=φ⁴/α=2 = selekcja w klasie na gęstości, NIE derywacja micro→macro (krok zakłada reprezentację amplitudową); cross-ref op-A3 |

## §5 — Build verification
- **`main.tex` (rdzeń, zawiera P3a status_map + P3b dodatekH):** `pdflatex` **exit 0, 553 strony** (2 przebiegi).
  **Brak NOWYCH dangling refs** (edycje użyły prozy, nie `\ref`). Pre-existing residual #32
  (`ax:substrat`, `ssec:disformal`, `app:A-aksjomaty`, `para:basin-stability`) — NIE z tego cyklu.
- **Papers (letter/companion/core):** moje edycje α=2 **składniowo poprawne**.
- **🟢 FINDING „papers" RESOLVED (2026-06-22, user autoryzacja „naprawa makr papers"):** pre-existing
  fatal errors makr (niezwiązane z α=2) **naprawione**, wszystkie trzy kompilują się do PDF (0 błędów):
  - **Korzeń `\gone`:** `\newcommand{\gone}{g_0^e}` kończy się indeksem górnym `^e` ⟹ każde `\gone^*`,
    `\gone^2`, `\gone^{(0)}` = *double superscript*. Fix globalny: opakowanie w nawiasy
    `\newcommand{\gone}{{g_0^e}}` (companion + letter) — naprawia **wszystkie** użycia bez dotykania
    dziesiątek miejsc.
  - **letter:** brak `\ZZ` → dodane `\newcommand{\ZZ}{\mathbb{Z}}`; `\gone` braces. ⟹ `tgp_letter.pdf` 4 str.
  - **companion:** brak `\usepackage[T1]{fontenc}` (l.527 polskie „1-p\k{e}tlowe" `\k{e}` w OT1) → dodane;
    `\gone` braces. ⟹ `tgp_companion.pdf` 11 str.
  - **core:** `$\Z2$` (l.374, `\Z` niezdefiniowane) → `$\Zp$` (zdef. `\mathbb{Z}_2`). ⟹ `tgp_core.pdf` 11 str.
  - Jeden błąd wprowadzony przeze mnie wcześniej (`\Zp` w letter) — naprawiony (`\mathbb{Z}_2`).
  - Pozostałe: undefined citations (wymaga rozwiązania bibliografii — patrz niżej).
- **🟢 BIBLIOGRAFIA + ODNOŚNIKI DOMKNIĘTE (2026-06-22, user „tak działaj"):** wszystkie trzy papers
  używają **wbudowanego `\begin{thebibliography}`** (ręczna lista `\bibitem`), więc **bibtex zbędny** —
  wystarczają 2–3 przebiegi `pdflatex`. Stan końcowy (0 błędów fatalnych, **0 undefined citations,
  0 undefined references**):
  - `tgp_letter.pdf` — **4 str.** (czysto po dodaniu `\ZZ`).
  - `tgp_core.pdf` — **12 str.**; naprawiony dangling `\ref{prop:substrate-action}` (etykieta tylko w
    Dodatku B pełnego manuskryptu, brak w standalone) → zamieniony na tekst „the substrate-action proposition".
  - `tgp_companion.pdf` — **14 str.**; dwa szerokie `table*` (`tab:roadmap`, `tab:all-predictions`) były
    **gubione** (float dwukolumnowy `[t]` za wysoki / kolejka deferred zablokowana). Fix: relaksacja
    parametrów floatów dwukolumnowych w preambule (`\dbltopfraction=0.95`, `\dblfloatpagefraction=0.05`,
    `dbltopnumber=4`) + `[tp]` + `\footnotesize` na dużej tabeli ⟹ obie tabele **umieszczone, oba `\ref`
    rozwiązane**. Residual: 1× benign revtex „Deferred float stuck during \clearpage" (kosmetyczne — float
    i tak umieszczony, co potwierdza 0 undefined references).

## §6 — Status końcowy
**🟢 CLOSED — DO-POPRAWY naprawione (P1+P2+P3+P4).** Łuk **#31** (derywacja obalona, sympy 5/5)
→ **#32** (gęstość kanoniczna) → **#36** (propagacja statusu): rdzeń + warstwa submission + meta-rdzeń
teraz **spójne** ze statusem „α=2 = selekcja aksjomatyczna na gęstości, nie derywacja z substratu".
Nagłówek „40 z 3 inputów" uczciwie zakwalifikowany (+ structural selection axioms). **Osobny finding:**
3 papers mają pre-existing fatal errors makr (do osobnej decyzji).
