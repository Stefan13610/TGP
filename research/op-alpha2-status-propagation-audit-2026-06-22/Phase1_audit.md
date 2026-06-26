---
title: "Phase 1 AUDIT — op-alpha2-status-propagation-audit (inwentaryzacja + klasyfikacja)"
date: 2026-06-22
type: phase_audit
status: 🟡 Phase 1 — tabela trafień + klasyfikacja (value-blind)
phase: 1
cycle: op-alpha2-status-propagation-audit-2026-06-22
anti_lakatos_lock: ACTIVE
---

# Phase 1 — AUDYT (klasyfikacja value-blind wg reguły §3 LOCK)

> Metoda: grep całego `TGP_v1/**/*.{md,tex}` po `α=2 / alpha2 / theorem / derived /
> algebraic / wyprowadzenie z substratu`, następnie klasyfikacja każdego trafienia
> w rdzeniu i dokumentach głównych/publikacyjnych. Pliki badawcze (`research/**`)
> i historyczne (`dodatek_pivot_history`, `audyt/**`) POZA zakresem (opisują własny czas).

## §1 — Tabela trafień (rdzeń + dokumenty główne/submission)

| # | Plik : linia | Cytat (skrót) | Klasa | Uzasadnienie wg reguły |
|---|---|---|---|---|
| H1 | `README.md` (Highlights) | „Kinetic coupling **α = 2: algebraic theorem, not a fit**" | **DO-POPRAWY** | Claim derywacji bez kwalifikatora selekcji; sprzeczny z F. |
| H2 | `README.md` (Status v2 pivot) | α = 2 „**established through** prop:substrate-action + classification theorem thm:alpha2" | **DO-POPRAWY** | „established through substrate-action" = derywacja; F mówi: selekcja na gęstości. |
| H3 | `tgp_letter.tex:44` | „The kinetic coupling α = 2 is **derived as an algebraic theorem**" | **DO-POPRAWY** | Submission (PRL). Jawna derywacja. |
| H4 | `tgp_companion.tex:55` | „α = 2 is **proved as an algebraic theorem from** [substrate]" | **DO-POPRAWY** | Submission (PRD). „from substrate" = dokładnie krok C3 obalony w op-A3. |
| H5 | `tgp_companion.tex:219` | „conformal argument **fixes** K = g^{2α} with α = 2" | **NIEJASNE** | Poprawne jako selekcja w klasie (R3), ale „fixes" bez cross-ref → ryzyko odczytu „derywacja". Lekki cross-ref. |
| H6 | `papers_external/tgp_core_paper/tgp_core.tex:64` | „a unique second-order kinetic operator with coefficient **α = 2 as an algebraic theorem**" | **DO-POPRAWY** | Submission (core paper). Derywacja bez kwalifikatora. |
| H7 | `tgp_core.tex:336` | „(C3) **fixes α = 2 through the** [substrate Hamiltonian]" | **DO-POPRAWY** | To jest dosłownie krok C3 z thm:D-uniqueness, którego micro/macro op-A3 obalił. |
| H8 | `tgp_core.tex:304` `thm:alpha2` | „Uniqueness of D_kin within class (C1)–(C3), **α = 2**" | **SPÓJNE** | Jednoznaczność *w klasie* — chronione R3. (Ale H6/H7 ramują ją jako „theorem from substrate" → tam DO-POPRAWY.) |
| H9 | `tgp_core.tex:139–141` | „a rigorous theorem. Numerical block-spin MC and FRG LPA' … to the α = 2 fixed point" | **NIEJASNE** | Hedge istnieje (MC/FRG numeryczne), ale „rigorous theorem" + brak cross-ref do op-A3 micro/macro caveat. |
| H10 | `core/_meta_latex/status_map.tex:72` | „K(φ) = K_geo φ⁴ ⇒ K(0)=0; **wyprowadzenie z substratu**" | **DO-POPRAWY** | Rdzeń-meta. „wyprowadzenie z substratu" sprzeczne z F (sek08/sek10 już mówią: selekcja). |
| H11 | `core/formalizm/dodatekH_lancuch_wyprowadzen.tex:113` | „Status N0-4: [AX]→[AN] (**wyprowadzony z A2**)"; prop:K0-from-substrate | **DO-POPRAWY** | Łańcuch wyprowadzeń prezentuje K=φ⁴/α=2 jako „wyprowadzony", nie selekcję. |

## §2 — Referencja docelowa spójności (SPÓJNE — potwierdzone, NIE ruszać)

| # | Plik : linia | Cytat (skrót) | Klasa |
|---|---|---|---|
| S1 | `sek08:1057+` `rem:alpha2-pivot-status-pl` | „thm:D-uniqueness jest wynikiem [jednoznaczności w klasie]… NIE wyprowadzeniem z substratu… α=2 ma status [selekcji]" | **SPÓJNE** (kotwica) |
| S2 | `sek10:208–216` | „Kanoniczny α=2 jest **postulowany na gęstości**… selekcja C1–C3, nie derywacja. Substrat NIE derywuje α=2 (dałby α=½)" | **SPÓJNE** (kotwica) |
| S3 | `TGP_FOUNDATIONS.md:60` | „K(φ) = K_geo φ⁴ (**z α=2 selection**)" | **SPÓJNE** |

## §3 — Podsumowanie liczbowe
- **SPÓJNE:** 4 (H8 + S1–S3).
- **DO-POPRAWY:** 7 (H1, H2, H3, H4, H6, H7, H10, H11) → *uwaga: 8 trafień; H1–H11 numerowane, policz: H1,H2,H3,H4,H6,H7,H10,H11 = **8 twardych**.*
- **NIEJASNE:** 3 (H5, H9; + status_map/tabela_epistemiczna etykieta [AN] do lekkiego cross-ref).

**Korekta licznika (value-blind):** DO-POPRAWY = **8** (H1, H2, H3, H4, H6, H7, H10, H11);
NIEJASNE = **2** (H5, H9). SPÓJNE = **4** (H8, S1, S2, S3).

## §4 — Obserwacja krytyczna (asymetria rdzeń vs submission)
Rdzeń manuskryptu (sek08/sek10) został naprawiony w #31/#32 → „selekcja na gęstości".
**Ale naprawa NIE spropagowała** do warstwy submission (letter, companion, core paper)
ani do meta-rdzenia (status_map, dodatekH). To jest **dokładnie ten sam wzorzec** co #35
(κ_E): cykl rozstrzygający status nie zaktualizował dokumentów pochodnych. Tu stawka wyższa,
bo trafienia H3/H4/H6 są w plikach przeznaczonych do recenzji zewnętrznej.

## §5 — Werdykt wstępny (Phase 1)
**DO-POPRAWY** (8 trafień twardych + 2 niejasne). Wynik NIE jest negatywny.
Propagacja statusu „α = 2 = selekcja na gęstości (nie derywacja)" z #31/#32 była **NIEPEŁNA**.
