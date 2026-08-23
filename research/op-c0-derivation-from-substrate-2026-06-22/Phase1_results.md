---
title: "Phase 1 — op-c0-derivation-from-substrate. Rekoncyliacja analityczna: c₀ = normalizacja TEGO SAMEGO operatora ∂ŝ∂ŝ co C_σ ⟹ dziedziczy rozbieżność liniową UV (F1). Werdykt wyliczony (LOCK §4 wiersz 1): c₀ = WOLNY PARAMETR UV. Sympy 5/5 PASS."
date: 2026-06-22
type: phase_results
phase: 1
cycle: op-c0-derivation-from-substrate-2026-06-22
status: 🟢 Phase 1 COMPLETE — werdykt wyliczony (LOCK §4 wiersz 1), sympy 5/5 PASS
anti_lakatos_lock: ACTIVE
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
---

# Phase 1 — wyniki (value-blind)

## §0 — Executive summary
**Sympy 5/5 PASS. Werdykt wyliczony z reguły LOCK §4 = wiersz 1: `c₀` jest normalizacją
tego samego operatora złożonego $\partial\hat s\,\partial\hat s$ w kanale spin-2, którego
współczynnik $p^2$ (= $C_\sigma$) udowodniono (F1, #33) jako UV-czuły ⟹ `c₀` dziedziczy
status WOLNEGO PARAMETRU UV.** Dwa niezależne argumenty zbieżne (UV + proweniencja matching/kalibracja).

## §1 — Walidacja silnika (forbidden #4)
| Check | Wynik | Cel | Status |
|---|---|---|---|
| V1 | $\int_{-1}^{1}(1-\mu^2)^2(4\mu^2-1)d\mu=-16/35$ | −16/35 (F1, #33) | ✅ PASS |
| V2 | $\int_{-1}^{1}P_2(\mu)d\mu=0$ | 0 (F2, #34) | ✅ PASS |

Silnik reprodukuje **oba** kanoniczne wyniki cykli nadrzędnych przed jakimkolwiek nowym twierdzeniem.

## §2 — Argument A (UV): c₀ żyje w tym samym kanale i ma tę samą rozbieżność
- **C1 (sympy):** projekcja $\int \mu^2 P_2\,d\mu = 4/15\neq0$ ⟹ operator dwupochodny
  $\partial\hat s\partial\hat s$ **realnie zasila kanał spin-2** — tam właśnie żyje sprzężenie
  $\sigma^{ij}C(\psi)$, a `c₀ = C(ψ=1)` jest jego normalizacją (F4).
- **C3 (sympy):** modelując współczynnik $p^2$ jako $C(\Lambda)=C_{\rm fin}+\kappa_{\rm lin}\Lambda$
  (F1: $\kappa_{\rm lin}\propto-16/35\neq0$), mamy $\partial_\Lambda C=\kappa_{\rm lin}\neq0$ ⟹
  **brak scheme-independent continuum**. Normalizacja tego operatora — **niezależnie czy nazwiemy ją
  $C_\sigma$ czy $c_0$** — jest UV-czuła, chyba że ustali ją **zewnętrzna tożsamość** (Ward / residuum).
- **F2/#34:** takiej tożsamości brak — **nie ma bieguna spin-2** (∫P₂=0), więc nie ma residuum on-shell
  do ustalenia normalizacji. ⟹ kanał wiersza 2 §4 (stosunek/tożsamość anulująca `D`) **nie zademonstrowany**.

## §3 — Argument B (proweniencja): „c₀=4π" z 2026-05-09 NIE jest wyprowadzeniem
Niezależnie od UV, źródło „4π" w cyklu 2026-05-09 (F5) to:
- $\xi_{\rm eff}=4\pi G\Phi_0^2$ z OP-7 T3.4 = **thm:amplitude-matching** = **warunek dopasowania**
  (R3, zalockowany w #35 jako „matching ≠ predykcja"). **Matching nie jest derywacją normalizacji.**
- mnożnik $\xi/G\approx1.06$ = **kalibracja do GW150914** (fit do danych) — wykluczony jako dowód (LOCK §7 FALSE-POSITIVE-ii).
- **C2 (sympy):** iloczyn $4\pi\cdot\frac{1}{3\pi}=\frac{4}{3}$ jest **algebraicznie trywialny** — każdy
  $c_0$ i $\kappa_\sigma$ o iloczynie 4/3 go odtworzy; „remarkable alignment" jest **konstruowalny post-hoc**
  (sam cykl 2026-05-09 §3.3 to przyznał). Wykluczony jako dowód (LOCK §7 FALSE-POSITIVE-iii).

⟹ Również **bez** argumentu UV: cykl 2026-05-09 nie dostarczył wyprowadzenia `c₀` — dostarczył
**matching + kalibrację**, czyli dokładnie status, który #35 zasądził dla amplitudy ($\kappa_E$).

## §4 — Werdykt wyliczony (LOCK §4)
Warunek wiersza 1 spełniony: `c₀` = współczynnik $p^2$ tego samego operatora co $C_\sigma$ **i** dziedziczy `D`
(Argument A); kanał wyjścia (wiersz 2/3) **nie zademonstrowany** (brak Warda/bieguna; matching/kalibracja
wykluczone). ⟹

> **`c₀` = WOLNY PARAMETR UV** (nie predykcja). Sektor grawitacyjno-radiacyjny TGP po falsyfikacji:
> **okno „c₀·κ_σ=4/3" to PODWÓJNE strojenie** — c₀ wolne ∧ κ_E ($=C_\sigma\sigma_0^2$) wolne (#34/#35).
> **Sektor NIE dostarcza falsyfikowalnej predykcji amplitudy/fazy GW.** Niezmienione (R2): liczba modów
> (2 TT + breathing), c_GW=c, m_σ²=2m_s² — to pozostają predykcjami **strukturalnymi**, niezależnymi od c₀.

**Pewność:** wysoka. Argument A (UV) + Argument B (proweniencja) są niezależne i zbieżne. Domknięcie
do FINAL nie wymaga Phase 2 (oba kanały wyjścia §4 zamknięte: brak tożsamości Warda — F2; matching/kalibracja
wykluczone — R3/§7). Ewentualny Phase 2 mógłby jedynie *szukać* nieznanej tożsamości Warda dla σ_ab
(track wieloletni, niski priorytet, analog „ontologia σ_ab" z #34 zakończenia).

## §5 — Anti-Lakatos (checklist)
- ✓ Kryteria (LOCK §4) i lean zalockowane przed liczbami; werdykt **wyliczony**, nie wybrany.
- ✓ Silnik zwalidowany (V1, V2) przed nowym twierdzeniem.
- ✓ `c₀` **NIE sfabrykowane** do 4/3 — przeciwnie, pokazano trywialność tej zgodności (C2).
- ✓ Wynik **negatywny zgłoszony wprost**. ✓ R1/R2/R3/R4/R5 nieruszone. ✓ Budżet stałych 0.
- ✓ Rdzeń / cykl 2026-05-09 / FOUNDATIONS **NIE edytowane** — rekomendacje w Phase_FINAL §3 (czekają na autoryzację).
