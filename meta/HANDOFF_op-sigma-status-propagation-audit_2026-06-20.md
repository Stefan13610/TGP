---
title: "HANDOFF / PROMPT — op-sigma-status-propagation-audit (nowa sesja): audyt spójności propagacji statusu κ_E = wolny parametr po sesjach #33/#34"
date: 2026-06-20
type: handoff-prompt
status: REGISTERED-QUEUED (oczekuje „działaj")
parent_cycles:
  - "[[../research/op-CG4-substrate-closure-2026-06-20/Phase_FINAL_close.md]]"
  - "[[../research/op-sigma-ab-pole-residue-2026-06-20/Phase_FINAL_close.md]]"
coordination: "[[../STATE.md]] (single-source)"
---

# PROMPT dla nowego agenta — op-sigma-status-propagation-audit

> Skopiuj treść poniższego bloku jako prompt nowej sesji (agent startuje bez kontekstu).

---

## Prompt (do wklejenia)

```
Jesteś ekspertem fizyki teoretycznej pracującym w projekcie TGP_v1 (Obsidian vault).
Zadanie: audyt spójności propagacji statusu po sesjach #33/#34.

KONTEKST (przeczytaj najpierw [[STATE.md]] sesje #33, #34):
W sesjach #33 (op-CG4-substrate-closure) i #34 (op-sigma-ab-pole-residue) dowiedziono,
że C_σ (≡ κ_E = 8πG₀C_σσ₀²/c³) — stała kinetyczna tensorowego nośnika GW σ_ab — jest
NIEREDUKOWALNYM WOLNYM PARAMETREM, nie predykcją. Dwa niezależne dowody:
 (a) wsp. p² operatora złożonego ∂ŝ∂ŝ w kanale TT spin-2 ma niezerową rozbieżność
     liniową UV (wsp. kątowy −16/35; analit.=num do 4 cyfr) → brak continuum scheme-indep.;
 (b) brak izolowanego bieguna spin-2 — kontakt φ⁴ ma ∫P₂dx=0 (s-wave nie wiąże d-wave)
     → brak residuum on-shell; „M²=2m_s²" = coeff OPE, nie masa bieguna.
Status zapisano w rdzeniu: sek08 rem:sigma-params + nowa rem:sigma-Csigma-free;
dodatekQ Q.5 (CG-4 status); PREDICTIONS_REGISTRY wiersz GW7 (FREE-PARAMETER).
Substrat NIE jest problemem (M0 φ⁴ Isinga = niepatologiczny; patologia była w bondzie M1).

CEL AUDYTU (value-blind, anti-Lakatos):
Znaleźć WSZYSTKIE miejsca w rdzeniu (.tex) i dokumentach głównych, gdzie sektor
tensorowy GW / σ_ab / κ_E / C_σ jest wciąż opisywany jako „przewidywany / wyprowadzony /
wyznaczony / derived / predicted / determined" w sposób NIESPÓJNY z nowym statusem
(κ_E = wolny parametr). Sklasyfikować każde trafienie: SPÓJNE / DO-POPRAWY / NIEJASNE.

ZAKRES (konkretne kotwice do sprawdzenia — NIE wyczerpujące):
 1. README.md / abstract main: framing „40 predictions from 3 inputs" oraz LICZBA
    PARAMETRÓW. Czy „3 inputs" (g₀, Ω_Λ, N) myli się z „3 parametry fenomenologiczne"
    (2 skalar + 1 tensor wolny, rem:sigma-params)? Sprawdzić spójność narracji.
 2. sek08: thm:amplitude-matching (ξ_eff matching), ssec:tensor-substrate, tabela
    rem:sigma-logic (mapa statusu sektora tensorowego) — czy status σ_ab/amplitudy
    spójny z rem:sigma-Csigma-free?
 3. GW4 (m_σ²/m_s²=2, registry, LOCKED): to OPE-invariant, ZOSTAJE LOCKED — UWAGA:
    NIE mylić go z C_σ. Potwierdzić, że audyt go nie „downgraduje" błędnie (m_σ² to
    masa, nie stała kinetyczna). Rozróżnienie jawne.
 4. TGP_FOUNDATIONS.md: status sektora radiacyjnego (CL-7 UNDERDETERMINED z sesji #29)
    — zaktualizować do „κ_E = wolny parametr (dowiedzione #33/#34)" jeśli niespójne.
 5. tgp_letter.tex / tgp_companion.tex / papers_external/*: czy gdzieś GW/σ_ab
    reklamowane jako predykcja bez warunkowości na κ_E?
 6. PREDICTIONS_REGISTRY: czy GW2/GW3/GW5/GW6 (polaryzacje, h_b=h_L, no-vector,
    masslessness) są niezależne od κ_E (powinny być — to liczba modów/symetria, nie
    amplituda) — potwierdzić, że audyt ich nie rusza.

REGUŁY (anti-Lakatos, OBOWIĄZKOWE):
 - Phase 0 LOCK przed jakąkolwiek edycją: lista trafień + klasyfikacja, zalockowana.
 - Phase0_balance.md (gate) jeśli dotyka registry/statusów.
 - Edycje rdzenia TYLKO addytywne/korekcyjne, z autoryzacją; po edycji core:
   pdflatex main.tex → exit 0, nowe cross-refy rozwiązane (build verification).
 - Rozróżnić: masa σ (m_σ²=2m_s², LOCKED, OPE) vs stała kinetyczna C_σ (wolny parametr).
   NIE downgradować GW4. NIE fabrykować.
 - Budżet nowych stałych 0. Higiena: usunąć __pycache__ / zagnieżdżone błędne katalogi.
 - Wynik negatywny („wszystko już spójne") jest dopuszczalny i wartościowy.

DELIVERABLE:
 - Cykl research/op-sigma-status-propagation-audit-2026-06-20/ (README + Phase0_lock +
   Phase1_audit + Phase_FINAL_close) z tabelą trafień i klasyfikacją.
 - Lista poprawek rdzenia (jeśli są) — zastosować po autoryzacji, z build verification.
 - Aktualizacja STATE.md (sesja #35).
```

---

## Notatka koordynacyjna (dla autora, nie część promptu)
- Łuk #33/#34 zamknięty; rdzeń zaktualizowany (build 553 str., exit 0).
- Ten audyt = domknięcie spójnościowe (propagacja statusu), niska–średnia złożoność, 1 sesja.
- Po audycie: rozważyć opcję (a)-bis (ontologia σ_ab / pole-residue) jako wieloletni research-track,
  LUB zostawić κ_E jako jawny wolny parametr (opcja b, już zapisana).
