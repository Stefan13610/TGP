---
title: "Phase 0 LOCK — op-sigma-status-propagation-audit"
date: 2026-06-22
type: phase_lock
status: 🔒 LOCKED (przed jakąkolwiek edycją rdzenia)
phase: 0
cycle: op-sigma-status-propagation-audit-2026-06-20
anti_lakatos_lock: ACTIVE
---

# Phase 0 — LOCK (kryteria zalockowane przed edycją)

> Reguła nadrzędna: kryteria klasyfikacji i zakres są ustalone **przed** zastosowaniem
> jakiejkolwiek poprawki rdzenia. Werdykt każdego trafienia jest **wyliczany** z reguły,
> nie wybierany. Wynik negatywny („wszystko już spójne") dopuszczalny i wartościowy.

## §1 — Fakt referencyjny (zalockowany, z #33/#34)
`F`: **κ_E ≡ 8πG₀C_σσ₀²/c³ (≡ stała kinetyczna C_σ) = nieredukowalny WOLNY PARAMETR UV.**
- Dowód (a): rozbieżność liniowa UV wsp. p² operatora ∂ŝ∂ŝ w kanale TT spin-2 (wsp. −16/35 ≠ 0).
- Dowód (b): brak izolowanego bieguna spin-2 (∫P₂dx=0; s-wave nie wiąże d-wave) → brak residuum.
- Survival (κ_E=5/6) osiągalne wyłącznie przez **wybór** parametru; M911-* warunkowe na κ_E.
- Domknięcie sektora = **opcja (b)** (przyjęcie κ_E jako wolnego parametru, uczciwe param-counting).

## §2 — Rozróżnienia zalockowane (NIE downgradować)
- `R1`: **m_σ² = 2 m_s²** (GW4) = masa / coeff OPE-invariant, **LOCKED**. NIE jest C_σ. NIE rusza audyt.
- `R2`: **GW2/GW3/GW5/GW6** (3 DOF polaryzacyjne, h_b=h_L, brak modów wektorowych, bezmasowość)
  = liczba modów / symetria — **niezależne od κ_E**. NIE rusza audyt.
- `R3`: **thm:amplitude-matching** (ξ_eff = 4πG₀σ₀Φ₀) = **warunek dopasowania**, nie predykcja.
  Spójny z F (matching ≠ pinowanie C_σ). NIE jest overclaimem sam w sobie.
- `R4`: **c_GW = c₀** (GW1) = własność tła geometrycznego — niezależne od κ_E.

## §3 — Reguła klasyfikacji (wyliczana)
Dla każdego trafienia T opisującego sektor tensorowy GW / σ_ab / κ_E / C_σ:
- **SPÓJNE**: T opisuje (i) mechanizm/strukturę/liczbę modów/symetrię/matching jako hipotezę,
  propozycję lub warunek; LUB (ii) explicite κ_E/C_σ jako wolny parametr; LUB (iii) wielkość
  z R1–R4. Brak konfliktu z F.
- **DO-POPRAWY**: T twierdzi, że sektor tensorowy (amplituda / C_σ / κ_E / „sektor GW") jest
  „przewidziany / wyprowadzony / wyznaczony / zamknięty / pinowalny z substratu" **bez**
  warunkowości na κ_E; LUB globalna liczba parametrów TGP jest podana jako „2" (pomija
  nieredukowalny C_σ) bez zakresowania do sektora skalarnego; LUB status LIVE sektora
  radiacyjnego mówi „domknięcie = pinowanie κ_E" (już rozstrzygnięte negatywnie).
- **NIEJASNE**: T jest poprawne w swoim (wąskim) zakresie, ale brak jawnego krzyżowego
  odniesienia do rem:sigma-Csigma-free może wprowadzać w błąd value-blind czytelnika.
  Rekomendacja: lekki cross-ref; niski priorytet.

## §4 — Zakres (kotwice; NIE wyczerpujące)
1. README.md / abstract — framing „40 predictions / 3 inputs" + liczba parametrów.
2. sek08 — thm:amplitude-matching, ssec:tensor-substrate, rem:sigma-params, rem:sigma-logic,
   rem:parsimony, rem:sigma-Csigma-free (referencja docelowa spójności).
3. GW4 (m_σ²/m_s²=2) — potwierdzić, że NIE zostaje downgradowany (R1).
4. TGP_FOUNDATIONS.md — CL-7 / status sektora radiacyjnego.
5. tgp_letter.tex / tgp_companion.tex / papers_external/* — GW/σ_ab jako predykcja bez warunkowości?
6. PREDICTIONS_REGISTRY — GW2/3/5/6 niezależne od κ_E (R2); GW7 = FREE-PARAMETER.
7. core/_meta_latex/status_map.tex, tabela_epistemiczna.tex — etykiety statusu σ_ab + liczba param.
8. sek07_predykcje.tex, dodatekQ — status sektora GW + CG-4.

## §5 — Forbidden moves (anti-Lakatos)
1. NIE downgradować GW4 (R1) ani GW2/3/5/6 (R2).
2. NIE fabrykować — żadnego „κ_E przewidziany" nie wprowadzać.
3. Edycje rdzenia TYLKO addytywne/korekcyjne, z autoryzacją usera, po Phase0_balance.
4. Budżet nowych stałych = 0.
5. Po edycji core: `pdflatex main.tex` → exit 0 + nowe cross-refy rozwiązane (build verification).
6. Wynik negatywny zgłosić wprost, jeśli zachodzi.
7. NIE mylić zakresu: „2 parametry skalarne (Warstwa II)" ≠ „2 parametry TGP" (z C_σ = 3).

## §6 — Falsyfikator audytu (dwustronny)
- Audyt FALSE-NEGATIVE jeśli pominie trafienie spełniające regułę DO-POPRAWY (§3).
- Audyt FALSE-POSITIVE jeśli zaklasyfikuje jako DO-POPRAWY trafienie chronione R1–R4 (§2).
