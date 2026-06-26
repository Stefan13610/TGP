---
title: "Phase 0 LOCK — op-alpha2-status-propagation-audit"
date: 2026-06-22
type: phase_lock
status: 🔒 LOCKED (przed jakąkolwiek edycją rdzenia / papers)
phase: 0
cycle: op-alpha2-status-propagation-audit-2026-06-22
anti_lakatos_lock: ACTIVE
parent_cycles:
  - "[[../op-A3-alpha-resolution-2026-06-14/Phase_FINAL_close.md]]"
  - "[[../op-amplitude-density-global-audit-2026-06-16/Phase_FINAL_close.md]]"
template: "[[../op-sigma-status-propagation-audit-2026-06-20/Phase0_lock.md]]"
---

# Phase 0 — LOCK (kryteria zalockowane przed edycją)

> Reguła nadrzędna: kryteria klasyfikacji i zakres są ustalone **przed** zastosowaniem
> jakiejkolwiek poprawki. Werdykt każdego trafienia jest **wyliczany** z reguły, nie
> wybierany. Wynik negatywny („wszystko już spójne") dopuszczalny i wartościowy.

## §1 — Fakt referencyjny (zalockowany, z #31 op-A3 + #32 audit gęstości)
`F`: **α = 2 NIE jest wyprowadzone z substratu pod ontologią Φ = ⟨ŝ²⟩ (gęstość).**
- Dowód (op-A3, sympy 5/5, value-blind): relacja EL α = p/2 dla K = φ^p jest **poprawna**;
  ale substrat K(u) ∝ u⁴ w polu **mikro** u + Φ = ⟨ŝ²⟩ ∝ u² daje po poprawnej zmianie
  zmiennych K_eff(ψ) ∝ ψ (liniowy) → **α = ½**, nie 2. α = 2 pojawia się w thm:D-uniqueness C3
  tylko przez **pominięcie** transformacji Φ ∝ u² (traktuje makro φ = Φ/Φ₀ jak mikro u).
- Niespójność potwierdzona z DWÓCH miejsc rdzenia (thm:D-uniqueness vs sek10 eq:kinetic_macro).
- **Status zasądzony (#32, w rdzeniu):** Φ = ⟨ŝ²⟩ pozostaje polem kanonicznym; **α = 2 ma status
  *selekcji aksjomatycznej na gęstości* (C1–C3), NIE derywacji** — sek08 `rem:alpha2-pivot-status-pl`,
  sek10:208–216. Domknięcie = **opcja C** (α = 2 jako postulat fenomenologiczny; thm:D-uniqueness =
  argument spójnościowy/selekcja w klasie, nie wyprowadzenie z substratu).

## §2 — Rozróżnienia zalockowane (NIE downgradować, NIE nadinterpretować)
- `R1`: **Relacja EL α = p/2 (sek08 eq:K-ode)** = **POPRAWNA**, niepodważona. Audyt NIE rusza formalizmu EL.
- `R2`: **α = 2 pozostaje fenomenologicznie wiarygodne** (PPN γ=β=1, profil solitonu, widmo mas,
  Koide). Audyt **NIE** twierdzi „α ≠ 2". Audyt dotyczy WYŁĄCZNIE statusu *epistemicznego*
  (selekcja/postulat vs derywacja).
- `R3`: **thm:D-uniqueness / thm:alpha2 jako selekcja w klasie C1–C3** = ważne twierdzenie
  o jednoznaczności *w obrębie założonej klasy*. Overclaimem jest tylko ramowanie
  „wyprowadzone / proved from the substrate Hamiltonian" bez kwalifikatora selekcji.
- `R4`: **Φ = ⟨ŝ²⟩ (gęstość) = pole kanoniczne** (zasądzone #32). Audyt **NIE** reaktywuje
  Opcji B (amplituda jako pole kanoniczne) — to byłby forbidden drift.

## §3 — Reguła klasyfikacji (wyliczana)
Dla każdego trafienia T opisującego status α = 2 / K(φ) = φ⁴ / thm:alpha2 / D-uniqueness:
- **SPÓJNE**: T opisuje α = 2 jako (i) *selekcję* / *postulat na gęstości* / aksjomat
  fenomenologiczny / wymóg spójności (PPN, masy); LUB (ii) jednoznaczność *w klasie C1–C3*
  bez claimu derywacji z substratu; LUB (iii) explicite cytuje rem:alpha2-pivot-status / op-A3.
  Brak konfliktu z F.
- **DO-POPRAWY**: T twierdzi, że α = 2 (lub K = φ⁴ na gęstości) jest **„wyprowadzone /
  derived / proved / algebraic theorem / fixed from the substrate Hamiltonian"** — bez
  kwalifikatora selekcji/postulatu, sprzecznie z F; LUB liczy „3 inputy" prezentując α = 2
  jako *wynik* (a nie niezależny postulat) w sposób, który ukrywa, że selekcja C1–C3 jest
  wejściem aksjomatycznym.
- **NIEJASNE**: T jest poprawne w wąskim zakresie (np. „α = 2 algebraic theorem" rozumiane
  jako jednoznaczność w klasie), ale brak jawnego cross-ref do rem:alpha2-pivot-status może
  wprowadzić value-blind czytelnika w błąd „derywacja z substratu". Rekomendacja: lekki cross-ref.

## §4 — Zakres (kotwice; NIE wyczerpujące)
1. README.md — „Kinetic coupling α = 2: algebraic theorem, not a fit"; „established through
   prop:substrate-action + classification theorem thm:alpha2"; nagłówek „40 z 3 inputów".
2. tgp_letter.tex — „α = 2 is derived as an algebraic theorem".
3. tgp_companion.tex — „α = 2 is proved as an algebraic theorem from [substrate]"; „conformal
   argument fixes K = g^{2α} with α = 2".
4. papers_external/tgp_core_paper/tgp_core.tex — „α = 2 as an algebraic theorem"; thm:alpha2
   „Uniqueness … α = 2"; „(C3) fixes α = 2 through the [substrate]"; hedge l.139–141 (MC/FRG fixed-point).
5. core/_meta_latex/status_map.tex — „K = K_geo φ⁴ ⇒ K(0)=0; wyprowadzenie z substratu".
6. core/formalizm/dodatekH_lancuch_wyprowadzen.tex — „Status N0-4: [AX]→[AN] (wyprowadzony z A2)";
   prop:K0-from-substrate.
7. sek00_summary.tex / tabela_epistemiczna.tex — etykieta statusu α = 2 ([AN]/[AX]) + liczba inputów.
8. Referencja docelowa spójności (SPÓJNE, NIE ruszać): sek08 rem:alpha2-pivot-status-pl,
   sek10:208–216, TGP_FOUNDATIONS.md („α = 2 selection").

## §5 — Forbidden moves (anti-Lakatos)
1. NIE podważać relacji EL (R1) ani fenomenologicznej wiarygodności α = 2 (R2).
2. NIE reaktywować Opcji B (amplituda = pole kanoniczne) — Φ = ⟨ŝ²⟩ zasądzone (R4).
3. NIE fabrykować derywacji — żadnego „α = 2 jednak wyprowadzone" nie wprowadzać.
4. Edycje rdzenia/papers TYLKO addytywne/korekcyjne, z autoryzacją usera, po Phase0_balance.
5. Po edycji core: `pdflatex main.tex` → exit 0 + nowe cross-refy rozwiązane (build verification).
6. Wynik negatywny zgłosić wprost, jeśli zachodzi.
7. Budżet nowych stałych = 0. Rozróżnić: „selekcja w klasie C1–C3" ≠ „derywacja z substratu".

## §6 — Falsyfikator audytu (dwustronny)
- Audyt FALSE-NEGATIVE jeśli pominie trafienie spełniające regułę DO-POPRAWY (§3).
- Audyt FALSE-POSITIVE jeśli zaklasyfikuje jako DO-POPRAWY trafienie chronione R1–R4 (§2)
  — np. potraktuje poprawne „jednoznaczność w klasie C1–C3" jako overclaim derywacji.

## §7 — Pytanie wtórne (parameter-counting; zalockowane jako OSOBNY deliverable)
Jeśli α = 2 jest selekcją aksjomatyczną (nie derywacją), to **czy należy ją liczyć jako
niezależne wejście** obok (g₀ᵉ, Ω_Λ, N=3)? LOCK: to pytanie jest *zgłaszane*, nie rozstrzygane
w tym cyklu jednostronnie — wymaga jawnej decyzji rdzenia (czy C1–C3 to „darmowe" założenia
strukturalne czy input). Audyt tylko **ujawnia** zależność nagłówka „3 inputy" od tej decyzji.
