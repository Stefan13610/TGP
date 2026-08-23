---
title: "NEEDS — proponowane edycje core po op-nonlinear-charge-constraint (user-gated)"
date: 2026-07-04
type: needs
tgp_owner: research/op-nonlinear-charge-constraint-2026-07-03
status: "OPEN — zero edycji core wykonanych; wszystkie propozycje addytywne, user-gated"
---

# NEEDS — edycje core wynikające z op-nonlinear-charge-constraint (user-gated)

Kontekst: cykl przetestował NIELINIOWĄ/ŁADUNKOWĄ wersję hipotezy
budżetowej autora (po obaleniu wersji liniowej w #62). Werdykty:
**V1 NEGATYWNY dla M0** (żaden budżetopodobny ładunek nie jest
zachowany), **V2 NEGATYWNY dla μ/τ w M1** (VK slope-positive na całych
gałęziach + kontinuum tachioniczne przy każdym ω + próżnia ω duchowiona
od ω≈0,2935), **V3 kierunek (i)** (niestabilność μ potwierdzona
nieliniowo w M0-f_ε: runaway do g→0 w skończonym czasie). Szczegóły i
liczby: [[README.md]].

## N1 (HIGH) — `sek08b` dopisek do `rem:wall-dynamics-2026-07-03` (addytywny, datowany)

Remark kończy się statusem „hipoteza wymaga więzu
nieliniowego/ładunkowego (typ Q-ball) — OPEN". Proponowany dopisek
(2026-07-04): wersja ładunkowa przetestowana i NEGATYWNA w pełnym
zakresie zamkniętej listy:
1. **M0 nie posiada budżetopodobnego ładunku zachowanego** — pełny
   inwentarz C1–C5 (sympy exact, test operatorem Eulera): 0/9
   zachowanych (energia = kontrola dodatnia, nie budżet).
2. **Rozszerzenie U(1) (M1, model-extension) nie stabilizuje:** ładunek
   Noether Q istnieje (exact), redukcja Q-ball φ_ω z
   W_eff = W − (ω²/2)f φ² działa i gałęzie istnieją dla ω ≤ 0,25
   (dalej kolaps na ścianę), ALE (i) VK: dQ/dω > 0 (slope-positive =
   gałąź niestabilna) dla μ/τ wszędzie; (ii) deflacja ładunkowa nie
   usuwa modów głębokich L₊; (iii) krawędź kontinuum
   c(ω) = −1 − 7ω² + O(ω⁴) OBNIŻA się z ω (nie podnosi), a próżnia
   ω-zależna φ_∞(ω) = 1 − 3ω² + … przecina ścianę duchową g* przy
   ω_gh = 0,2935 — stabilizacja ładunkowa w M1 nierozstrzygalna z
   powodu tła, a tam gdzie rozstrzygalna: negatywna.
3. **Nieliniowa ewolucja (M0-f_{ε=0,2}, jedyny N-zbieżny punkt):**
   perturbacja modu głębokiego μ (a = ±0,01; ±0,03) rośnie wykładniczo
   (σ zgodne z teorią liniową) BEZ saturacji i pole opuszcza dziedzinę
   modelu (g→0) w czasie t* ≈ 3,6 — dynamiczny odpowiednik statycznego
   kolapsu τ z W2b: miękka ściana nie trzyma pola.
   Status hipotezy budżetowej: OPEN wyłącznie POZA przetestowaną klasą
   (dyskretność substratu; inna symetria niż U(1) na fazie; sektor
   sprzężony F-A). Cytować: research/op-nonlinear-charge-constraint-2026-07-03.

## N2 (MED) — lepton paper, Limitations T-OP4 (dopisek)

Po zdaniu o „nonlinear/charge-type constraint deferred" dopisać wynik:
- "The charge-type route has now been tested and fails: the canonical
  embedding admits no conserved budget-like charge (closed inventory,
  exact); in the minimal U(1) complexification the Q-ball branches are
  Vakhitov–Kolokolov slope-positive (unstable) for mu/tau over the
  entire existence window (omega <= 0.25), the tachyonic continuum edge
  moves DOWN with omega (c = -1 - 7 omega^2 + ...), and the
  omega-vacuum crosses the ghost wall at omega ~ 0.29. Nonlinear
  evolution in the smooth-wall model confirms the instability
  dynamically (finite-time exit of the field from the model domain,
  g -> 0). Stability of mu/tau, if restorable at all, requires
  mechanisms outside the tested class (substrate discreteness, a
  different symmetry, or the gravitationally coupled F-A sector)."

## N3 (LOW) — `sek08b` przy `rem:spectral-CP7` (pointer)

Jednozdaniowy pointer: łańcuch CP-7 → op-wall-dynamics →
op-nonlinear-charge-constraint zamyka NEGATYWNIE wszystkie trzy
przetestowane ścieżki stabilizacji μ/τ (deflacja liniowa; więz
ładunkowy/VK; stabilizacja nieliniowa) — stabilność korony pozostaje
OPEN i wymaga mechanizmu spoza klasy pól gładkich z odbiciem ad-hoc.

## N4 (RESEARCH, nie-core) — możliwe następne cykle (propozycja, NIE zobowiązanie)

Wszystkie trzy ścieżki wewnątrz gładkiej teorii pola zamknięte
negatywnie. Pozostałe kierunki pierwszoprincypialne:
(a) **ściana z dyskretności substratu** — poza rodziną f_ε; wymaga
    modelu mikroskopowego (sieć/automat), nie regularyzacji;
(b) **inna symetria nośna budżetu** — nie U(1) na fazie zanurzenia,
    lecz symetria substratowa (np. objętościowa/permutacyjna); pytanie
    czysto strukturalne, przed jakąkolwiek numeryką;
(c) **sektor sprzężony F-A** — stabilność korony jako własność pełnego
    układu z grawitacją (dualizm L04): mody głębokie mogą być
    artefaktem traktowania korony w izolacji;
(d) **reinterpretacja termodynamiczna** — metastabilność + czas życia
    (tunelowanie/rozpad) zamiast stabilności absolutnej; wymaga
    oszacowania czasu życia μ/τ względem skal obserwacyjnych.

---

## Log

- 2026-07-04 (sesja #63, agent-implementator): utworzone po zamknięciu
  Phase 1–3. Zero edycji core wykonanych przed user-gate.
