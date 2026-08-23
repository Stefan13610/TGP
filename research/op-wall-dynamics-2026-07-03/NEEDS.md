---
title: "NEEDS — proponowane edycje core po op-wall-dynamics (user-gated)"
date: 2026-07-03
type: needs
tgp_owner: research/op-wall-dynamics-2026-07-03
status: "N1–N3 EXECUTED 2026-07-03 (user-gate przyznany: 'NEEDS N1–N3'); N4 ROZPISANE 2026-07-03 jako cykl op-nonlinear-charge-constraint-2026-07-03 (PHASE0-LOCKED, realizacja: następna sesja)"
---

# NEEDS — edycje core wynikające z op-wall-dynamics (wszystkie addytywne, user-gated)

## N1 (HIGH) — `sek08b` `rem:ghost-artifact-scope-CP7` (uzupełnienie addytywne)

Remark deklaruje hipotezę roboczą autora (ściana = próg budżetu) i
odsyła do tego cyklu. Proponowane uzupełnienie (nowa uwaga po
istniejącej lub dopisek datowany):
1. **Wersja liniowa hipotezy obalona**: więzy budżetowe K1–K4 i pary
   nie zerują indeksu μ/τ (minimum: μ→1, τ→1; mod rezydualny nie jest
   kierunkiem rodziny profili). Hipoteza wymaga więzu
   nieliniowego/ładunkowego (typ Q-ball z symetrią, nie liniowa
   deflacja) — status: OPEN.
2. **Gładka regularyzacja EL nie zastępuje odbicia**: w rodzinie
   f_ε=½[f+√(f²+ε²)] profil τ kolabuje dla każdego ε∈{0,2…0,02}
   (jak substrat α=1); spektra μ/τ niezbieżne przy ε→0 — odbicie
   pozostaje jedynym znanym nośnikiem mechanizmu gen-3.
3. **Nowy fakt strukturalny**: górny próg g_crit=8/5 (H7) pokrywa się
   numerycznie (0,71%) z progiem pierwszego kontaktu profilu ze ścianą
   dolną g*=e^{−1/4} (g₀_wall=1,6114) — dwa progi = dwie projekcje
   jednego mechanizmu ściennego. (Cytować: research/op-wall-dynamics-2026-07-03,
   Phase3.)

## N2 (MED) — `papers_external/paper_lepton_masses` Limitations T-OP4 (doprecyzowanie)

T-OP4 (spectral stability μ/τ OPEN) — dopisać wynik tego cyklu
(zgodnie z dyspozycją LOCK-a: „jeśli negatywne → utrzymać OPEN +
doprecyzować"):
- "Linear budget-type constraints (four families + pairs, pre-registered)
  do not restore stability: the minimum attainable saddle index is 1
  for both mu and tau, and the surviving mode is not the
  profile-family direction (overlap < 0.01)."
- "The mass-ratio mechanism is sensitive to the wall model: within the
  smooth regularization family f_eps the tau profile ceases to exist
  (collapse) for every eps tested, and the mu-based r21 drifts by
  +1.9%…+23% vs the hard-wall baseline (allowed 0.1%). Stability and
  wall-model independence remain OPEN (nonlinear/charge-type
  constraint deferred)."

## N3 (LOW) — `sek08b` okolice H7/H8 (lub notatka przy `cor:ghost-artifact`)

Addytywna uwaga: g_crit=8/5 z H7 („próg kolapsu") zinterpretowany
operacyjnie: to amplituda centralna, przy której profil PO RAZ PIERWSZY
dotyka ściany kinetycznej g* (zmierzone g₀_wall=1,6114; guard g*+0,005;
zgodność 0,71%); konieczny warunek analityczny (bez tarcia):
W(g₀) ≤ W(g*) ⇒ g₀ ≥ 1,1696 (tarcie 2/r przesuwa próg wyżej).
Wzmacnia H7/H8 strukturalnie; nie zmienia żadnych wartości.

## N4 (RESEARCH, nie-core) — następny cykl (propozycja, NIE zobowiązanie)

Hipoteza budżetowa po tym cyklu wymaga wersji NIELINIOWEJ:
(a) właściwy więz Q-ball: ładunek zachowany z symetrii substratu
    (jaka symetria niesie „budżet tworzonej przestrzeni"? — pytanie
    pierwszoprincypialne, nie numeryczne);
(b) stabilność orbitalna na poziomie nieliniowym (energia przy stałym
    ładunku, kryterium Vakhitova–Kolokolova zamiast indeksu siodłowego);
(c) model ściany z dyskretności substratu (poza rodziną f_ε);
(d) obserwacja K4 (mody głębokie ⊂ podprzestrzeń budżetu rdzenia)
    jako wskazówka co do NOŚNIKA właściwego ładunku.

---

## Log

- 2026-07-03 (sesja #62, agent-implementator): utworzone po zamknięciu
  Phase 1–3. Zero edycji core wykonanych przed user-gate.
- 2026-07-03 (sesja #62, user-gate przyznany: „NEEDS N1–N3"):
  - N1+N3 ✅ `sek08b` — NOWY `rem:wall-dynamics-2026-07-03` po
    `rem:ghost-artifact-scope-CP7` (3 punkty: W1 liniowy negatywny +
    fakt K4; W2 kolaps τ w f_ε + dryf r₂₁; W3a g_crit=8/5 ⟺ próg
    aktywacji ściany, warunek konieczny g₀≥1,1696) + pointer w
    `rem:spectral-CP7` pkt 3 („wykonany 2026-07-03…").
  - N2 ✅ lepton paper — Limitations T-OP4 doprecyzowane (linear
    constraints insufficient; wall-model sensitivity; hard-wall
    reflection = structural input; OPEN na poziomie nieliniowym).
  - Build-gate'y: `main.tex` exit 0, 7 undefined refs = identyczny
    zbiór pre-existing (#32), **0 nowych**, 0 błędów;
    `tgp_lepton_masses.tex` exit 0, 0 undefined refs, 0 błędów.
  - N4: pozostaje jako propozycja (decyzja użytkownika).
- 2026-07-03 (sesja #62, autoryzacja: „rozpisz cykl badawczy N4 dla
  nowego agenta"): N4 → **ROZPISANE** jako
  [[../op-nonlinear-charge-constraint-2026-07-03/Phase0_balance.md]]
  (PHASE0-LOCKED; realizacja: następna sesja, osobny agent).
