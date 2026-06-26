---
title: "op-amplitude-density-global-audit — Phase 0 LOCK"
date: 2026-06-16
type: research-cycle
status: LOCKED
cycle: op-amplitude-density-global-audit
phase: 0
---

# Phase 0 — LOCK (protokół value-blind)

## Cel cyklu

Po integracji **Opcji B** (sesja #31): rdzeń wprowadził rozróżnienie
**amplituda** $\varphi$ (kanoniczne pole kinetyczne, $K=\varphi^4$, $\alpha=2$)
vs **gęstość** $\Phi=\langle\hat s^2\rangle\propto\varphi^2$ (faza metryczna, osobna).
Edycje zweryfikowano tylko pod kątem **kompilacji** (pdflatex 0 błędów),
NIE pod kątem **spójności fizyki w dół łańcucha**.

**Pytanie cyklu (value-blind):** czy każde użycie pola ($\Phi$ / $\varphi$ / $\langle\hat s^2\rangle$)
w sekcjach downstream jest spójne z konwencją LOCKED Opcji B?
W szczególności: czy gdzieś równanie pola / człon kinetyczny / metryka traktuje
**gęstość** $\Phi$ jako kanoniczne pole $\alpha=2$ (co byłoby błędem — w gęstości $\alpha=\tfrac12$)?

To jest **audyt odpowiedzialnościowy**: domknięcie luki, którą sam wprowadziłem edycją rdzenia.

## Konwencja LOCKED (z sek08 `rem:amplitude-vs-density-alpha`)

| Wielkość | Symbol | Rola | $K$ | $\alpha$ |
|---|---|---|---|---|
| **Amplituda** | $\varphi=\hat s/s_0$ (cg) | **kanoniczne pole kinetyczne** | $K(\varphi)=K_{\rm geo}\varphi^4$ | **2** |
| **Gęstość** | $\Phi=\langle\hat s^2\rangle\propto\varphi^2$ | faza metryczna (osobna) | $K_{\rm eff}(\psi)\propto\psi$ | $\tfrac12$ |

**Klucz:** $\alpha$ NIE jest niezmiennikiem redefinicji $\varphi\to\Phi\propto\varphi^2$.
Kanoniczne równanie pola TGP i $\alpha=2$ definiuje się **w amplitudzie** $\varphi$.
$\Phi$ i $\varphi$ to ta sama fizyka w różnych zmiennych — **nie** dwa różne pola fizyczne.

## Kotwice rdzenia (Opcja B, sesja #31)
- sek08 `rem:amplitude-vs-density-alpha` (l.1076) — kotwica rozstrzygnięcia
- sek10 `rem:canonical_g4` (l.179) + `rem:K_to_f_amplitude` (l.203)
- dodatekQ2 `rem:A3-correction-alpha` (l.343)
- (istniejące) sek08 `rem:alpha2-pivot-status-pl` (l.1056); paper `rem:alpha2-pivot-status`

## Zakres audytu (sekcje downstream)
1. **sek08c** metryka z substratu — $g_{\rm eff}(\Phi)$: czy faza metryczna = gęstość (poprawnie) czy myli z amplitudą?
2. **sek07 / predykcje** (cień BH, PPN $\gamma=\beta=1$) — w jakiej zmiennej wyprowadzone równanie pola?
3. **soliton / atom_from_soliton** — równanie solitonu $\nabla^2\varphi+2(\nabla\varphi)^2/\varphi$: amplituda?
4. **widmo mas / mass_scaling_k4** — skalowanie $K\propto\varphi^4$ (amplituda): spójne?
5. **sek09 / równanie pola makro** — eq:kinetic_macro (flagowane w sek10 §K_to_f jako density-frame).
6. **most_gamma_phi (dodatekQ/Q2)** — most $\gamma\leftrightarrow\Phi$: zmienna.

## Gate'y value-blind (LOCKED — kryteria orzekane PRZED audytem)

Dla każdego znalezionego użycia pola w kontekście **kinetycznym / równania pola / metryki**:

- **G-CONSISTENT**: użycie jawnie zgodne z konwencją —
  (a) kontekst kinetyczny/$\alpha=2$ → zmienna = amplituda $\varphi$ z $K=\varphi^4$; LUB
  (b) kontekst metryczny (faza, $g_{\rm eff}$) → zmienna = gęstość $\Phi=\langle\hat s^2\rangle$; LUB
  (c) konwersja $\Phi\leftrightarrow\varphi$ jawnie wykonana ($\Phi\propto\varphi^2$, eq:Phi_phi_id).
- **G-INCONSISTENT**: kontekst kinetyczny/$\alpha=2$ używa **gęstości** $\Phi$ jako pola kanonicznego
  BEZ konwersji (tj. implikuje $\alpha=2$ w gęstości — sprzeczne, bo tam $\alpha=\tfrac12$);
  LUB metryka traktuje amplitudę jako fazę bez konwersji.
- **G-AMBIGUOUS**: symbol $\Phi$/$\varphi$ użyty bez jawnej deklaracji ról, gdzie obie interpretacje
  dają różną fizykę (wymaga doprecyzowania, ale nie jest twardym błędem).
- **G-NEUTRAL**: użycie poza zakresem rozróżnienia (np. $\Phi$ jako ogólny symbol w kontekście
  termodynamicznym/statystycznym, gdzie amplituda vs gęstość nie wpływa na wynik).

## Reguła werdyktu (LOCKED — wyliczana, nie wybierana)

- **DERIVED-CONSISTENT**: 0× G-INCONSISTENT ORAZ wszystkie użycia kinetyczne/metryczne sklasyfikowane
  jako G-CONSISTENT (G-AMBIGUOUS dopuszczalne jako rekomendacje doprecyzowania, nie blokery).
  ⟹ Opcja B globalnie spójna; edycja rdzenia kompletna.
- **PARTIAL**: 0× G-INCONSISTENT, ale ≥1 G-AMBIGUOUS w miejscu wpływającym na wynik fizyczny
  (PPN/masy/metryka) ⟹ lista doprecyzowań, bez sprzeczności.
- **INCONSISTENT**: ≥1× G-INCONSISTENT ⟹ konkretna lista miejsc do poprawy (rdzeń);
  poprawki wymagają osobnej autoryzacji (chyba że trywialne notacyjne, jak w Opcji B).

## Anti-Lakatos (LOCKED)
- ✓ Werdykt **wyliczany** z reguły powyżej, nie wybierany pod z góry założony wynik „spójne".
- ✓ Każde G-INCONSISTENT **zgłoszone jawnie**, nawet jeśli oznacza, że moja edycja #31 była niepełna.
- ✓ NIE „naprawiam" przez przemianowanie symboli tak, by ukryć rozbieżność fizyczną.
- ✓ Rdzeń NIE edytowany w tym cyklu bez osobnej autoryzacji (poza zgłoszeniem listy).
- ✓ Rozróżnienie testu: NIE proponuję „obserwacyjnego testu $\alpha=2$ vs $\tfrac12$" — to ta sama
  fizyka w różnych zmiennych (źle postawione); audyt dotyczy **wewnętrznej spójności notacyjnej**.

## Plan faz
- **Phase 1** — inwentaryzacja: zmapować wszystkie użycia $\Phi/\varphi/\langle\hat s^2\rangle$
  w kontekście kinetycznym/metrycznym w sekcjach downstream (Explore, równolegle).
- **Phase 2** — klasyfikacja value-blind każdego użycia wg gate'ów; tabela.
- **Phase FINAL** — agregat → werdykt; lista G-INCONSISTENT/AMBIGUOUS; aktualizacja STATE.md.

**STATUS: LOCKED 2026-06-16.** Kryteria zamrożone przed audytem.
