---
title: "SCOPING — op-disformal-hamiltonian-viability: formalnie rozstrzygnąć werdykt sektora radiacyjnego via DISFORMAL VIABILITY (sygnatura/hiperboliczność g_eff + DOF count slaved-TT + niezależność od O12 + skaling screeningu) — zastępuje nierobustny argument induced-TT z op-disformal-stability Phase 2"
date: 2026-06-14
type: meta-scoping
status: "PRE-PHASE-0 NOTE (wymaga własnego Phase 0 + 'działaj'; zero werdyktów)"
origin: "User 2026-06-14 (sesja #26): audyt op-disformal-stability → 'C'. AUDIT_verdict wykazał: argument Phase 2 (induced-TT) nierobustny; konkluzja BROKEN prawdopodobnie poprawna via viability g_eff. Ten cykl ma to potwierdzić SOLIDNIE."
parent_cycles:
  - "[[../research/op-disformal-stability-2026-06-14/]] (Phase 1 ✓ B<0; Phase 2 argument nierobustny — NIE lockować)"
  - "[[../research/op-disformal-stability-2026-06-14/AUDIT_verdict_2026-06-14.md]] (reviewer audit: BROKEN via viability)"
resolves: "Werdykt sektora radiacyjnego (cofa UNDERDETERMINED rodzica) — SOLIDNYM argumentem viability, nie induced-TT"
anti_lakatos_note: "Liczby poprzedników LOCKED. Pre-derywacja §2 = oczekiwanie (BROKEN via viability), rachunek nadrzędny. Cykl może też dać NOT-BROKEN, jeśli viability da się ocalić."
tags:
  - scoping
  - disformal-viability
  - signature-hyperbolicity
  - DOF-count
  - O12-independence
  - fast-kill
---

# SCOPING — viability disformalna jako właściwy test sektora radiacyjnego

## §1 — Pytanie (sformułowanie)

AUDIT_verdict ([[../research/op-disformal-stability-2026-06-14/AUDIT_verdict_2026-06-14.md]])
wykazał EXACT: $g_{\rm eff}=\mathrm{diag}(-A,A+bG^2,A,A)$, $\det=-A^4(1+u)$ — wartość własna radialna
$A(1+u)$ flipuje przy $|u|=1$ ($=r_V$) dla $B<0$. Stąd trylemat: silne ekranowanie ($|u|\gg1$) jest
albo geometrycznie niedopuszczalne ($B<0$, flip sygnatury), albo ghostowe ($B>0$), a viable+zdrowy
reżim ($B<0,|u|<1$) nie ekranuje ($\Rightarrow$ PR-025 konforemny stoi).

> **Czy intersekcja {$g_{\rm eff}$ Lorentzowska} ∩ {skalar bez ghost/instability} ∩ {nietrywialne
> ekranowanie $\dot P_b$} jest PUSTA dla każdego $B(\Phi)$, znaku i magnitudy (niezależnie od O12)?
> Jeśli TAK ⇒ sektor radiacyjny BROKEN (solidnie, via viability). Jeśli istnieje okno ⇒ NOT-BROKEN
> (sign+magnitude-pinned), werdykt rodzica pozostaje UNDERDETERMINED/węższy.**

## §2 — Pre-derywacja (oczekiwanie: BROKEN via viability; rachunek nadrzędny)

Trzy kroki do potwierdzenia formalnie:
1. **Sygnatura/hiperboliczność (rdzeń):** $1+(B/A)X>0$ wymagane dla Lorentzowskości; dla $B<0$, $X>0$
   ⇒ $|u|<1$. Sprawdzić na pełnej metryce (nie tylko diagonalnej reprezentacji), w tym strefę falową
   i bliską; potwierdzić $|u|=1\equiv r_V$.
2. **DOF count (potwierdzenie audytu §2):** jawnie pokazać, że δg jest slaved do δΦ (metryka emergentna,
   eq:disformal-perturbation — wszystkie człony ∝ δΦ/∂δΦ) ⇒ induced-TT NIE jest propagującym DOF ⇒
   $c_T^2<0$ z `prop:cT` to artefakt, nie niestabilność (uzasadnia odrzucenie argumentu Phase 2).
3. **Niezależność od O12:** dla DOWOLNEGO $B$: albo $|u|<1$ (viable, brak screeningu → PR-025), albo
   $|u|>1$ (screening → patologia). Trylemat domyka się **bez** znajomości $B(\Phi)$ — to czyni werdykt
   robustnym (silniejszym niż UNDERDETERMINED, które czekało na O12).

Warunkowo (gdyby viability ocalała): skaling ekranowania $S(u)$ — czy pełny rachunek DRW potwierdza
$S\!\sim\!1/|1-u|$ (silne ekranowanie wymaga $|u|>1$), czy istnieje kanał ekranujący przy $|u|<1$.

## §3 — Falsyfikatory (szkic do LOCK)

| ID | Test | Reguła |
|---|---|---|
| F-VIA-A | Sygnatura $g_{\rm eff}$ (eigenvalues / det) w funkcji $u$, oba znaki B, strefa bliska+falowa | flip nieusuwalny dla screening-reżimu ⟹ wkład BROKEN |
| F-VIA-B | DOF count: czy induced-TT propaguje niezależnie (jest własny człon kinetyczny) | NIE (slaved) ⟹ argument Phase 2 formalnie odrzucony; TAK ⟹ re-otwiera |
| F-VIA-C | Niezależność od O12: trylemat dla dowolnego B/znaku/magnitudy | intersekcja pusta ∀B ⟹ BROKEN robustny; okno ⟹ NOT-BROKEN |
| F-VIA-D | (warunkowy) skaling $S(u)$ — czy screening wymaga $\|u\|>1$ | tak ⟹ trylemat domyka; nie ⟹ przelicz |
| F-VIA-E | Agregat z flag | (A∧C) ⟹ **BROKEN-via-viability**; okno ⟹ **NOT-BROKEN** (sign+mag-pinned) |

## §4 — Forbidden moves (szkic)

1. Liczby PR-025/survival/parent LOCKED.
2. Zakaz powrotu do induced-TT jako argumentu (audyt wykazał nierobustność) — wolno cytować jako
   „błędna ścieżka", nie jako dowód.
3. Zakaz strojenia B/M_* (O12) pod wynik; trylemat ma być O12-niezależny (F-VIA-C).
4. Zafiksować konwencję sygnatury + znak X PRZED rachunkiem (mostly-plus, statyczne X>0).
5. Bilans energii/skaling screeningu z op-disformal-radiation-resolution (heurystyka) — oznaczyć jako
   takie; pełny DRW poza zakresem (warunkowo F-VIA-D).
6. Budżet nowych stałych 0.

## §5 — Format i koszt

Fast-kill strukturalny, 1–2 fazy, sympy (det/eigenvalues g_eff; DOF count), 0 danych. Reuse:
operator $Z^{\mu\nu}$ + g_eff z AUDIT_verdict_sympy (EXACT). **Win-win:** BROKEN-via-viability ⟹
solidne domknięcie sektora grawitacyjnego (z poprawnym dowodem, zastępuje induced-TT); NOT-BROKEN ⟹
istnieje viable okno ekranujące (sensacja — wymaga sign+magnitude pin z O12).

## §6 — Rejestracja

**`op-disformal-hamiltonian-viability`** — REGISTERED (NOT activated). Cykl-rodzic op-disformal-stability
**pozostaje z wstrzymaną Phase FINAL** dopóki ten cykl nie dostarczy poprawnego argumentu; wtedy
op-disformal-stability domyka się werdyktem **via viability** (nie induced-TT).
