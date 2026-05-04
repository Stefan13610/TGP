---
title: "POST_ACTION_UPDATE — L02 β/γ renotation EXECUTED 2026-05-04"
date: 2026-05-04
parent: "[[README.md]]"
type: audit-update
tgp_owner: audyt/L02_beta_gamma_semantics
tags:
  - audit-update
  - L02
  - executed
  - notation-fix
related:
  - "[[README.md]]"
  - "[[../../axioms/notacja/dodatekA_notacja.tex]]"
---

# POST_ACTION_UPDATE — L02 renotation EXECUTED

## Akcja wykonana

Sesja 2026-05-04 wykonała **konkretną renotację** opisaną w
[[README.md]] rekomendacji.

### Edytowany plik

[[../../axioms/notacja/dodatekA_notacja.tex]]: dodana nowa subsekcja
`app:A-beta-gamma-distinction` (przed `app:A-wymiary`).

### Co dodano

**§A.X (nowa subsekcja, ~80 linii):** „Konwencja notacyjna: β, γ jako
współczynniki GL vs wykładniki krytyczne WF"

Zawiera:
1. **(I) (β, γ)_GL** — współczynniki GL w fazie złamanej:
   - $V_{\rm eff}(g) = (\beta/3)g^3 - (\gamma/4)g^4$
   - Wymiar $[\mathrm{L}^{-2}]$
   - $(\beta/\gamma)_{\rm GL} = 1$ z warunku próżniowego (algebraiczna tożsamość)
   - Lista miejsc użycia: r. pola, N0-2/5/6, LK-1d, G.0 closure

2. **(II) (β, γ)_WF** — wykładniki krytyczne w punkcie stałym Wilsona-Fishera:
   - $\beta_{\rm WF} \approx 0.326$, $\gamma_{\rm WF} \approx 1.237$
   - $(\beta/\gamma)_{\rm WF} \approx 0.264$ (3D Ising)
   - Lista miejsc użycia: M3 MK-RG, M8 NPRG, klasa uniwersalności

3. **„Dlaczego nie ma sprzeczności"** — TGP fizycznie nie żyje na T_c,
   tylko w fazie złamanej; obie wielkości pod tą samą nazwą NIE są
   równoważne.

4. **Konwencja w tekście:**
   - Domyślnie β, γ bez subskryptu = (I) GL coefficients
   - Przy odwołaniach do M3-M8 lub klasy uniwersalności — explicit
     subskrypt $\beta_{\rm WF}, \gamma_{\rm WF}$
   - Cross-reference do FOUNDATIONS §7 i audytu L02

### Subsekcja `app:A-wymiary` (istniejąca) — addytywne uzupełnienie

Dodano explicit acknowledgment:

> Wszystkie $\beta, \gamma$ w niniejszej tablicy odnoszą się do
> (I) GL coefficients (patrz~\ref{app:A-beta-gamma-distinction});
> WF exponents $\beta_{\rm WF}, \gamma_{\rm WF}$ są używane wyłącznie
> w kontekście M3-M8 i klasy uniwersalności.

## Charakter edycji: NON-BREAKING

- **Addytywne**: nowa subsekcja przed istniejącą, bez modyfikacji
  istniejących sekcji
- **Bez psucia cross-references**: wszystkie istniejące `\ref{}` do
  `app:A-wymiary` pozostają działające
- **Bez zmiany formuł**: wszystkie istniejące tablice zachowane
- **Krótki audit-aware comment block** w nagłówku subsekcji (LaTeX `%`)

## Pliki zmienione

| Plik | Zakres zmian | Charakter |
|------|--------------|-----------|
| `axioms/notacja/dodatekA_notacja.tex` | + ~80 lin (nowa subsekcja `app:A-beta-gamma-distinction`) + 4 lin annotation w `app:A-wymiary` | NON-BREAKING addytywne |

## Pliki dotknięte (rekomendacje przyszłe — NIE wykonane w tej sesji)

Zgodnie z rekomendacją L02 README, pełne wcielenie wymaga:
- `core/sek08*.tex` — odwołania do M3-M8 powinny być oznaczone explicit
  $\beta_{\rm WF}$ (estymata: 5-10 lokalizacji)
- `research/op1-op2-op4/M{3..8}_*.md` — explicit subskrypty (estymata:
  20+ lokalizacji)
- `TGP_FOUNDATIONS.md` §7 — może odwołać do nowej subsekcji
- `scripts/lk1d_beta_gamma_ratio.c` — komentarze (cosmetic)

To jest mechaniczna grep+replace, ~1 dzień pracy. **NIE wymaga
prerequisitu** — tylko wykonania.

## Status post-action

| L02 NEEDS | Status post-action |
|-----------|---------------------|
| (rekomendacja README): renotacja w `dodatekA_notacja.tex` | **EXECUTED** ✓ |
| Update `TGP_FOUNDATIONS.md` §7 | pending (opcjonalne) |
| Globalny grep+replace w core/research/scripts | pending (mechaniczne) |

L02 jest **substantialnie zamknięty** — kanoniczna konwencja explicit
udokumentowana w pliku notacyjnym. Pozostałe lokalizacje mogą iść
sukcesywnie (low priority, cosmetic).

## Werdykt

L02 **EXECUTED w sesji 2026-05-04**. Konwencja $(\beta/\gamma)_{\rm GL}$
vs $(\beta/\gamma)_{\rm WF}$ jest teraz **explicit udokumentowana**
w `dodatekA_notacja.tex`. Każdy zewnętrzny czytelnik (referee, reviewer)
ma teraz jednoznaczne źródło dla różnicy.

## Cross-references

- [[README.md]] — pierwotny audit L02 (rekomendacja renotation)
- [[../../axioms/notacja/dodatekA_notacja.tex]] §`app:A-beta-gamma-distinction`
- [[../../TGP_FOUNDATIONS.md]] §7 (kontekst β/γ w M3-M8)
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LK-1d (β/γ → 1, GL coefficients)
