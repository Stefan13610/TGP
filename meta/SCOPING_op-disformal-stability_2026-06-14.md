---
title: "SCOPING — op-disformal-stability: audyt ghost/niestabilności operatora fluktuacji disformalnego (W-DRR-1). Czy istnieje znak/zakres B(Φ) jednocześnie (i) ghost-free, (ii) gradient-stabilny, (iii) ekranujący Ṗ_b, (iv) zgodny ze statycznym Vainshteinem rdzenia — czy reżim ekranowania wymusza patologię ⇒ sektor radiacyjny BROKEN przez stabilność?"
date: 2026-06-14
type: meta-scoping
status: "PRE-PHASE-0 NOTE (nie pre-rejestracja; nie cykl; zero werdyktów; wymaga własnego Phase 0 + 'działaj')"
origin: "User 2026-06-14 (sesja #25): analiza wyników op-disformal-radiation-resolution (UNDERDETERMINED) → W-DRR-1 (znak B → Z^{rr}<0) zidentyfikowane jako najtańsza droga do DEFINITYWNEGO werdyktu. User: 'A + B' → rozpisać scoping stabilności."
parent_cycle: "[[../research/op-disformal-radiation-resolution-2026-06-13/]] CLOSED-RESOLVED UNDERDETERMINED (W-DRR-1 zarejestrowane jako MED, otwarte)"
resolves_doubt: "W-DRR-1 (op-disformal-radiation-resolution Phase1 §8 / FINAL §4)"
anti_lakatos_note: "Mapa terenu / pytanie strukturalne. Żadne twierdzenie nie nabywa statusu przez tę notę. Werdykty poprzedników (PR-025, survival, op-disformal-radiation-resolution) LOCKED. Pre-derywacja §2 = oczekiwanie, nie próg."
tags:
  - scoping
  - disformal-stability
  - ghost-gradient-instability
  - k-essence
  - W-DRR-1-followup
  - fast-kill
  - anti-Lakatos-LOCKED
---

# SCOPING — stabilność operatora disformalnego: czy ratunek jest ghostem?

## §1 — Pytanie (sformułowanie, doprecyzowane)

op-disformal-radiation-resolution ustalił, że disformalny Vainshtein **tłumi** strumień
$\dot P_b$ czynnikiem $1/u$ ($u=bX_{\rm bg}/A$, $b=B/M_*^4$), ale magnituda zależy od
niewyprowadzonego $B(\Phi)$ (O12), a κ_E niepinowane → **UNDERDETERMINED**. Zarejestrował
jednak DOUBT **W-DRR-1**: część operatora fluktuacji $Z^{rr}=2A-6bX<0$ dla $u>1/3$ — możliwa
**niestabilność gradientowa / ghost**, zależna od znaku $B$.

> **Czy istnieje znak (i zakres) $B(\Phi)$ taki, że disformalny operator fluktuacji jest
> JEDNOCZEŚNIE: (i) ghost-free (zdrowy znak członu czasowego), (ii) gradient-stabilny
> ($c_s^2\ge0$, brak Laplacian instability), (iii) ekranujący $\dot P_b$ (czynnik $1/|1-u|<1$),
> (iv) zgodny ze statycznym Vainshteinem rdzenia ($r_V(\odot)$, $\gamma\to10^{-12}$) — czy też
> reżim silnego ekranowania ($u\gg1$, gdzie tłumienie jest znaczące) WYMUSZA patologię,
> co falsyfikuje sektor radiacyjny przez STABILNOŚĆ (ostrzej niż przez strumień)?**

To **potencjalnie najtańszy decydujący falsyfikator** całego łańcucha: rozstrzyga znakiem
$B$, **nie** wymaga rozwiązania całego O12 ani pinowania κ_E. Wynik ostry w obie strony:
BROKEN (patologia nieusuwalna) lub sign-pin $B$ (zawęża O12, sektor zostaje UNDERDETERMINED
ale węższy).

## §2 — Pre-derywacja analityczna (R1 #18 — zapisana PRZED Phase 0; OCZEKIWANIE, nie próg)

Operator fluktuacji k-essence (LOCKED z op-disformal-radiation-resolution Phase 1, EXACT):
$$Z^{\mu\nu}=2L'(X)\eta^{\mu\nu}+4L''(X)\partial^\mu\bar\phi\,\partial^\nu\bar\phi,\quad L'=A-bX,\;L''=-b.$$
Dla statycznego tła (∂_tφ̄=0, gradient radialny; konwencja mostly-plus, $X=(\nabla\bar\phi)^2>0$),
z $u\equiv bX/A$:
- **Czasowy (ghost):** $Z^{00}\propto -(A-bX)=-A(1-u)$ — znak zdrowy gdy $1-u>0$.
- **Radialny (gradient, wzdłuż ∇φ̄):** $Z^{rr}=2A-6bX=2A(1-3u)$ — zdrowy gdy $1-3u>0$.
- **Poprzeczny:** $Z^{\perp}\propto(A-bX)=A(1-u)$.
- **Prędkość dźwięku (Garriga–Mukhanov):** $c_s^2=\dfrac{L'}{L'+2XL''}=\dfrac{A-bX}{A-3bX}=\dfrac{1-u}{1-3u}.$

**Struktura znaków (oczekiwana):**

| reżim | $1-u$ | $1-3u$ | $c_s^2$ | diagnoza |
|---|---|---|---|---|
| $u<1/3$ | + | + | + | zdrowy, ale ekranowanie słabe ($1/|1-u|\sim1$) |
| $1/3<u<1$ | + | − | **<0** | **niestabilność gradientowa** |
| $u>1$ | − | − | + (ale $L'<0$) | **ghost** (zły znak kinetyki) |

**Napięcie centralne (pre-flagowane):** silne tłumienie wymaga $|u|\gg1$ — a dla $b>0$
(tj. $X>0$, $B>0$) to dokładnie reżim ghost/niestabilności. **Hipoteza robocza (do testu,
NIE claim):** zdrowy + ekranujący wymaga $B<0$ — wtedy $u<0$, $1-u=1+|u|>0$, $1-3u=1+3|u|>0$,
$c_s^2=(1+|u|)/(1+3|u|)\in(0,1)$ zdrowe, a tłumienie $1/(1+|u|)<1$ wciąż działa. **Jeśli tak —
werdykt zależy od (iv): czy $B<0$ jest zgodne z tym, do czego rdzeń używa $B$** (polaryzacje GW
`prop:disformal-polarization`, statyczny $r_V$/$\gamma$, ograniczenia O12). Konflikt znaku ⇒ BROKEN;
zgodność ⇒ sign-pin (zawęża O12). **Rachunek nadrzędny nad tą hipotezą; konwencja sygnatury/znaku
X = pierwsza rzecz do zafiksowania w Phase 0.**

## §3 — Falsyfikatory (szkic do LOCK w Phase 0)

| ID | Test | Reguła (kandydat) |
|---|---|---|
| F-STA-A | No-ghost: znak członu czasowego $Z^{00}$ / kinetyki $L'$ w reżimie ekranowania $|u|\gtrsim1$ | ghost dla wszystkich znaków B dających ekranowanie ⟹ wkład do BROKEN |
| F-STA-B | Gradient: $c_s^2=\frac{1-u}{1-3u}\ge0$ i $Z^{rr},Z^{\perp}\ge0$ w reżimie ekranowania | $c_s^2<0$ nieusuwalne ⟹ wkład do BROKEN |
| F-STA-C | Zgodność sign($B$)_zdrowy z (iv): rdzeń `prop:disformal-polarization` + statyczny $r_V$/$\gamma$ + O12 `rem:B-constraints` | konflikt znaku ⟹ BROKEN; zgodność ⟹ sign-pin (UNDERDETERMINED węższy) |
| F-STA-D | Agregat (z flag) | (A∨B = patologia nieusuwalna) ∨ (C = konflikt) ⟹ **BROKEN** (sektor radiacyjny sfalsyfikowany przez stabilność); inaczej ⟹ **SIGN-PINNED** (wkład do domknięcia O12, werdykt sektora pozostaje UNDERDETERMINED) |

Werdykt klasyfikacyjny, wyliczany z flag; 0 hardcoded; 0 danych obserwacyjnych (czysto strukturalny).

## §4 — Format i koszt (wzór fast-kill)

1–2 fazy merytoryczne, sympy (znaki $Z^{\mu\nu}$, $c_s^2$ jako funkcje $u$ i sign($B$);
audyt zgodności znaku z rdzeniem), 0 nowych stałych, 0 danych. Reuse: operator $Z^{\mu\nu}$
EXACT z op-disformal-radiation-resolution Phase 1 (LOCKED). **Win-win:** BROKEN ⟹ definitywne
domknięcie sektora grawitacyjnego (cofa UNDERDETERMINED ku falsyfikacji — przez stabilność,
nie strumień); SIGN-PINNED ⟹ pierwsze twarde ograniczenie na $B(\Phi)$, realny wkład do O12.

## §5 — Forbidden moves (szkic)

1. Zakaz rewizji PR-025 / survival / op-disformal-radiation-resolution (LOCKED).
2. Zakaz strojenia $B$/$M_*$/κ_E do danych lub do „zdrowego" wyniku (anti-Lakatos; znak B
   wynika z rachunku + zgodności z rdzeniem, NIE jest dobierany pod werdykt).
3. Zakaz mylenia ekranowania statycznego (rdzeń $r_V$, $\gamma$) z radiacyjnym — oba muszą
   być spójne z tym samym znakiem B (to istota F-STA-C).
4. Zafiksować konwencję sygnatury + znak $X$ dla statycznego tła PRZED analizą znaków (źródło
   pomyłek; mostly-plus z op-disformal-radiation-resolution jako domyślne, jawnie zadeklarowane).
5. Zakaz domykania O12 w całości (to osobny, większy problem) — tu tylko sign-pin B.
6. Budżet nowych stałych: 0.

## §6 — Rejestracja follow-upu

**`op-disformal-stability`** — REGISTERED jako kandydat (NOT activated; wymaga własnego Phase 0
+ user „działaj"). [[../research/op-disformal-stability-2026-06-14/]] (README). §2–§5 = szkic
zakresu. Rozstrzyga W-DRR-1 z op-disformal-radiation-resolution. Niezależny od pełnego O12
(tylko sign-pin) i od pinowania κ_E (osobna droga domknięcia).

> **Relacja do warunków domknięcia (op-disformal-radiation-resolution §6):** ten cykl atakuje
> warunek 2 (B(Φ)) wyłącznie od strony ZNAKU/stabilności — najtańszy fragment. Warunki 1 (pin κ_E)
> i 3 (mikro-derywacja M_*) pozostają osobne. Pełna falsyfikowalność sektora wymaga wszystkich trzech.
