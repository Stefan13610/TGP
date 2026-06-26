---
title: "Phase 0 — pre-rejestracja LOCK: rozstrzygnięcie W-DRR-1. Czy istnieje znak/zakres B(Φ) taki, że disformalny operator fluktuacji jest JEDNOCZEŚNIE (i) ghost-free, (ii) gradient-stabilny (c_s²≥0), (iii) ekranujący strumień Ṗ_b, (iv) zgodny ze statycznym Vainshteinem rdzenia + c_T (GW170817) — czy reżim silnego ekranowania WYMUSZA patologię ⇒ sektor radiacyjny BROKEN przez STABILNOŚĆ. Werdykt: BROKEN / SIGN-PINNED."
date: 2026-06-14
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-14 (before any computation)
cycle: op-disformal-stability-2026-06-14
scoping: "[[../../meta/SCOPING_op-disformal-stability_2026-06-14.md]]"
parent_cycle: "[[../op-disformal-radiation-resolution-2026-06-13/Phase_FINAL_close.md]] CLOSED-RESOLVED UNDERDETERMINED"
resolves_doubt: "W-DRR-1 (op-disformal-radiation-resolution Phase1 §1/§8 + Phase_FINAL §4): Z^{rr}=2A-6bX<0 dla u>1/3; znak B niewyprowadzony (O12)"
authorization: "User 2026-06-14 (sesja #25): 'zająć się cyklem badawczym op-disformal-stability' (aktywacja cyklu zarejestrowanego w SCOPING §6 + README)"
cycle_category: "STABILITY-AUDIT (fast-kill; potencjalnie definitywny falsyfikator przez stabilność, niezależny od pełnego O12 i pinowania κ_E)"
parent_verdicts_LOCKED:
  - "PR-025 PSR J0737−3039 TRIGGERED (13 227σ skalar masywny / 2 646σ σ bezmasowy) — LOCKED nietykalne"
  - "op-gravitational-sector-survival INDETERMINATE — LOCKED"
  - "op-disformal-radiation-resolution UNDERDETERMINED (F-DRR-1 PARTIAL / F-DRR-2 UNPINNED / F-DRR-3 POSTULATE) — LOCKED"
  - "Operator fluktuacji EXACT: Z^{μν}=2(A−bX)η^{μν}−4b ∂^μφ̄ ∂^νφ̄ (Phase1 sympy, match=True) — LOCKED"
  - "Strumień: Ṗ_φ^LIVE=(1/u)·(1/6)P_GR, u=bX_bg/A (Phase1) — LOCKED (kontekst (iii))"
independent_of: "pełne O12 (tylko sign-pin B), pinowanie κ_E = C_σσ₀² (osobna droga domknięcia), op-nucleation-dimensionality"
anti_lakatos_lock: "INHERITED (survival → op-disformal-radiation-resolution); aktywny od Phase 0"
tags:
  - pre-registration
  - disformal-stability
  - ghost-gradient-instability
  - k-essence
  - W-DRR-1-followup
  - sign-of-B
  - fast-kill
  - anti-Lakatos-LOCKED
---

# Phase 0 — pre-rejestracja LOCK

> **Status:** Phase 0 = LOCK przed jakimkolwiek rachunkiem. ZERO sympy, ZERO danych
> obserwacyjnych w tym pliku. Każda kolejna faza (1, [2], FINAL) wymaga osobnego user „działaj".

## §1 — Pytanie

Cykl-rodzic ([[../op-disformal-radiation-resolution-2026-06-13/Phase_FINAL_close.md]]) orzekł
**UNDERDETERMINED**: disformalny Vainshtein **tłumi** strumień $\dot P_b$ czynnikiem kinetycznym
$1/u$ ($u=bX_{\rm bg}/A$, $b=B/M_*^4$), ale magnituda zależy od **niewyprowadzonego** $B(\Phi)$
(O12). Zarejestrował przy tym DOUBT **W-DRR-1**: radialna część operatora fluktuacji
$Z^{rr}=2A-6bX<0$ dla $u>1/3$ — możliwa **niestabilność gradientowa / ghost**, zależna od znaku $B$.

> **Pytanie wiodące:** Czy istnieje znak (i zakres) $B(\Phi)$ taki, że disformalny operator
> fluktuacji jest **JEDNOCZEŚNIE**:
> - **(i)** ghost-free (zdrowy znak członu czasowego $Z^{00}$ / kinetyki $L'$),
> - **(ii)** gradient-stabilny ($c_s^2\ge0$; brak Laplacian instability: $Z^{rr},Z^{\perp}$ zdrowe),
> - **(iii)** ekranujący strumień $\dot P_b$ (czynnik $1/|1-u|<1$ realny, nietrywialny),
> - **(iv)** zgodny ze statycznym Vainshteinem rdzenia ($r_V(\odot)$, $\gamma\to10^{-12}$) **oraz**
>   z prędkością tensorów $c_T$ (`prop:cT`, GW170817) i zawartością polaryzacyjną
>   (`prop:disformal-polarization`) — wszystko przy **tym samym** znaku $B$,
>
> — czy też reżim silnego ekranowania ($|u|\gtrsim1$, gdzie tłumienie jest znaczące) **WYMUSZA
> patologię**, co falsyfikuje sektor radiacyjny przez **STABILNOŚĆ** (ostrzej i taniej niż przez
> strumień)?

**Dwa pełnoprawne wyniki, oba ostre:**
- **BROKEN:** żaden znak $B$ nie jest jednocześnie zdrowy + ekranujący + zgodny z (iv) ⇒ sektor
  radiacyjny **sfalsyfikowany przez stabilność** (cofa UNDERDETERMINED ku falsyfikacji — przez
  patologię operatora, nie przez strumień).
- **SIGN-PINNED:** istnieje zdrowy znak $B$ ⇒ **pierwsze twarde ograniczenie na $B(\Phi)$**
  (wkład do domknięcia O12, warunek 2 z op-disformal-radiation-resolution §6); werdykt sektora
  pozostaje **UNDERDETERMINED**, ale węższy (przedział/znak $B$ zawężony).

Cykl **czysto strukturalno-rachunkowy**: zero danych obserwacyjnych w fazie wyprowadzania
(ograniczenia GW170817/PPN cytowane LOCKED z `rem:B-constraints`, **nie przeliczane**; znak $B$
wynika z rachunku znaków + zgodności z rdzeniem, NIE z dopasowania do danych).

## §2 — Konwencja (ZAFIKSOWANA PRZED analizą znaków — forbidden #4 scoping)

Źródło pomyłek znakowych; fiksujemy jawnie, raz, przed jakimkolwiek rachunkiem:

- **Sygnatura:** **mostly-plus** $\eta_{\mu\nu}=\mathrm{diag}(-,+,+,+)$ (dziedziczone z
  op-disformal-radiation-resolution Phase 1, gdzie operator $Z^{\mu\nu}$ wyprowadzono EXACT).
- **Tło statyczne:** $\partial_t\bar\phi=0$, gradient **radialny** $\partial_i\bar\phi=\hat r\,\partial_r\bar\phi$.
- **Znak $X$:** $X\equiv(\nabla\bar\phi)^2=g^{ij}\partial_i\bar\phi\,\partial_j\bar\phi>0$
  (przestrzenny gradient, dodatni w mostly-plus). **To jest reżim statyczny — różny od
  kosmologicznego $\dot{\bar\phi}\ne0$, gdzie $X<0$.** Rozróżnienie krytyczne (forbidden #3).
- **Parametr Vainshteina:** $u\equiv bX/A$ (statyczny, radialny). Znak $u$ = znak $b$ = znak $B$
  (bo $X>0$, $A=e^{2\delta\Phi/\PhiZero}>0$). **Silne ekranowanie ⇒ $|u|\gtrsim1$.**
- **Lagranżjan k-essence (LOCKED):** $L=A\,X-\tfrac{b}{2}X^2$, $L'=A-bX$, $L''=-b$.

## §3 — Wejścia LOCKED (read-only; modyfikacja któregokolwiek = nowy PR/O12, nie rescue)

| # | Input | Wartość / forma | Źródło LOCKED |
|---|---|---|---|
| S1 | Operator fluktuacji (EXACT) | $Z^{\mu\nu}=2(A-bX)\eta^{\mu\nu}-4b\,\partial^\mu\bar\phi\,\partial^\nu\bar\phi$, $L''=-b$ | op-disformal-radiation-resolution Phase1 §1 (sympy, match=True) |
| S2 | Składowe znaków na tle statycznym (do potwierdzenia sympy w Phase 1) | $Z^{00}\propto-(A-bX)=-A(1-u)$; $Z^{rr}=2A-6bX=2A(1-3u)$; $Z^{\perp}\propto A(1-u)$ | wyprowadzenie z S1 (pre-derywacja §6; rachunek nadrzędny) |
| S3 | Prędkość dźwięku skalara (Garriga–Mukhanov) | $c_s^2=\dfrac{L'}{L'+2XL''}=\dfrac{A-bX}{A-3bX}=\dfrac{1-u}{1-3u}$ | standardowa k-essence; do potwierdzenia sympy z S1 |
| S4 | Czynnik tłumienia strumienia (kontekst (iii)) | $\dot P_\phi^{\rm LIVE}=(1/|1-u|)\cdot\tfrac16 P_{\rm GR}$ (skaling $1/u$ dla $u\gg1$) | op-disformal-radiation-resolution Phase1 §3 (LOCKED) |
| S5 | Metryka disformalna LIVE | $g_{\mu\nu}=A(\Phi)\eta_{\mu\nu}+\dfrac{B(\Phi)}{M_*^4}\partial_\mu\Phi\,\partial_\nu\Phi$, $A=e^{2\delta\Phi/\PhiZero}$ | sek08 `hyp:disformal` (eq:disformal-metric-formalizm) |
| S6 | Prędkość tensorów (kontekst (iv)) | $c_T^2=\dfrac{A\,c_0^2}{A+\frac{B}{M_*^4}(\partial\bar\Phi)^2_\perp}$ — tło kosmiczne $c_T=c_0$ EXACT; lokalnie zależy od znaku $B$ | sek08 `prop:cT` (eq:cT) |
| S7 | Zawartość polaryzacyjna (kontekst (iv)) | 6 modów $E(2)$; mody tensorowe $h_+,h_\times$ z $B\,\partial_i\bar\Phi\,\delta\Phi$ — istnienie zależy od $B\ne0$ | sek08 `prop:disformal-polarization` |
| S8 | Ograniczenia obserwacyjne na $B$ (cytowane, NIE przeliczane) | GW170817: $|B|\dot{\bar\Phi}^2/M_*^4\lesssim10^{-15}$; PPN (Cassini/LLR); spójność $\bar\Phi=\PhiZero\Rightarrow$ człon znika | sek08 `rem:B-constraints` (O12, priorytet 1) |
| S9 | Status $B(\Phi)$ | **niewyprowadzony — otwarty problem O12** (znak i magnituda); LOCKED z op-disformal-radiation-resolution FINAL §6 warunek 2 | sek08 `rem:B-constraints`; FINAL §6 |
| S10 | Skala Vainshteina statyki | $r_V(\odot)\sim5\times10^6$ AU; $\gamma_{\rm PPN}$ ekranowane do $\sim10^{-12}$ | recon §3; op-disformal-radiation-resolution Phase0 R12 |

**Uwaga reuse:** zero nowych stałych. Wszystkie symbole ($A,B,M_*,X,b,u,c_0,\PhiZero$) z rdzenia /
cyklu-rodzica. $C_\sigma,\sigma_0,\kappa_E$ **NIE** wchodzą (to osobna droga domknięcia — warunek 1).

## §4 — Zbiór testów CLOSED (zamknięty w Phase 0; żaden nowy mid-cycle)

| Ti | Test | Obiekt rachunku | A priori zagrożony wynik |
|---|---|---|---|
| **T-STA-A** | **No-ghost.** Znak członu czasowego $Z^{00}$ / kinetyki $L'=A-bX$ w reżimie ekranowania $|u|\gtrsim1$, jako funkcja sign($B$). | $Z^{00}$, $L'$ vs $u$ i sign($B$) | $u>1$ (B>0, silne ekran.) ⇒ $L'<0$ ⇒ **ghost** |
| **T-STA-B** | **Gradient.** $c_s^2=\frac{1-u}{1-3u}\ge0$ ORAZ $Z^{rr}=2A(1-3u)\ge0$, $Z^{\perp}=A(1-u)\ge0$ w reżimie ekranowania, jako funkcja sign($B$). | $c_s^2$, $Z^{rr}$, $Z^{\perp}$ vs $u$ i sign($B$) | $1/3<u<1$ ⇒ $c_s^2<0$ ⇒ **niestab. gradientowa** |
| **T-STA-C** | **Zgodność znaku.** Czy znak $B$ dający (i)+(ii)+(iii) jest **zgodny** ze znakiem, którego rdzeń wymaga dla (iv): $c_T^2$ (S6, brak nadświetlności/Laplaciana mody tensorowej), polaryzacje (S7), statyczny $r_V$/$\gamma$ (S10). | porównanie sign($B$)_zdrowy-radiacyjny vs sign($B$)_wymagany-rdzeń | **konflikt znaku** (np. radiacja chce $B<0$, $c_T$ chce $B>0$) |
| **T-STA-D** | **Agregat** (WYLICZONY z flag A∧B∧C) | reguła §5.1 | — |

**Trzy reżimy (zamknięte; mapa terenu, NIE werdykt):** $u<1/3$ (zdrowy, ekran. słabe) /
$1/3<u<1$ (gradient instab.) / $u>1$ (ghost). Pytanie: czy istnieje **znak/zakres** $B$ unikający
reżimów patologicznych **przy** zachowaniu (iii)+(iv). Zbiór CLOSED — żaden nowy test/reżim po LOCK.

## §5 — Falsyfikatory LOCKED (reguły IMMUTABLE; werdykt WYLICZONY z flag)

| ID | Test | Reguła decyzyjna (IMMUTABLE) |
|---|---|---|
| **F-STA-A** | T-STA-A — no-ghost. Sympy: znak $Z^{00}$/$L'$ jako funkcja $u$, sign($B$). | Flaga ∈ {`HEALTHY` (istnieje znak/zakres $B$ ekranujący z $L'>0$, $Z^{00}$ zdrowe), `GHOST-FORCED` (każdy znak $B$ dający (iii) wymusza $L'<0$ w reżimie ekran.)}. |
| **F-STA-B** | T-STA-B — gradient. Sympy: $c_s^2$, $Z^{rr}$, $Z^{\perp}$ jako funkcja $u$, sign($B$). | Flaga ∈ {`HEALTHY` (istnieje znak/zakres $B$ ekranujący z $c_s^2\ge0$ i $Z^{rr},Z^{\perp}\ge0$), `GRADIENT-FORCED` (każdy ekranujący znak $B$ wymusza $c_s^2<0$ lub $Z^{rr}<0$)}. |
| **F-STA-C** | T-STA-C — zgodność znaku z rdzeniem (iv). Sympy/analiza: sign($B$)_zdrowy (z A∧B) vs sign($B$)_rdzeń (S6 $c_T$, S7 polaryzacje, S10 $r_V$/$\gamma$). | Flaga ∈ {`CONSISTENT` (ten sam znak $B$ spełnia (i)+(ii)+(iii) **i** (iv)), `SIGN-CONFLICT` (znak zdrowy-radiacyjny ≠ znak wymagany przez rdzeń)}. |
| **F-STA-D** | Agregat (WYLICZONY z F-STA-A ∧ -B ∧ -C). | Reguła §5.1. Próg liczbowy NIE dotyczy (klasyfikacyjne); zero hardcoded T_pass. |

### §5.1 — Reguła agregatu F-STA-D (IMMUTABLE)

```
patologia_nieusuwalna := (F-STA-A = GHOST-FORCED) ∨ (F-STA-B = GRADIENT-FORCED)
konflikt_znaku        := (F-STA-C = SIGN-CONFLICT)

broken     := patologia_nieusuwalna ∨ konflikt_znaku
sign_pinned := (F-STA-A = HEALTHY) ∧ (F-STA-B = HEALTHY) ∧ (F-STA-C = CONSISTENT)

if broken:        ⇒ op-disformal-stability → BROKEN
                    (sektor radiacyjny sfalsyfikowany przez STABILNOŚĆ;
                     cofa UNDERDETERMINED ku FALSIFIED-AS-AXIOMATIZED — przez patologię operatora)
elif sign_pinned: ⇒ op-disformal-stability → SIGN-PINNED
                    (zdrowy znak/zakres B wyprowadzony ⇒ pierwsze twarde ograniczenie na B(Φ),
                     wkład do O12 warunek 2; werdykt sektora pozostaje UNDERDETERMINED, węższy)
else:             ⇒ UNDERDETERMINED + R1 flag (niedomknięcie analizy znaków; ewentualny spawn)
```

- **BROKEN** (HONEST_NEGATIVE = pełnoprawny PASS): żaden znak $B$ nie jest zdrowy+ekranujący+zgodny.
- **SIGN-PINNED** (wymaga 3/3: A=HEALTHY ∧ B=HEALTHY ∧ C=CONSISTENT; 2/3 NIE wystarcza): istnieje
  zdrowy znak, zgodny z rdzeniem ⇒ sign-pin $B$ (np. „$B<0$ na tle statycznym" lub „$B>0$") —
  realny wkład do domknięcia O12, ale sektor wciąż UNDERDETERMINED (magnituda + κ_E + M_* otwarte).

## §6 — Pre-derywacja analityczna (OCZEKIWANIE, NIE próg — rachunek nadrzędny)

Zapisane PRZED Phase 0; do obalenia rachunkiem. Struktura znaków (z S2/S3, tło statyczne, $X>0$):

| reżim | $1-u$ | $1-3u$ | $c_s^2$ | $L'$ | diagnoza |
|---|---|---|---|---|---|
| $u<1/3$ ($B>0$ słabe) | + | + | + | + | zdrowy, ale ekranowanie słabe ($1/|1-u|\sim1$) |
| $1/3<u<1$ ($B>0$) | + | − | **<0** | + | **niestabilność gradientowa** |
| $u>1$ ($B>0$ silne) | − | − | + | **−** | **ghost** ($L'<0$) |
| $u<0$ ($B<0$, dowolna magnituda) | + | + | $\frac{1+|u|}{1+3|u|}\in(0,1)$ | + | **zdrowy + tłumienie $1/(1+|u|)<1$ działa** |

**Napięcie centralne (pre-flagowane, do testu — NIE claim):**
1. Dla $b>0$ ($B>0$) **silne ekranowanie** ($u\gtrsim1$) leży dokładnie w reżimie ghost/gradient.
2. **Hipoteza robocza:** zdrowy + ekranujący operator skalara wymaga **$B<0$** (wtedy $u<0$, oba
   znaki zdrowe, $c_s^2\in(0,1)$, a tłumienie $1/(1+|u|)<1$ wciąż realne).
3. **ALE** rdzeń `prop:cT` (S6): $c_T^2=A c_0^2/\big(A+\tfrac{B}{M_*^4}(\partial\bar\Phi)^2_\perp\big)$
   z $(\partial\bar\Phi)^2_\perp>0$ ⇒ dla $B<0$ mianownik $<A$ ⇒ $c_T^2>c_0^2$ (**nadświetlna** moda
   tensorowa) lub utrata hiperboliczności; dla $B>0$ ⇒ $c_T^2\le c_0^2$ (subluminal). Tj. rdzeń
   **preferuje $B>0$** dla zdrowej propagacji tensorów — **przeciwnie** do hipotezy radiacyjnej (2).
4. **Stąd centralne ryzyko F-STA-C = `SIGN-CONFLICT`** (stabilny skalar radiacyjny chce $B<0$;
   zdrowy $c_T$ tensorów chce $B>0$ — różne projekcje $X$, **ten sam** znak $B$) ⇒ kandydat BROKEN.

**Zastrzeżenia (rachunek może obalić oczekiwanie — raportować wprost):**
- $c_T$ używa składowej **transwersalnej** $(\partial\bar\Phi)^2_\perp$ (silnie tłumionej w typowej
  geometrii detekcji — S6 (iii)), a $Z^{rr}$ — składowej **radialnej**. Czy napięcie znaku jest
  **nieusuwalne** (te same $B$, różne projekcje), czy reżimy się **rozłączają geometrycznie** —
  rozstrzyga Phase 1/2 rachunkiem, nie tą tabelą.
- Magnituda $|u|$ na orbicie J0737 vs reżim, w którym tłumienie jest „znaczące" — czy zdrowy zakres
  $B<0$ daje **nietrywialne** tłumienie (iii), czy tylko $1/(1+|u|)\approx1$ (ekranowanie iluzoryczne).

**Jeśli rachunek da `CONSISTENT` (sign-pin) — to silny wynik pozytywny (pierwsze ograniczenie O12).
Jeśli `SIGN-CONFLICT`/`*-FORCED` — pełnoprawna falsyfikacja przez stabilność. Oba pre-deklarowane.**

## §7 — Forbidden moves (IMMUTABLE; INHERITED + cyklowe)

1. Zakaz rewizji werdyktów LOCKED: PR-025 (13 227σ/2 646σ), survival INDETERMINATE,
   op-disformal-radiation-resolution UNDERDETERMINED + flagi (F-DRR-1/2/3), operator $Z^{\mu\nu}$ EXACT,
   det J=2. Zmiana = osobny PR/O12.
2. **Zakaz strojenia $B$/$M_*$/κ_E do danych LUB do „zdrowego" wyniku** (anti-Lakatos): znak $B$
   wynika z rachunku znaków + zgodności z rdzeniem, **NIE** jest dobierany pod werdykt (ani BROKEN, ani SIGN-PINNED).
3. **Zakaz mylenia ekranowania statycznego** (rdzeń $r_V$, $\gamma$, tło $X>0$) **z radiacyjnym**
   oraz tła **statycznego** ($\dot{\bar\phi}=0$, $X>0$) z **kosmologicznym** ($\dot{\bar\phi}\ne0$, $X<0$);
   oba reżimy muszą być spójne z **tym samym** znakiem $B$ (istota F-STA-C).
4. **Zafiksować konwencję sygnatury + znak $X$ PRZED analizą znaków** — wykonane §2 (mostly-plus,
   $X=(\nabla\bar\phi)^2>0$ statyczne). Żaden znak nie zmienia się mid-cycle.
5. **Zakaz domykania O12 w całości** (to osobny, większy problem) — tu **wyłącznie sign-pin $B$**;
   magnituda $B(\Phi)$, $C_\sigma$/κ_E, mikro-derywacja $M_*$ **poza zakresem** (forbidden no-mix).
6. Budżet nowych stałych: **0**. Zero danych obserwacyjnych w fazie wyprowadzania (S8 cytowane LOCKED,
   nie przeliczane); ewentualne porównanie z danymi tylko w oznaczonej fazie, jeśli w ogóle.
7. Zakaz hardcoded T_pass; flagi F-STA-A/B/C/D wyliczane z jawnego rachunku znaków (lekcja GST #19 /
   op-disformal-radiation-resolution forbidden #7).
8. **Symetryczna dyscyplina:** zakaz rescue-by-tuning (#2) **ORAZ** zakaz przedwczesnego negatywu —
   BROKEN dopuszczalny tylko po jawnym rachunku (znaki $Z^{\mu\nu}$, $c_s^2$, zgodność $c_T$), nie z
   samej tabeli pre-derywacji §6.
9. Zakaz miękkiego domknięcia: BROKEN / SIGN-PINNED — oba pełnoprawne; żadnego łagodzenia ku
   „prawie uratowane" ani ku „prawie sfalsyfikowane".
10. `SIGN-PINNED` wymaga 3/3 (A=HEALTHY ∧ B=HEALTHY ∧ C=CONSISTENT); 2/3 = non-pinned ⇒ broken lub R1.
11. Zakaz altitude-creep do budowy v2 ani do domknięcia O12; ten cykl = werdykt W-DRR-1 (sign-pin B)
    + (jeśli BROKEN) wkład do warunków brzegowych. v2/O12 = osobny program.

## §8 — Risk register

- **R-sign-projection-conflation (HIGH):** centralna pułapka — utożsamienie projekcji radialnej
  ($Z^{rr}$, gradient skalara) z transwersalną ($c_T$, $(\partial\bar\Phi)^2_\perp$). Czy to **ten sam**
  znak $B$ w obu (⇒ konflikt nieusuwalny), czy reżimy rozłączają się geometrycznie. Mitygacja: jawny
  rachunek obu projekcji z **jednej** metryki S5; F-STA-C porównuje sign($B$), nie magnitudy.
- **R-static-cosmo-conflation (HIGH):** znak $X$ różny na tle statycznym ($X>0$) vs kosmologicznym
  ($X<0$) ⇒ znak $u$ i diagnoza odwracają się. Mitygacja: konwencja §2 zafiksowana; każdy reżim
  oznaczony jawnie (forbidden #3).
- **R-rescue-temptation (MED):** po UNDERDETERMINED rodzica pokusa ogłoszenia SIGN-PINNED (znak $B<0$
  „ratuje"). Mitygacja: F-STA mechaniczne z flag; SIGN-PINNED wymaga 3/3 **łącznie z** C=CONSISTENT
  (zgodność z $c_T$ rdzenia), nie tylko zdrowia radiacyjnego.
- **R-premature-negative (MED, symetryczne):** łatwy BROKEN z samej tabeli §6 bez rachunku zgodności
  geometrycznej. Mitygacja: forbidden #8; F-STA-C musi policzyć/porównać znaki, nie zacytować napięcia.
- **R-screening-triviality (MED):** zdrowy zakres $B<0$ może dawać tłumienie $1/(1+|u|)\approx1$
  (ekranowanie iluzoryczne) ⇒ (iii) formalnie spełnione, fizycznie puste. Mitygacja: T-STA wymaga
  **nietrywialnego** $1/|1-u|<1$ w reżimie istotnym dla Ṗ_b (warunek (iii) z magnitudą, nie tylko znakiem).
- **R-O12-scope-creep (MED):** pokusa „przy okazji" wyprowadzić magnitudę $B$ lub domknąć O12.
  Mitygacja: forbidden #5/#11; wyłącznie sign-pin.

## §9 — Anti-Lakatos checklist (Phase 0)

- ✓ Phase 0 LOCKED przed jakimkolwiek rachunkiem (zero sympy, zero danych w tym pliku).
- ✓ Konwencja (sygnatura mostly-plus + znak $X$ statyczny) zafiksowana §2 PRZED analizą znaków.
- ✓ Zbiór testów {T-STA-A,B,C,D} CLOSED; trzy reżimy ($u<1/3$ / $1/3<u<1$ / $u>1$) + gałąź $u<0$ zamknięte; żaden nowy mid-cycle.
- ✓ Flagi F-STA-A/B/C/D wyliczane z flag; 0 hardcoded T_pass; reguła agregatu §5.1 IMMUTABLE.
- ✓ Dwa wyniki (BROKEN / SIGN-PINNED) pre-deklarowane jako pełnoprawne; wymóg 3/3 dla SIGN-PINNED (#10).
- ✓ Pre-derywacja §6 = oczekiwanie, nie próg; napięcie znaku ($B<0$ radiacja vs $B>0$ $c_T$) jawnie pre-flagowane jako ryzyko, NIE werdykt.
- ✓ Werdykty poprzedników LOCKED i nietykalne (forbidden #1); operator $Z^{\mu\nu}$ EXACT reused, nie wyprowadzany od nowa.
- ✓ Budżet nowych stałych: 0; zakres = wyłącznie sign-pin $B$ (O12 jako całość poza zakresem, #5).
- ✓ Symetryczna ochrona: anty-rescue-by-tuning (#2/#10) ORAZ anty-przedwczesny-negatyw (#8).
- ✓ Rozróżnienie statyka/radiacja oraz statyka/kosmologia (znak $X$) zafiksowane jako reguła (#3).

## §10 — Plan faz (każda wymaga osobnego user „działaj")

- **Phase 0** — ten plik (LOCK). ✅
- **Phase 1 — T-STA-A + T-STA-B (no-ghost + gradient; F-STA-A, F-STA-B):** sympy — z operatora $Z^{\mu\nu}$
  (S1) na tle statycznym: znaki $Z^{00}$, $Z^{rr}$, $Z^{\perp}$, $L'$ i $c_s^2$ (S3) jako funkcje $u$ i
  sign($B$); identyfikacja zdrowego znaku/zakresu $B$ dla (i)+(ii) **z** warunkiem ekranowania (iii).
  Output: flagi F-STA-A, F-STA-B.
- **Phase 2 — T-STA-C (zgodność znaku z rdzeniem; F-STA-C):** sympy/analiza — sign($B$)_zdrowy (z Phase 1)
  vs sign($B$)_rdzeń z $c_T$ (S6, brak nadświetlności/Laplaciana mody tensorowej), polaryzacji (S7),
  statyki $r_V$/$\gamma$ (S10). Czy ten sam znak $B$ spełnia (iv). Output: flaga F-STA-C.
  *(Jeśli Phase 1 da `*-FORCED` ⇒ broken już na A/B; Phase 2 może być skrócona do potwierdzenia.)*
- **Phase FINAL** — agregat F-STA-D (§5.1, z flag); closure; (jeśli ratyfikowane) propagacja:
  aktualizacja W-DRR-1 w op-disformal-radiation-resolution FINAL §4, FOUNDATIONS §3.6.10.6,
  `rem:B-constraints`/O12 (sign-pin lub falsyfikacja), STATE. Dyspozycja PR — wyłącznie user.

## §11 — Numer PR

RESERVED — przydzieli Phase FINAL. Scenariusze:
- **BROKEN** → adnotacja **FALSIFIED-BY-STABILITY** sektora radiacyjnego (cofa UNDERDETERMINED;
  adnotacja przy op-disformal-radiation-resolution + FOUNDATIONS; bez nowego PR obserwacyjnego).
- **SIGN-PINNED** → adnotacja strukturalna **„sign($B$) pinned"** w `rem:B-constraints`/O12 (pierwsze
  twarde ograniczenie na $B(\Phi)$); sektor pozostaje UNDERDETERMINED (węższy); ewentualny kandydat PR
  — wyłącznie decyzja user (anty-scope-creep).
