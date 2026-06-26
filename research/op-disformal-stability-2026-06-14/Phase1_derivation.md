---
title: "Phase 1 — T-STA-A (no-ghost) + T-STA-B (gradient): zdrowy (ghost-free + c_s²≥0) ORAZ ekranujący operator fluktuacji disformalnej ISTNIEJE i wymaga JEDNOZNACZNIE B<0 (u<0). Dla B>0 ekranowanie wymaga u>2 = głęboki ghost, a region zdrowy u<1/3 WZMACNIA strumień (S>1). F-STA-A = HEALTHY, F-STA-B = HEALTHY; kandydat sign-pin: B<0. Werdykt cyklu przenosi się na F-STA-C (zgodność B<0 z c_T rdzenia)."
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-disformal-stability-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14: 'działaj z Phase 1'"
phase0_lock: "[[./Phase0_balance.md]]"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "10/10 PASS; 0 hardcoded T_pass; 0 danych; OBA znaki B testowane; werdykt z flag. EXACT: Z^{μν} match=True (reuse S1); Z^{00}=-2A(1-u), Z^{rr}=2A(1-3u), Z^{⊥}=2A(1-u); c_s²=(1-u)/(1-3u)."
flag_F_STA_A: "HEALTHY (istnieje zakres B<0 ekranujący z L'>0, Z^{00}/Z^{⊥} zdrowe)"
flag_F_STA_B: "HEALTHY (istnieje zakres B<0 ekranujący z c_s²∈(1/3,1), Z^{rr},Z^{⊥}>0)"
sign_pin_candidate: "B<0 (u<0) — JEDNOZNACZNY zakres zdrowy+ekranujący; warunkowy na F-STA-C"
feeds: "F-STA-D basin: NOT-(*-FORCED) (patologia nieusuwalna OBALONA — zdrowy znak istnieje) ⇒ werdykt cyklu rozstrzyga F-STA-C (Phase 2): zgodność B<0 z rdzeniem (c_T, polaryzacje, r_V)"
anti_lakatos_lock: PRESERVED
carries_tension: "B<0 (zdrowy radiacyjnie) vs prop:cT preferuje B>0 (subluminal tensory) — SIGN-CONFLICT kandydat (Phase0 §6, R-sign-projection-conflation HIGH); rozstrzyga Phase 2"
---

# Phase 1 — T-STA-A (no-ghost) + T-STA-B (gradient)

## §0 — Verdict at a glance

| Element | Wynik |
|---|---|
| **Operator** (reuse S1, EXACT) | $Z^{\mu\nu}=2(A-bX)\eta^{\mu\nu}-4b\,\partial^\mu\bar\phi\,\partial^\nu\bar\phi$ — match=True; tło statyczne $X=G^2>0$ |
| **Składowe (S2 potwierdzone)** | $Z^{00}=-2A(1-u)$, $Z^{rr}=2A(1-3u)$, $Z^{\perp}=2A(1-u)$; $u=bX/A$ |
| **No-ghost** (T-STA-A) | $L'=A(1-u)>0\iff u<1$; dla $u>1$ ($B>0$ silne): $L'<0\Rightarrow$ **ghost** |
| **Gradient** (T-STA-B) | $c_s^2=\frac{1-u}{1-3u}$; stabilny $\iff u<1/3$; dla $1/3<u<1$: $Z^{rr}<0,c_s^2<0\Rightarrow$ **niestab. gradientowa** |
| **Ekranowanie** $S=1/\|1-u\|$ | $S<1$ (genuine suppression) $\iff u<0$ **lub** $u>2$; dla $0<u<2$: $S>1$ (**wzmocnienie**!) |
| **health $\wedge$ screening** | $\{u<1/3\}\cap\{u<0\,\cup\,u>2\}=\boxed{u<0}$ — **jednoznacznie $B<0$** |
| **gałąź $B>0$** | health$\wedge$screening = **∅** (ekranowanie $B>0$ wymaga $u>2$ = głęboki ghost) |
| **gałąź $B<0$** | health$\wedge$screening = **cały $u<0$** ($c_s^2\in(\tfrac13,1)$ subluminal, $S=1/(1+\|u\|)<1$) |
| **F-STA-A / F-STA-B** | **HEALTHY / HEALTHY** — patologia NIE jest wymuszona; zdrowy znak istnieje ($B<0$) |
| **Werdykt cyklu** | **nierozstrzygnięty na A/B** ⇒ przenosi się na **F-STA-C** (zgodność $B<0$ z rdzeniem) |

**sympy: 10/10 PASS.** Werdykt z flag; 0 hardcoded; oba znaki $B$ testowane (nie założone).

## §1 — Operator i składowe znakowe (T-STA-A1, EXACT — reuse S1)

Z zalockowanego $L=A\,X-\tfrac{b}{2}X^2$ (S1), różniczkując po symbolicznych $\partial_\mu\phi$ i
podstawiając tło statyczne $\partial_\mu\bar\phi=(0,G,0,0)$ ($\partial_t\bar\phi=0$, gradient radialny;
mostly-plus, **forbidden #4 — konwencja zafiksowana przed rachunkiem**):

$$X=\eta^{\mu\nu}\partial_\mu\bar\phi\,\partial_\nu\bar\phi=+G^2>0\quad(\text{reżim statyczny; }X<0\text{ byłoby kosmologiczne}).$$
$$Z^{\mu\nu}=2(A-bX)\eta^{\mu\nu}-4b\,\partial^\mu\bar\phi\,\partial^\nu\bar\phi\quad(\text{match}=\text{True}).$$

Składowe (z $u=bX/A$; **potwierdzają S2**, sympy):
$$Z^{00}=-2A(1-u),\qquad Z^{rr}=2A(1-3u),\qquad Z^{\perp}=2A(1-u).$$

To **reuse** operatora EXACT z op-disformal-radiation-resolution Phase 1 — nie nowy obiekt, zero nowych stałych.

## §2 — No-ghost (T-STA-A2) i gradient (T-STA-B1)

**No-ghost** (znak głównego członu kinetycznego, konwencja k-essence projektu $\text{no-ghost}\iff L'>0$):
$$L'=A-bX=A(1-u)>0\iff u<1.$$
Dla $u>1$ ($B>0$, silne ekranowanie) $L'<0$ ⇒ **ghost** (zła kinetyka).

**Gradient** (prędkość dźwięku Garriga–Mukhanov, sympy potwierdza S3):
$$c_s^2=\frac{L'}{L'+2XL''}=\frac{A-bX}{A-3bX}=\frac{1-u}{1-3u}.$$
Warunki stabilności gradientowej (sympy `solveset` nad $\mathbb R$):
- $c_s^2\ge0$: $u<1/3$ **lub** $u\ge1$;
- $Z^{rr}=2A(1-3u)\ge0$: $u\le1/3$ (mod radialny, wzdłuż gradientu);
- $Z^{\perp}=2A(1-u)\ge0$: $u\le1$ (mod transwersalny).

Przecięcie ⇒ **gradient-stabilny $\iff u<1/3$**. W oknie $1/3<u<1$: $Z^{rr}<0$ i $c_s^2<0$ ⇒
**niestabilność gradientowa**. (Region $u\ge1$ z $c_s^2\ge0$ jest wykluczony przez ghost $L'<0$.)

## §3 — Ekranowanie: gdzie strumień jest naprawdę TŁUMIONY (T-STA-B2)

Czynnik z op-disformal-radiation-resolution Phase 1 §3 (LOCKED, S4): $\dot P^{\rm LIVE}/\dot P^{\rm unscr}=S(u)=1/|1-u|$.
**Genuine suppression** ($S<1$) $\iff|1-u|>1\iff$ (sympy):
$$u<0\ \ (B<0)\qquad\text{lub}\qquad u>2\ \ (B>0,\ \text{głęboki ghost}).$$
**Krytyczne:** dla $0<u<2$ ($B>0$ słabe/pośrednie) $S>1$ — strumień jest **WZMACNIANY**, nie tłumiony
(jeszcze gorzej niż brak ekranowania). Czyli $B>0$ daje suppression **tylko** w reżimie ghost.

## §4 — Egzystencja zdrowego znaku (T-STA-C prep — rdzeń wyniku)

Pytanie Phase 0: **czy istnieje znak/zakres $B$ jednocześnie (i) ghost-free + (ii) gradient-stabilny + (iii) ekranujący?**
$$\text{full\_health}=\{u<1\}\cap\{u<\tfrac13\}=\{u<\tfrac13\},\qquad
\text{screening}=\{u<0\}\cup\{u>2\}.$$
$$\boxed{\text{health}\wedge\text{screening}=\{u<\tfrac13\}\cap(\{u<0\}\cup\{u>2\})=\{u<0\}.}$$

Rozbicie wg znaku $B$ (sympy):
- **Gałąź $B>0$ ($u>0$):** $\text{health}\wedge\text{screening}=\varnothing$. Ekranowanie $B>0$ wymaga $u>2$,
  ale tam $L'<0$ (ghost); region zdrowy $u<1/3$ daje $S>1$ (wzmocnienie). **Brak zdrowego ekranującego $B>0$.**
- **Gałąź $B<0$ ($u<0$):** $\text{health}\wedge\text{screening}=\{u<0\}$ (całość). Weryfikacja (T-STA-Bneg):
  dla $u=-w,\ w>0$: $c_s^2=\frac{1+w}{1+3w}\in(\tfrac13,1)$ (**subluminal**, zdrowy), $Z^{rr}=2A(1+3w)>0$,
  $Z^{\perp}=2A(1+w)>0$, $S=1/(1+w)<1$ (ekranuje). **Zdrowe ∀ magnitudy $w$.**

**Wniosek: zdrowy (no-ghost ∧ gradient) ORAZ ekranujący operator ISTNIEJE i wymaga JEDNOZNACZNIE $B<0$.**
Patologia (ghost/gradient) **NIE jest nieusuwalna** — jest wymuszana tylko, jeśli z góry narzucić $B>0$.
Ponieważ znak $B$ jest swobodny (niewyprowadzony, O12 — S9), **istnienie zdrowej gałęzi $B<0$ jest faktem strukturalnym**.

## §5 — Flagi (wyliczone z rachunku)

- **F-STA-A = HEALTHY:** istnieje zakres $B<0$ ekranujący z $L'>0$ (i $Z^{00},Z^{\perp}$ zdrowe). NIE `GHOST-FORCED`
  (ghost wymuszony byłby tylko gdyby każdy ekranujący znak dawał $L'<0$ — obalone gałęzią $B<0$).
- **F-STA-B = HEALTHY:** ta sama gałąź $B<0$ daje $c_s^2\in(\tfrac13,1)$, $Z^{rr},Z^{\perp}>0$. NIE `GRADIENT-FORCED`.

Reguła agregatu Phase 0 §5.1: $\text{patologia\_nieusuwalna}=(\text{A}=\text{GHOST-FORCED})\vee(\text{B}=\text{GRADIENT-FORCED})=\textbf{False}$.
⇒ **broken NIE jest jeszcze ustalone**; $\text{sign\_pinned}$ wymaga dodatkowo $\text{F-STA-C}=\text{CONSISTENT}$.
**Werdykt cyklu przenosi się na Phase 2 (F-STA-C).**

## §6 — Napięcie przeniesione do F-STA-C (Phase 2) — pre-flagowane, NIE werdykt

Zdrowy znak radiacyjny to **$B<0$**. Ale rdzeń `prop:cT` (S6):
$$c_T^2=\frac{A\,c_0^2}{A+\frac{B}{M_*^4}(\partial\bar\Phi)^2_\perp},\qquad(\partial\bar\Phi)^2_\perp>0.$$
Dla $B<0$ mianownik $<A\Rightarrow c_T^2>c_0^2$ (**moda tensorowa nadświetlna** / utrata hiperboliczności);
rdzeń **preferuje $B>0$** (subluminal). To dokładnie **`SIGN-CONFLICT` kandydat** z Phase 0 §6:
- skalar radiacyjny (mod **radialny** $Z^{rr}$, gradient) chce **$B<0$**;
- tensory rdzenia (składowa **transwersalna** $(\partial\bar\Phi)^2_\perp$) chcą **$B>0$**;
- **ten sam** znak $B$, różne projekcje $X$.

**Otwarte do Phase 2 (R-sign-projection-conflation, HIGH):** czy napięcie jest **nieusuwalne** (ta sama
metryka, oba warunki na sign $B$) ⇒ **BROKEN** (falsyfikacja przez stabilność), czy reżimy **rozłączają się
geometrycznie** (radial vs transwersal; $(\partial\bar\Phi)^2_\perp$ silnie tłumione w geometrii detekcji, S6(iii))
⇒ **SIGN-PINNED $B<0$**. Rozstrzyga jawny rachunek Phase 2, **nie** ta nota (forbidden #8).

## §7 — DOUBTS / uwagi

- **W-STA-1 (MED, NOWY):** zdrowy zakres $B<0$ daje $c_s^2\in(\tfrac13,1)$ — **skalar subluminal**, ale `prop:cT`
  sugeruje **tensory nadświetlne** dla tego samego $B<0$. Współistnienie subluminal-skalar / superluminal-tensor
  na tym samym tle = potencjalny problem przyczynowości (Phase 2 F-STA-C musi to rozstrzygnąć, nie wygładzić).
- **W-STA-2 (LOW):** „ekranowanie" $B<0$ jest w sensie **znaku** ($S<1$ ∀$u<0$); **magnituda** nietrywialnego
  tłumienia ($|u|\gtrsim1$ na orbicie) zależy od $|B|$/$M_*$ = O12 (poza zakresem, forbidden #5). Tu rozstrzygamy
  wyłącznie kierunek (suppression vs enhancement), nie wielkość.

## §8 — Anti-Lakatos compliance (Phase 1)

✓ Werdykt z flag (10/10 PASS); 0 hardcoded T_pass. ✓ **Oba znaki $B$ przetestowane** (nie założono $B<0$ —
wyszedł z przecięcia health∧screening). ✓ Operator $Z^{\mu\nu}$ EXACT **reused** (S1), nie wyprowadzany od nowa.
✓ Konwencja (mostly-plus, $X>0$ statyczne) zafiksowana w Phase 0 §2, użyta jawnie (forbidden #4). ✓ Rozróżnienie
statyka/kosmologia (znak $X$) i suppression/enhancement zachowane (forbidden #3). ✓ Zakres = **wyłącznie sign-pin**;
magnituda $B$/$M_*$/$\kappa_E$ poza (forbidden #5). ✓ **Anty-przedwczesny-negatyw:** patologia NIE ogłoszona —
rachunek pokazał zdrową gałąź $B<0$, więc broken-przez-A/B obalone (forbidden #8). ✓ Napięcie $c_T$ pre-flagowane
jako ryzyko F-STA-C, NIE jako werdykt. ✓ 0 danych obserwacyjnych. ✓ Budżet nowych stałych: 0.

---

**Phase 1 COMPLETE — F-STA-A = HEALTHY, F-STA-B = HEALTHY; kandydat sign-pin B<0.**
Werdykt cyklu (BROKEN vs SIGN-PINNED) rozstrzyga **Phase 2 / F-STA-C**: zgodność $B<0$ z rdzeniem
(`prop:cT` $c_T$, `prop:disformal-polarization`, statyczny $r_V$/$\gamma$). **Wymaga osobnego user „działaj".**
