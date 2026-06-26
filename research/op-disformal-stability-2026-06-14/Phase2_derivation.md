---
title: "Phase 2 — T-STA-C (zgodność znaku B z rdzeniem, warunek (iv)): SIGN-CONFLICT. Zdrowy radiacyjnie znak B<0 (Phase 1) jest NIEZGODNY z rdzeniem: prop:cT wymaga B≥0 (subluminal+hiperboliczny tensor). Przecięcie znaków = ∅. Nieusuwalne w strefie falowej (k̂∥∇Φ̄ ⇒ u_⊥=u; próg niestabilności tensora |u|=1 = r_V; wewnątrz r_V, J0737, B<0 ⇒ c_T²<0). Implikowany agregat F-STA-D = BROKEN."
type: phase_result
status: PHASE2_COMPLETE
phase: 2
cycle: op-disformal-stability-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14: 'działaj z Phase 2'"
phase0_lock: "[[./Phase0_balance.md]]"
phase1_result: "[[./Phase1_derivation.md]] (F-STA-A=HEALTHY, F-STA-B=HEALTHY; sign-pin kandydat B<0)"
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
sympy: "6/6 PASS; 0 hardcoded T_pass; 0 danych; OBA znaki B testowane; werdykt z flag. Kluczowe: tensor zdrowy (denom>0 ∧ subluminal) ⟺ b≥0; scalar zdrowy+ekranujący ⟺ b<0; przecięcie ∅; u_⊥=u w strefie falowej (k̂∥∇Φ̄); próg |u|=1=r_V."
flag_F_STA_C: "SIGN-CONFLICT (B<0 radiacja vs B≥0 tensor rdzenia; nieusuwalne; brak okna ekranuje∧tensor-zdrowy)"
implied_aggregate_F_STA_D: "BROKEN (broken = patologia_nieusuwalna ∨ konflikt_znaku = False ∨ True = True; reguła Phase0 §5.1) — formalny LOCK w Phase FINAL"
anti_lakatos_lock: PRESERVED
defers_to_FINAL: "formalny agregat F-STA-D + closure ceremony + propagacja (FOUNDATIONS, rem:B-constraints/O12, parent W-DRR-1, STATE) — osobny user „działaj""
---

# Phase 2 — T-STA-C: zgodność znaku B z rdzeniem (warunek (iv))

## §0 — Verdict at a glance

| Element | Wynik |
|---|---|
| **$c_T^2$ z rdzenia** (S6, reconstr.) | $c_T^2=\dfrac{c_0^2}{1+u_\perp}$, $u_\perp=u\cos^2\theta$, $\theta=\angle(\hat k,\nabla\bar\Phi)$ |
| **Tensor zdrowy** (hiperbol. $\wedge$ subluminal) | $\iff b\ge0$ ($B\ge0$); **subluminalność sama nie wystarcza** (zalicza też $c_T^2<0$) — wymóg denom$>0$ jawnie |
| **$B<0$ ⇒ tensor chory** | nadświetlny ($c_T^2>c_0^2$) dla $\|u\|<1$; **$c_T^2\to\infty$ przy $\|u\|=1$**; $c_T^2<0$ (niestab. gradientowa) dla $\|u\|>1$ |
| **Skalar zdrowy+ekranujący** (Phase 1) | $\iff b<0$ ($B<0$) |
| **Przecięcie znaków** | $\{B<0\}\cap\{B\ge0\}=\varnothing$ ⇒ **SIGN-CONFLICT** |
| **Strefa falowa** | $\hat k\parallel\nabla\bar\Phi$ ⇒ $\cos^2\theta=1$ ⇒ $u_\perp=u$ (ta sama projekcja co skalar) |
| **Próg niestab. tensora** | $\|u\|=1$ = $r_V$ (def. S10); wewnątrz $r_V$ (J0737, $\|u\|\gg1$): $B<0$ ⇒ $c_T^2<0$ |
| **Okno** (ekranuje $\wedge$ tensor zdrowy) | $\varnothing$; ucieczka $\|u\|\ll1$ ⇒ $S\to1$ (ekranowanie trywialne) — nie ratuje |
| **F-STA-C** | **SIGN-CONFLICT** |
| **Implikowany F-STA-D** (§5.1) | **BROKEN** (formalny LOCK: Phase FINAL) |

**sympy: 6/6 PASS.** Werdykt z flag; 0 hardcoded; oba znaki $B$ testowane.

## §1 — Prędkość tensorów i parametr kierunkowy (T-STA-C1)

Z rdzenia `prop:cT` (S6): $c_T^2=A c_0^2/\big(A+b\,q^2\big)$, $q^2=(\hat k\cdot\nabla\bar\Phi)^2\ge0$,
$b=B/M_*^4$. Parametryzując $q^2=X\cos^2\theta$ ($X=(\nabla\bar\Phi)^2$, $\theta=\angle(\hat k,\nabla\bar\Phi)$):
$$\boxed{c_T^2=\frac{c_0^2}{1+u_\perp},\qquad u_\perp=u\cos^2\theta,\qquad u=\frac{bX}{A}.}$$
- $\theta=\tfrac\pi2$ (propagacja **transwersalna**): $u_\perp=0$ ⇒ $c_T=c_0$ (rdzeń case-iii) — ale to **nie** geometria własnej radiacji układu.
- $\theta=0$ (propagacja **radialna**, $\hat k\parallel\nabla\bar\Phi$): $u_\perp=u$ (maksymalna projekcja) — **geometria strefy falowej** (§3).

## §2 — Zdrowy tensor wymaga $B\ge0$ (T-STA-C2)

„Zdrowy" mod tensorowy = **hiperboliczny** (denom $A+bq^2>0$, dodatnia kinetyka) **ORAZ subluminal** ($c_T^2\le c_0^2$).
**Pułapka (jawnie obsłużona):** sama nierówność $c_T^2\le c_0^2$ jest spełniona również w regionie **niestabilnym**
($c_T^2<0\le c_0^2$) — dlatego wymóg denom$>0$ musi wejść jawnie. sympy `solveset` (reprezentanci dodatni, wynik
zależny tylko od znaku $b$, skala-niezmienny — sprawdzone dla dwóch zestawów):
$$\{\text{denom}>0\}\cap\{c_T^2\le c_0^2\}=\{b\ge0\}.$$
Dla **$B<0$** ($b<0$): nadświetlny dla $|u_\perp|<1$ (np. $c_T^2=2c_0^2$), **$c_T^2\to\infty$ przy $|u_\perp|=1$**
(denom$=0$), $c_T^2<0$ (niestabilność gradientowa tensora) dla $|u_\perp|>1$. **Tensor rdzenia zdrowy $\iff B\ge0$.**

## §3 — Geometria strefy falowej: konflikt nieusuwalny (T-STA-C3)

W strefie falowej tło $\bar\Phi$ (dominujący monopol) ma gradient **radialny**; GW z układu propagują **radialnie**
na zewnątrz ⇒ $\hat k\parallel\nabla\bar\Phi$ ⇒ $\cos^2\theta=1$ ⇒
$$u_\perp=u\quad(\text{IDENTYCZNA projekcja, jaką widzi skalar}).$$
Rdzeniowy „case-iii" ($\theta=\tfrac\pi2$, $c_T=c_0$) opisuje GW z **odległego** źródła mijające **inną** masę
poprzecznie — **nie** własną radiację układu. Próg niestabilności tensora $|u_\perp|=1$ pokrywa się **dokładnie**
z definicją $r_V$ (S10: $u(r_V)=1$):
- na zewnątrz $r_V$ ($|u|<1$): $B<0$ ⇒ tensor nadświetlny (skończony);
- **wewnątrz $r_V$ ($|u|>1$, gdzie J0737 leży głęboko)**: denom $=A(1-|u|)<0$ ⇒ $c_T^2<0$ ⇒ **mod tensorowy
  gradient-niestabilny** w samej strefie, gdzie ekranowanie ma działać.

## §4 — SIGN-CONFLICT i brak okna ratunku (T-STA-C4, C5)

$$\text{scalar\_healthy\_sign (Phase 1)}=\{b<0\},\qquad
\text{tensor\_healthy\_sign (§2)}=\{b\ge0\}.$$
$$\boxed{\{b<0\}\cap\{b\ge0\}=\varnothing\ \Rightarrow\ \text{SIGN-CONFLICT.}}$$
Ten sam $B(\Phi)$ (jeden znak w danym tle) nie może być jednocześnie $<0$ (zdrowy skalar) i $\ge0$ (zdrowy tensor).

**Nieusuwalność (T-STA-C5):** czy istnieje *okno* {skalar ekranuje nietrywialnie} $\cap$ {tensor zdrowy} dla $B<0$?
- skalar ekranuje nietryw. ($S=1/(1+|u|)\le\tfrac12$) ⟺ $|u|\ge1$ (wewnątrz $r_V$);
- tensor $B<0$ nie-niestabilny ⟺ $|u|<1$ (na zewnątrz $r_V$);
- przecięcie = $\varnothing$.

**Jedyna ucieczka** — słabe sprzężenie $|u|\ll1$ — daje $S\to1$ (ekranowanie **trywialne**: $\dot P_b$ niezmienione),
więc **nie ratuje** falsyfikowalności sektora; a tensor i tak pozostaje strukturalnie **nadświetlny**. Brak reżimu,
w którym ekranowanie radiacyjne działa **i** tensor jest zdrowy.

## §5 — Co NIE pinuje znaku (T-STA-C6)

- **Polaryzacje** (`prop:disformal-polarization`, S7): mody $h_+,h_\times\propto B\,\partial_i\bar\Phi\,\delta\Phi$
  istnieją dla **$B\ne0$** (oba znaki) ⇒ **nie** preferują znaku. Konflikt **nie** stąd.
- **Statyka** ($r_V$, $\gamma_{\rm PPN}$, S10): używa **tego samego** operatora $Z$ na tle statycznym ($X>0$),
  więc stabilność statycznego tła Vainshteina też preferuje $B<0$ — **spójne** ze skalarem radiacyjnym,
  **nie** dodaje nowego konfliktu.
- **Konflikt znaku pochodzi WYŁĄCZNIE z $c_T$ (sektor tensorowy, `prop:cT`).** Jest to pierwszy warunek (iv),
  który wymusza znak **przeciwny** do skalarowego.

## §6 — Implikowany agregat F-STA-D (formalny LOCK: Phase FINAL)

Reguła Phase 0 §5.1 (IMMUTABLE, pre-zarejestrowana):
```
F-STA-A = HEALTHY        (Phase 1)
F-STA-B = HEALTHY        (Phase 1)
F-STA-C = SIGN-CONFLICT  (ta faza)
patologia_nieusuwalna = (A=GHOST-FORCED) ∨ (B=GRADIENT-FORCED) = False
konflikt_znaku        = (C=SIGN-CONFLICT)                       = True
broken = patologia_nieusuwalna ∨ konflikt_znaku = True
⇒ op-disformal-stability → BROKEN
```
**Sektor radiacyjny TGP_v1 — SFALSYFIKOWANY PRZEZ STABILNOŚĆ.** Disformalny Vainshtein *tłumi* strumień
(parent: PARTIAL), ale jedyny znak $B$ czyniący to tłumienie zdrowym i nietrywialnym ($B<0$) **niszczy
hiperboliczność modu tensorowego** w strefie falowej (gdzie $\hat k\parallel\nabla\bar\Phi$, próg $=r_V$).
To rozstrzyga **W-DRR-1** ostrzej niż strumień: cofa parent UNDERDETERMINED ku falsyfikacji — **przez patologię
operatora, nie przez strumień**.

> **Uwaga zasięgu (forbidden #1):** to NIE rewizja liczb poprzedników. PR-025 (13 227σ/2 646σ), survival
> INDETERMINATE, parent flagi F-DRR-1/2/3 — nietknięte. Zmienia się **agregatowy status sektora radiacyjnego**:
> z UNDERDETERMINED (parent, droga strumienia) na **BROKEN-BY-STABILITY** (ta droga, W-DRR-1 rozstrzygnięte).
> Formalny LOCK agregatu + propagacja = Phase FINAL.

> **Brak sprzeczności z OP-7:** OP-7 orzekło $c_{\rm GW}=c_0$ EXACT w **decoupling na tle KOSMOLOGICZNYM**
> ($\nabla\bar\Phi=0$ ⇒ człon disformalny znika, $c_T=c_0$ dla dowolnego $B$ — rdzeń case-i). Konflikt pojawia
> się dopiero w **lokalnym tle z silnym gradientem** (strefa falowa, $|u|\gtrsim1$) — reżim, którego OP-7 nie
> analizowało dla modu disformalnego. To **nowy** wynik tego cyklu, nie korekta OP-7.

## §7 — DOUBTS / uwagi

- **W-STA-1 (z Phase 1, ROZSTRZYGNIĘTY):** subluminal-skalar / superluminal-tensor dla tego samego $B<0$ —
  potwierdzone jako rdzeń konfliktu; nie współistnieją zdrowo (różny wymagany znak $B$).
- **W-STA-3 (LOW, NOWY):** akceptacja nadświetlności jako „nie-fatalnej" (niektóre EFT, Babichev/Bruneton) mogłaby
  złagodzić $|u|<1$ do CONSISTENT — **ale** reżim ekranowania ($|u|\ge1$, gdzie sektor cokolwiek robi) daje
  $c_T^2<0$ (jawna **niestabilność gradientowa**, nie tylko nadświetlność), więc złagodzenie nie zmienia werdyktu
  w reżimie istotnym. Zapis dla protokołu; nie zmienia F-STA-C.
- **W-STA-2 (z Phase 1, podtrzymany):** magnituda $|B|$/$M_*$ (czy $|u|\gtrsim1$ realnie na J0737) = O12; tu
  rozstrzygamy strukturalnie (znak + geometria), a $r_V(\odot)\gg$ orbita (S10) lokuje J0737 w $|u|\gg1$.

## §8 — Anti-Lakatos compliance (Phase 2)

✓ Werdykt z flag (6/6 PASS); 0 hardcoded T_pass. ✓ **Oba znaki $B$** testowane; tensor-zdrowy ($b\ge0$) i
scalar-zdrowy ($b<0$) wyliczone z `solveset`, przecięcie $\varnothing$ — nie założone. ✓ **Pułapka
„subluminal ⊇ niestabilny" jawnie obsłużona** (wymóg denom$>0$), inaczej fałszywy CONSISTENT. ✓ **Anty-przedwczesny
-negatyw (forbidden #8):** sprawdzono ucieczki (transwersalna geometria case-iii; słabe sprzężenie; akceptacja
nadświetlności W-STA-3) — wszystkie zawodzą w reżimie ekranowania ⇒ BROKEN po jawnym rachunku, nie z analogii.
✓ **Anty-rescue (forbidden #2):** nie strojono $B$ pod CONSISTENT; konflikt zaraportowany wprost wbrew „wygodzie".
✓ $c_T$ z rdzenia EXACT (S6), reused. ✓ Liczby poprzedników nietknięte (forbidden #1); zmiana = agregatowy status.
✓ 0 danych; budżet nowych stałych 0; zakres = sign-pin (forbidden #5).

---

**Phase 2 COMPLETE — F-STA-C = SIGN-CONFLICT. Implikowany agregat F-STA-D = BROKEN.**
Formalny LOCK agregatu + closure ceremony + propagacja (parent W-DRR-1, FOUNDATIONS §3.6.10.6,
`rem:B-constraints`/O12, STATE) + dyspozycja PR = **Phase FINAL**, wymaga osobnego user „działaj".
