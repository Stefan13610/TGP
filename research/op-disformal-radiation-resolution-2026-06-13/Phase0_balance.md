---
title: "Phase 0 — pre-registration LOCK: rozstrzygnięcie D6. Czy disformalny Vainshtein LIVE tłumi STRUMIEŃ energii Ṗ_b z układu podwójnego (nie tylko amplitudę GW dalekiego pola), czy κ_E = ξ_eff/λ (kombinacja kinetyczna C_σσ₀² ⊥ warunkowi amplitudy) da się PRZYPIĄĆ z LIVE, i czy M_* jest WYPROWADZone — wynik domyka sektor radiacyjny TGP_v1: D6 → BROKEN / CLEAN / UNDERDETERMINED."
date: 2026-06-13
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-13 (before any computation)
cycle: op-disformal-radiation-resolution-2026-06-13
parent_cycle: "[[../op-gravitational-sector-survival-2026-06-13/Phase_FINAL_close.md]] CLOSED-RESOLVED INDETERMINATE (D6 = LIVE_UNRESOLVED)"
scope_source: "[[../op-gravitational-sector-survival-2026-06-13/Phase2_derivation.md]] §6 + [[../op-gravitational-sector-survival-2026-06-13/RECON_disformal_vainshtein_2026-06-13.md]] §5"
authorization: "User 2026-06-13: 'zająć się cyklem op-disformal-radiation-resolution' (aktywacja spawn survival §6)"
cycle_category: "MECHANISM-RESOLUTION (rozstrzyga D6 z cyklu-rodzica rachunkiem; wynik = werdykt sektora radiacyjnego)"
parent_verdicts_LOCKED:
  - "PR-025 PSR J0737−3039 TRIGGERED (13 227σ skalar masywny / 2 646σ σ bezmasowy) — branże konkretne falsyfikowane, LOCKED nietykalne"
  - "op-phi-radiative-dof-audit HONEST_NEGATIVE (δΦ JEST radiacyjnym DOF; kinetyka konforemna K₁=const)"
  - "op-gravitational-sector-survival INDETERMINATE (Phase1 8/8 no-go FP7 nad pod-teorią konforemną; Phase2 5/5: D6 disformal LIVE_UNRESOLVED)"
  - "Phase2 Filar I EXACT: L_kin^disformal = A·X − (b/2)X² (człon X² ⇒ premise FP7 złamana) — LOCKED"
  - "Phase2 Filar II EXACT: det J[(λ,ξ)→(λξ, ξ/λ)] = 2ξ/λ ≠ 0 (O_amp ⊥ O_flux ⇒ κ_E unpinned) — LOCKED"
independent_of: "op-nucleation-dimensionality (aktywny; RP² = inwentarz, zakaz mieszania zakresów), PR-022, recovery PR-020"
anti_lakatos_lock: "INHERITED od survival; aktywny od Phase 0"
tags:
  - pre-registration
  - disformal-Vainshtein
  - radiation-energy-balance
  - kappa-E-pinning
  - Mstar-status
  - mechanism-resolution
  - anti-Lakatos-LOCKED
---

# Phase 0 — pre-registration LOCK

> **Status:** Phase 0 = LOCK przed jakimkolwiek rachunkiem. ZERO sympy, ZERO danych
> obserwacyjnych w tym pliku. Każda kolejna faza (1, 2, FINAL) wymaga osobnego user „działaj".

## §1 — Pytanie

Cykl-rodzic ([[../op-gravitational-sector-survival-2026-06-13/Phase_FINAL_close.md]]) orzekł
**INDETERMINATE**: sektor grawitacyjny TGP_v1 nie jest ani sfalsyfikowany, ani uratowany,
ponieważ **D6 — disformalny kanał LIVE — pozostał LIVE_UNRESOLVED** (nie clean: M_* underived
+ κ_E unpinned; nie broken: zachowuje (a) α_i≡0 / (b) statykę 1/r / (c) §1 FOUNDATIONS). Ten
cykl ma **rozstrzygnąć D6 rachunkiem** i tym samym domknąć sektor radiacyjny.

> **Pytanie wiodące:** Czy disformalna struktura LIVE (natywny Vainshtein) sprowadza
> **strumień energii radiacyjnej Ṗ_b** z układu podwójnego do GR — przy M_* rozliczonym
> jako wyprowadzone-lub-postulat i κ_E rozliczonym jako pinowane-lub-swobodne — czy też
> sektor radiacyjny TGP_v1 jest jednak sfalsyfikowany, bądź strukturalnie niefalsyfikowalny?

Trzy pod-pytania (z survival §6 + recon §5), każde z osobną flagą:

1. **W-GSS-1 — strumień vs amplituda (HIGH, najważniejsze):** rdzeń sek07 `[ROZWIĄZANY]`
   + `rem:disformal-status` podają 18-rzędowe tłumienie **amplitudy GW dalekiego pola** z kanału
   disformalnego ($h^{\rm disf}_+\sim h_{\rm GR}/(kr)\sim10^{-40}$). Czy to samo tłumienie
   dotyczy **strumienia energii Ṗ_b** z orbity (bilans Isaacson/$T^{0r}$)? Czy też **sprzężenie
   konforemne $A(\Phi)$ do materii (NIEEKRANOWANE) wciąż wzbudza δΦ i daje $(1/6)P_{\rm GR}$**
   (PR-025), niezależnie od tłumienia amplitudy disformalnej? Vainshtein $r_V(\odot)\sim5\times10^6$ AU
   ≫ separacja J0737 — ale **tłumienie radiacji ≠ ekranowanie statyki**.
2. **W-GSS-3 — pinowanie κ_E (MED→decydujące dla UNDERDETERMINED):** czy istnieje LOCKED relacja
   LIVE pinująca kombinację kinetyczną kontrolującą strumień σ_ab (konkretnie $C_\sigma\sigma_0^2$,
   §2 R8) do wartości GR, **niezależnie** od warunku amplitudy $\xi_{\rm eff}=4\pi G_0\sigma_0\PhiZero$?
   Phase2 Filar II (det J = 2ξ/λ ≠ 0) dowiódł ortogonalności $O_{\rm amp}\perp O_{\rm flux}$ —
   pytanie czy LIVE dostarcza piątego warunku domykającego.
3. **W-GSS-2 — status M_* (MED):** $M_*=m_P$ — derywacja czy postulat? Rozstrzygnięcie
   niespójności klasyfikacji: sek08 tab. (l.171) „Warstwa III — wyprowadzone" vs
   `status_map` (l.214–217) „Propozycja; wymiarowanie $\PhiZero+\ell_P$; **brak mikro-derywacji**".
   Jeśli $M_*$ okaże się **fitowane** (np. do $r_V(\odot)$) ⇒ D6 = tuning-belt (forbidden) ⇒ ku BROKEN.

Cykl **strukturalno-rachunkowy**: zero dostępu do danych obserwacyjnych w fazie wyprowadzania
(K1–K5 cytowane LOCKED, nie przeliczane); dane pulsarowe wyłącznie w ewentualnej fazie
porównawczej, jeśli w ogóle.

## §2 — Wejścia LOCKED (read-only; modyfikacja któregokolwiek = nowy PR, nie rescue)

| # | Input | Wartość / forma | Źródło LOCKED |
|---|---|---|---|
| R1 | Metryka disformalna LIVE | $g_{\mu\nu}=A(\Phi)\eta_{\mu\nu}+\dfrac{B(\Phi)}{M_*^4}\partial_\mu\Phi\,\partial_\nu\Phi$ | sek08 `hyp:disformal` (eq:disformal-metric-formalizm) |
| R2 | Kinetyka disformalna (X-nieliniowa) | $\sqrt{-g}\,g^{\mu\nu}\partial_\mu\Phi\partial_\nu\Phi=A X-\tfrac{b}{2}X^2+\mathcal O(b^2)$, $b=B/M_*^4$ — **EXACT** | Phase2 Filar I (sympy, LOCKED) |
| R3 | Tłumienie amplitudy disformalnej | $h^{\rm disf}_+\sim h_{\rm GR}/(kr)\sim10^{-40}$ (18 rzędów, strefa falowa) — dot. **amplitudy**, nie wykazane dla strumienia | sek07 `[ROZWIĄZANY]`; `rem:disformal-status` |
| R4 | Akcja sektora σ_ab (radiator GW) | $S_\sigma=\int\!\sqrt{-g}\,[\tfrac{C_\sigma}{2}(\partial\sigma)^2-\tfrac{m_\sigma^2}{2}\sigma^2+\tfrac{\xi_{\rm eff}}{\PhiZero^2}\sigma_{ab}\partial^a\Phi\partial^b\Phi]$, $C_\sigma>0$, $m_\sigma\to0$ | sek08 `prop:sigma-eom` (eq:S-sigma) |
| R5 | Równanie ruchu σ_ab | $\Box\sigma_{ab}=S^{\rm TT}_{ab}=-(\xi_{\rm eff}/c_0^4)\Lambda_{ab,cd}T^{cd}$; far-field $\sigma_{ab}=\dfrac{\xi_{\rm eff}}{4\pi c_0^4 r}\ddot Q^{\rm TT}_{ab}$ | sek08 `prop:sigma-eom`, eq:box-sigma, eq:sigma-far-field |
| R6 | Most metryka↔σ | $h^{\rm TT}_{ij}=2\sigma_{ij}/(\sigma_0\PhiZero)$ (sek07/`thm:amplitude-matching`); por. `prop:metric-with-tensor` $2\sigma/\sigma_0$ — **niejednoznaczność czynnika $\PhiZero$ do rozliczenia w Phase 2** | sek07 l.999–1001; sek08 eq:h-TT-from-sigma |
| R7 | Warunek amplitudy (pinuje $\xi_{\rm eff}/\sigma_0$) | $\boxed{\xi_{\rm eff}=4\pi G_0\sigma_0\PhiZero}$ — dopasowuje **amplitudę** $h^{\rm TT}$ do GR | sek08 `thm:amplitude-matching` (eq:xi-matching) |
| R8 | Strumień energii σ (kanoniczny) | Isaacson/$T^{0r}_\sigma\propto C_\sigma\langle\dot\sigma\,\partial_r\sigma\rangle\Rightarrow\dot P_\sigma\propto C_\sigma\sigma_0^2\PhiZero^2\langle\dot h^2\rangle$ — zależy od $C_\sigma\sigma_0^2$, **NIE** od $\xi_{\rm eff}/\sigma_0$ | wyprowadzenie z R4+R6 (do potwierdzenia sympy w Phase 2); Phase2 Filar II struktura |
| R9 | Ortogonalność amp⊥flux | $\det J[(\lambda,\xi)\mapsto(\lambda\xi,\xi/\lambda)]=2\xi/\lambda\neq0$ — **EXACT** ⇒ $O_{\rm amp}\perp O_{\rm flux}$ | Phase2 Filar II (sympy, LOCKED) |
| R10 | Status parametrów σ | $\xi_{\rm eff}$ wyznaczone (R7); $m_\sigma\to0$ przyjęte; **$C_\sigma$ (równoważnie $\sigma_0$) = „dodatkowy parametr, w zasadzie z substratu, obecnie NIEZOBLICZONY"** | sek08 `rem:sigma-params`, `rem:param-counting` |
| R11 | Status $M_*$ | $M_*=m_P\approx1{,}2\times10^{19}$ GeV: `status_map` **„Propozycja; brak mikro-derywacji"**; sek08 tab. „Warstwa III" — **niespójność** | `status_map` l.214–217 vs sek08 l.171 |
| R12 | Skala Vainshteina | $r_V(\odot)\sim5\times10^6$ AU ≫ separacja J0737 (układ głęboko w reżimie ekranowania **statyki**) | recon §3; Phase2 §3 |
| R13 | Bilans radiacyjny konforemny | δΦ radiacyjny (nieekranowany) ⇒ $P_\phi=(1/6)P_{\rm GR}$; gałęzie 13 227σ (σ masywny ⇒ 1/6) / 2 646σ (σ bezmasowy ⇒ 7/6) | PR-025 (LOCKED) |
| R14 | Sprzężenie materia↔δΦ | konforemne $A(\Phi)=e^{2\delta\Phi/\PhiZero}$; źródło $-(q/\PhiZero)\rho\,\delta\Phi$ — **NIEEKRANOWANE pochodnie** (kontakt liniowy w ρ) | Phase0 survival §2; emergent-metric LOCKED |

**Definicja operacyjna κ_E (konkretyzacja Phase2 abstrakcyjnego $\xi/\lambda$):**
$O_{\rm amp}\equiv\xi_{\rm eff}/\sigma_0$ (pinowane przez R7 do $4\pi G_0\PhiZero$);
$O_{\rm flux}\equiv C_\sigma\sigma_0^2$ (kontroluje Ṗ przez R8; NIE pinowane przez R7).
$\kappa_E$ := stosunek $O_{\rm flux}$ do jego wartości GR-Isaacsona. R9 dowodzi $O_{\rm amp}\perp O_{\rm flux}$.
**Zadanie W-GSS-3 = sprawdzić, czy LIVE dostarcza piątego warunku pinującego $O_{\rm flux}$.**

## §3 — Zbiór testów CLOSED (zamknięty w Phase 0; żaden nowy mid-cycle)

Rozstrzygnięcie D6 = trzy testy rachunkowe na PEŁNEJ akcji LIVE (disformal R1–R2 + σ_ab R4).
Zbiór CLOSED — żaden nowy kanał/test nie wchodzi po LOCK. Kanały wkładu do Ṗ_b (recon §3):

| Kanał | Statyka | Radiacja — co rozstrzyga ten cykl |
|---|---|---|
| **C1 — konforemny** (mod oddechowy δΦ) | Newton + PPN ✓ | T1: czy disformal tłumi **strumień**, czy $A(\Phi)$ nieekranowane → $(1/6)P_{\rm GR}$ przeżywa |
| **C2 — disformalny** | dodatkowe polaryzacje | T1: czy 18-rzędowe tłumienie (R3) dotyczy **strumienia energii**, nie tylko amplitudy |
| **C3 — σ_ab (tensor)** | — (TT znika dla sferycznej statyki) | T2: czy $\dot P_\sigma$ = GR (pin $O_{\rm flux}$) czy swobodne (κ_E unpinned) |

| Ti | Test | Obiekt rachunku | A priori zagrożony wynik |
|---|---|---|---|
| **T1** | **Strumień vs amplituda (W-GSS-1).** Bilans energii Isaacson/$T^{0r}$ dla kanału konforemnego C1 (źródło $A(\Phi)$ nieekranowane, R14) ORAZ disformalnego C2 (R2) w reżimie Vainshteina. Czy strumień δΦ z orbity jest tłumiony jak amplituda (R3), czy nieekranowane sprzężenie konforemne wciąż daje $(1/6)P_{\rm GR}$? | $\dot P_\phi^{\rm LIVE}$ (sfera $T^{0r}$ w strefie falowej, akcja disformalna) vs $(1/6)P_{\rm GR}$ | **UNSCREENED_FLUX** ⇒ D6 broken na C1 (werdykt cofa ku FALSIFIED) |
| **T2** | **Pinowanie κ_E (W-GSS-3).** Próba wyprowadzenia $O_{\rm flux}=C_\sigma\sigma_0^2$ z LIVE (zliczanie warunków `rem:param-counting`: 4 warunki na $(\mu,m_0^2,\lambda_0,J)$ — czy domyka piąty na $C_\sigma$?). Czy strumień σ_ab = GR-Isaacson niezależnie od R7? | $\dot P_\sigma/\dot P_{\rm GR}$ jako funkcja $C_\sigma\sigma_0^2$; istnienie LOCKED relacji pinującej | **UNPINNED** ⇒ underdetermined; **PINNED_DEVIATION** ⇒ broken |
| **T3** | **Status M_* (W-GSS-2).** Rozstrzygnięcie niespójności R11. Czy istnieje mikro-derywacja $M_*=m_P$ z substratu (`prop:Mstar-from-substrate`), czy postulat wymiarowy? Czy gdziekolwiek $M_*$ fitowane do $r_V$/danych? | klasyfikacja $M_*$: derived / postulate / fitted | **FITTED** ⇒ tuning-belt forbidden ⇒ broken |

**Uwaga reuse/no-mix:** RP²/nielokalność (D5 cyklu-rodzica) NIE wchodzi — to inwentarz GAP
op-nucleation-dimensionality; zakaz mieszania zakresów (forbidden #9). Ten cykl wyłącznie
o disformal+σ_ab w klasie LIVE lokalnej.

## §4 — Falsyfikatory LOCKED (reguły IMMUTABLE; werdykt WYLICZONY z flag)

| ID | Test | Reguła decyzyjna (IMMUTABLE) |
|---|---|---|
| **F-DRR-1** | T1 — strumień vs amplituda. Sympy: bilans $T^{0r}$ na akcji disformalnej (R1–R2) + kanał konforemny (R14). | Flaga ∈ {`SCREENED_FLUX` (strumień δΦ tłumiony jak amplituda ⇒ człon $1/6$ usunięty z Ṗ_b), `UNSCREENED_FLUX` ($A(\Phi)$ nieekranowane ⇒ $(1/6)P_{\rm GR}$ przeżywa w strumieniu), `PARTIAL` (tłumiony, ale nie do GR)}. Werdykt z **rachunku $T^{0r}$**, nie z analogii amplitudowej. |
| **F-DRR-2** | T2 — pinowanie $O_{\rm flux}=C_\sigma\sigma_0^2$. Sympy: $\dot P_\sigma$ z R4+R6+R8; test czy 4 warunki substratu (R10) domykają piąty. | Flaga ∈ {`PINNED_GR` (LOCKED relacja LIVE wymusza $O_{\rm flux}$ = wartość GR-Isaacson — dowód jawny), `PINNED_DEVIATION` (wymusza ≠ GR), `UNPINNED` (strukturalnie swobodne; det J ≠ 0 R9)}. |
| **F-DRR-3** | T3 — status $M_*$. Audyt `prop:Mstar-from-substrate` + grep za fitem $M_*↔r_V$. | Flaga ∈ {`DERIVED` (mikro-derywacja z pierwszych zasad), `POSTULATE` (wymiarowy, underived, NIE fit), `FITTED` (strojony do $r_V$/danych)}. |
| **F-DRR-C** | Agregat (WYLICZONY z F-DRR-1 ∧ -2 ∧ -3) | Reguła w §4.1. Próg liczbowy NIE dotyczy (klasyfikacyjne); zero hardcoded T_pass. |

### §4.1 — Reguła agregatu F-DRR-C (IMMUTABLE)

```
broken := (F-DRR-1 = UNSCREENED_FLUX) ∨ (F-DRR-2 = PINNED_DEVIATION) ∨ (F-DRR-3 = FITTED)
clean  := (F-DRR-1 = SCREENED_FLUX) ∧ (F-DRR-2 = PINNED_GR) ∧ (F-DRR-3 = DERIVED)

if broken:                              ⇒ D6 → BROKEN
elif clean:                             ⇒ D6 → CLEAN
elif (F-DRR-2 = UNPINNED) ∧ ¬broken:    ⇒ D6 → UNDERDETERMINED
else:                                   ⇒ D6 → UNDERDETERMINED + R1 flag (luka; ewentualny spawn)
```

- **D6 → BROKEN** (HONEST_NEGATIVE = pełnoprawny PASS): disformal NIE tłumi Ṗ_b **lub** κ_E
  daje deviację **lub** $M_*$ fitowane ⇒ sektor radiacyjny jednak sfalsyfikowany; werdykt
  sektora cofa się ku **FALSIFIED-AS-AXIOMATIZED** (specyfikacja warunków brzegowych dla v2).
- **D6 → CLEAN** (wymaga 3/3 warunków clean — 2/3 NIE wystarcza): tłumienie strumienia
  WYPROWADZONE + $M_*$ derived + κ_E → GR bez strojenia ⇒ ratunek sektora; kandydat **v1.1**;
  nowy PR z testowalną predykcją.
- **D6 → UNDERDETERMINED:** $O_{\rm flux}$ strukturalnie swobodne ⇒ sektor radiacyjny **nie
  czyni falsyfikowalnej predykcji w obecnej formie** (osobny status: underdetermination-parametryczna,
  NIE INDETERMINATE-mechanizm); uczciwy wynik wymagający domknięcia teoretycznego (pin $C_\sigma$ z substratu).

## §5 — Pre-derywacja (oczekiwanie, NIE próg)

Oczekiwane a priori (do obalenia rachunkiem — rachunek nadrzędny):
- **T1:** ryzyko `UNSCREENED_FLUX` realne — sprzężenie konforemne $A(\Phi)$ do materii (R14) jest
  liniowe w ρ i **nie ma pochodnych do ekranowania Vainshteinem**; tłumienie 18-rzędów (R3) wykazane
  dla *amplitudy h dalekiego pola*, nie dla *strumienia skalarnego z orbity* (recon §4 explicite).
  To najpoważniejsze zagrożenie dla rescue — i dlatego T1 = Phase 1 (priorytet).
- **T2:** ryzyko `UNPINNED` wysokie — `rem:sigma-params` mówi wprost, że $C_\sigma$ jest
  „dodatkowym parametrem, obecnie niezobliczonym"; R9 (det J ≠ 0) dowiódł ortogonalności.
- **T3:** oczekiwane `POSTULATE` (status_map), NIE `FITTED` (recon §2: „nie fit $r_V$") —
  ale niespójność sek08 vs status_map musi być rozstrzygnięta jawnie, nie wygładzona.

**Jeśli rachunek da wynik przeciwny oczekiwaniu — raportować wprost.** `SCREENED_FLUX`+`PINNED_GR`+`DERIVED`
(rescue) byłby silnym wynikiem pozytywnym programu; `UNSCREENED_FLUX` (falsyfikacja) — pełnoprawnym negatywem.

## §6 — Forbidden moves (IMMUTABLE; INHERITED z survival + cyklowe)

1. Zakaz rewizji werdyktów LOCKED: PR-025 liczby (13 227σ/2 646σ), survival INDETERMINATE,
   Phase1 no-go FP7, Phase2 Filary I/II EXACT, radiative-dof-audit. Zmiana = osobny PR.
2. **Zakaz strojenia κ_E ($O_{\rm flux}=C_\sigma\sigma_0^2$) do GR/danych** (analog κ_E PR-025 §6.2);
   pin musi być WYPROWADZONY z LIVE lub raportowany jako UNPINNED.
3. **Zakaz strojenia $M_*$** do $r_V$/danych; $M_*=m_P$ traktowane jako **underived** dopóki
   mikro-derywacja nie pokazana (status quo status_map, R11).
4. **Bilans ENERGII (Isaacson/$T^{0r}$), NIE amplitudy** — lekcja PR-025 T5/recon §4: „$h^\sigma_{TT}=h^{GR}_{TT}\Rightarrow$ GR w Ṗ_b" jest *non sequitur* (R9). Każdy werdykt o Ṗ_b z jawnego strumienia.
5. Zakaz dostępu do danych obserwacyjnych w fazie wyprowadzania; dane pulsarowe wyłącznie
   w ewentualnej fazie porównawczej, oznaczonej.
6. Zakaz nowych stałych (budżet 0; $A,B,M_*,C_\sigma,\sigma_0,\xi_{\rm eff},\PhiZero,\lambda,q,m$ symboliczne, z rdzenia).
7. Zakaz hardcoded T_pass; flagi F-DRR-1/2/3/C wyliczane z jawnego rachunku (lekcja GST #19).
8. Symetryczna dyscyplina: zakaz rescue-by-tuning (#2/#3) **ORAZ** zakaz przedwczesnego
   negatywu — BROKEN dopuszczalny tylko po jawnym rachunku $T^{0r}$ (nie z analogii amplitudowej).
9. Zakaz mieszania zakresu z op-nucleation-dimensionality (RP²/D5 = inwentarz GAP tam, nie tu).
10. Zakaz miękkiego domknięcia: BROKEN / CLEAN / UNDERDETERMINED — wszystkie trzy pełnoprawne;
    żadnego łagodzenia ku „prawie uratowane".
11. `D6 → CLEAN` wymaga dowodu 3/3 (SCREENED_FLUX ∧ PINNED_GR ∧ DERIVED); 2/3 = non-clean.
12. Zakaz altitude-creep do budowy v2; ten cykl = werdykt D6 + (jeśli BROKEN) warunki brzegowe;
    v2 = osobny program.

## §7 — Risk register

- **R-flux-amplitude-conflation (HIGH):** centralna pułapka — założenie „tłumienie amplitudy (R3)
  ⇒ tłumienie strumienia". Sprzężenie konforemne (R14) jest nieekranowane. Mitygacja: jawny
  $T^{0r}$ kanał-po-kanale (C1/C2/C3), forbidden #4.
- **R-rescue-temptation (HIGH):** po INDETERMINATE rodzica pokusa ogłoszenia CLEAN. Mitygacja:
  F-DRR mechaniczne z flag; pin κ_E wymaga LOCKED relacji, nie wiarygodności; wymóg 3/3 (#11).
- **R-premature-negative (MED, symetryczne):** łatwy BROKEN bez pełnego rachunku strumienia.
  Mitygacja: forbidden #8; T1 musi policzyć $T^{0r}$, nie tylko zacytować R3.
- **R-Mstar-inconsistency (MED):** sek08 „Warstwa III" vs status_map „Propozycja". Musi być
  rozstrzygnięte jawnie (T3), nie wygładzone; jeśli nierozstrzygalne → POSTULATE (konserwatywnie).
- **R-Vainshtein-regime (MED):** $r_V(\odot)\sim5\times10^6$ AU ≫ orbita (R12) dotyczy ekranowania
  **statyki**; strefa radiacji jest na zewnątrz. Mitygacja: jawne rozróżnienie gdzie liczony strumień
  (bliskie pole w $r_V$ vs strefa falowa); nie utożsamiać „głęboko w Vainshteinie statycznie" z „radiacja wytłumiona".
- **R-kappa-E-grounding (MED):** abstrakcyjne $(\lambda,\xi)$ Phase2 muszą być ugruntowane w
  konkretach $C_\sigma,\sigma_0,\xi_{\rm eff}$ (R8/§2 def). Mitygacja: Phase 2 zaczyna od jawnej
  identyfikacji $O_{\rm flux}=C_\sigma\sigma_0^2$ z eq:S-sigma; rozliczenie niejednoznaczności $\PhiZero$ (R6).

## §8 — Anti-Lakatos checklist (Phase 0)

- ✓ Phase 0 LOCKED przed jakimkolwiek rachunkiem (zero sympy w tym pliku).
- ✓ Zbiór testów {T1,T2,T3} CLOSED; kanały {C1,C2,C3} zamknięte; żaden nowy mid-cycle.
- ✓ Werdykty F-DRR-1/2/3/C wyliczane z flag; 0 hardcoded T_pass.
- ✓ Trzy wyniki (BROKEN/CLEAN/UNDERDETERMINED) pre-deklarowane jako pełnoprawne.
- ✓ Pre-derywacja §5 = oczekiwanie, nie próg.
- ✓ Werdykty poprzedników LOCKED i nietykalne (forbidden #1); korekta zakresu survival respektowana.
- ✓ Budżet nowych stałych: 0.
- ✓ Symetryczna ochrona: anty-rescue-by-tuning (#2/#3/#11) ORAZ anty-przedwczesny-negatyw (#8).
- ✓ Bilans energii (nie amplitudy) zafiksowany jako reguła (#4) — lekcja PR-025 T5 wbudowana.
- ✓ Niespójność M_* (sek08 vs status_map) zarejestrowana jako do-rozstrzygnięcia (R11/T3), nie wygładzona.

## §9 — Plan faz (każda wymaga osobnego user „działaj")

- **Phase 0** — ten plik (LOCK). ✅
- **Phase 1 — T1 (strumień vs amplituda; F-DRR-1):** sympy — bilans energii Isaacson/$T^{0r}$
  na akcji disformalnej (R1–R2) dla kanału konforemnego C1 (źródło $A(\Phi)$ nieekranowane) i disformalnego C2;
  rozstrzygnięcie czy 18-rzędowe tłumienie (R3) dotyczy strumienia czy tylko amplitudy. Output: flaga F-DRR-1.
- **Phase 2 — T2 + T3 (κ_E + M_*; F-DRR-2, F-DRR-3):** sympy — $\dot P_\sigma$ z R4/R6/R8,
  identyfikacja $O_{\rm flux}=C_\sigma\sigma_0^2$, test pinowania (zliczanie warunków substratu);
  audyt `prop:Mstar-from-substrate` + grep za fitem $M_*↔r_V$. Output: flagi F-DRR-2, F-DRR-3.
- **Phase FINAL** — agregat F-DRR-C (z flag); closure; propagacja (FOUNDATIONS §3.6.10.6 CL-2
  domknięcie/aktualizacja, REALITY_CONTACT_AUDIT, STATE); dyspozycja PR (numer/adnotacja — wyłącznie user).

## §10 — Numer PR

RESERVED — przydzieli Phase FINAL. Scenariusze: **BROKEN** → adnotacja FALSIFIED-AS-AXIOMATIZED
sektora radiacyjnego (bez nowego PR obserwacyjnego, lub adnotacja PR-025 „scope domknięty");
**CLEAN** → nowy PR (kandydat PR-026) z predykcją sektora σ_ab v1.1; **UNDERDETERMINED** →
adnotacja strukturalna „sektor radiacyjny niefalsyfikowalny w obecnej formie; wymaga pinowania $C_\sigma$".
