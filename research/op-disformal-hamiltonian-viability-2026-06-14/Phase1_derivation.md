---
title: "Phase 1 — F-VIA-A..E: sektor radiacyjny TGP_v1 → BROKEN-via-viability. g_eff traci sygnaturę Lorentzowską przy |u|=1 (=r_V) dla B<0 (EXACT); induced-TT SLAVED (argument Phase2 op-disformal-stability formalnie void); trylemat {Lorentz}∩{skalar zdrowy}∩{screening}=∅ dla obu znaków B (O12-niezależnie). 5/5 PASS."
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-disformal-hamiltonian-viability-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14: 'działaj'"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "5/5 PASS; 0 hardcoded; werdykt z flag; EXACT: det g_eff=-A^3(A+bG^2)=-A^4(1+u); trylemat=EmptySet (sympy solveset/Intersection)."
flag_F_VIA_A: "SIGNATURE-FLIP (B<0,|u|>1: radial eig A(1+u)<0)"
flag_F_VIA_B: "SLAVED (δg∝δΦ; induced-TT nie-DOF; argument induced-TT Phase2 void)"
flag_F_VIA_C: "EMPTY (trylemat pusty B<0 i B>0)"
flag_F_VIA_D: "SCREENING-NEEDS-|u|>~1 (S=1/|1-u|<1/2 ⟺ |u|>1)"
flag_F_VIA_E: "BROKEN-via-viability"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — viability: sektor radiacyjny → BROKEN

## §0 — Verdict at a glance
| Falsyfikator | Wynik |
|---|---|
| **F-VIA-A** sygnatura $g_{\rm eff}$ | **SIGNATURE-FLIP** — $g_{\rm eff}=\mathrm{diag}(-A,A(1+u),A,A)$; radial eig $<0$ dla $\|u\|>1$ (B<0) ⇒ flip (drugi czas) |
| **F-VIA-B** DOF count | **SLAVED** — δg∝δΦ (metryka emergentna); induced-TT nie ma własnego kinetyka ⇒ nie-DOF ⇒ argument induced-TT (op-disformal-stability Phase 2) **formalnie void** |
| **F-VIA-D** skaling screeningu | **NEEDS-\|u\|>~1** — $S=1/\|1-u\|<\tfrac12\iff\|u\|>1$ |
| **F-VIA-C** trylemat | **EMPTY** — {Lorentz}∩{skalar zdrowy}∩{screening}=∅ dla **obu** znaków B |
| **F-VIA-E** agregat | **BROKEN-via-viability** |

**5/5 PASS. Werdykt z flag; induced-TT NIE użyte jako dowód (F-VIA-B je zamyka).**

## §1 — Twarda geometria (EXACT)
$$g_{\rm eff}=A\eta+b\,\partial\bar\Phi\partial\bar\Phi=\mathrm{diag}(-A,\;A+bG^2,\;A,\;A),\quad
\det g_{\rm eff}=-A^4(1+u),\ u=\tfrac{bX}{A}.$$
Wartość własna radialna $A(1+u)$ zmienia znak przy $u=-1$. **B<0, $\|u\|>1$:** radialna $<0$ ⇒
dwie osie czasowe ⇒ utrata Lorentzowskości (standardowy warunek disformal viability $1+(B/A)X>0$).
**B>0:** $1+u>0$ zawsze (Lorentz OK). Próg $\|u\|=1$ **pokrywa się z $r_V$** — to degeneracja $g_{\rm eff}$,
NIE „niestabilność tensora" (korekta op-disformal-stability Phase 2).

## §2 — induced-TT jest SLAVED (zamknięcie argumentu Phase 2)
Metryka emergentna $g_{\rm eff}=g_{\rm eff}(\Phi,\partial\Phi)$; perturbacja (eq:disformal-perturbation)
$\delta g_{\mu\nu}=A'\delta\Phi\,\eta+\tfrac{B}{M_*^4}[\partial\delta\Phi\,\partial\bar\Phi+\dots]$ —
**wszystkie człony liniowe w δΦ/∂δΦ; brak niezależnego członu kinetycznego dla $h^{TT}$.** Więc
$\delta g^{TT}$ jest niedynamiczne (slaved do δΦ), a jego formalne $c_T^2<0$ (`prop:cT`) **nie jest**
niestabilnością propagującego modu (`rem:GW-scope-2026`: „NIE niezależne tensorowe stopnie swobody").
Fizyczny skalar jest zdrowy dla B<0 ($(1-3u)/(1-u)>0$). **Argument BROKEN-via-induced-TT formalnie void.**

## §3 — Trylemat (rdzeń werdyktu; O12-niezależny)
Zbiory (sympy, EXACT): zdrowy skalar $=\{u<\tfrac13\}$; Lorentz $g_{\rm eff}=\{u>-1\}$;
nietrywialny screening $=\{|u|>1\}$ ($u<-1$ ∨ $u>3$).
$$\{u<\tfrac13\}\cap\{u>-1\}\cap(\{u<-1\}\cup\{u>3\})=\varnothing\quad(\text{B<0}:\varnothing,\ \text{B>0}:\varnothing).$$
- **B<0:** skalar zdrowy + screening współistnieją ($u<-1$), **ale tam $g_{\rm eff}$ flip** ⇒ wykluczone.
- **B>0:** $g_{\rm eff}$ OK, **ale screening ($u>3$) = ghost skalara** ($L'=A(1-u)<0$).
- **B<0, $\|u\|<1$:** $g_{\rm eff}$ OK + skalar zdrowy, **ale $S\to1$ — brak ekranowania** ⇒ PR-025 konforemny stoi.

**Dla DOWOLNEGO B (znak, magnituda): albo $\|u\|<1$ (brak ekranowania → PR-025), albo $\|u\|>1$
(patologia geometryczna lub ghost).** ⇒ niezależne od O12 ($B(\Phi)$, $M_*$ nieznane — werdykt i tak stoi).

## §4 — Werdykt F-VIA-E i jego status epistemiczny
> **D6 / sektor radiacyjny TGP_v1 → BROKEN-via-viability.** Disformalny mechanizm ratunkowy
> NIE może ekranować $\dot P_b$ w zdrowy sposób: silne ekranowanie wymaga dominującego członu
> disformalnego ($\|u\|\gtrsim1$), a ten albo łamie Lorentzowskość $g_{\rm eff}$ (B<0), albo daje
> ghost skalara (B>0). Konforemna falsyfikacja PR-025 **rozciąga się na pełne LIVE**.

**Twardość:** geometria (F-VIA-A flip, EXACT; ghost z $Z$ LOCKED) — twarda. **Jedyna nie-twarda
przesłanka:** „nietrywialne ekranowanie ⇒ $\|u\|\gtrsim1$" (jakościowo solidne: perturbacyjnie mały
disformal nie da O(1) supresji; specyficzny skaling $1/|1-u|$ dziedziczony jako heurystyka z
op-disformal-radiation-resolution — F-VIA-D). Gdyby pełny rachunek DRW dał silne ekranowanie przy
$\|u\|<1$, trylemat trzeba przeliczyć (mało prawdopodobne strukturalnie).

## §5 — Co to znaczy dla łańcucha
- **op-disformal-stability:** werdykt BROKEN **potwierdzony, ale POPRAWNYM argumentem** (viability),
  nie induced-TT. Jego Phase FINAL może domknąć się cytując ten cykl.
- **op-disformal-radiation-resolution (UNDERDETERMINED):** zaostrzony do BROKEN — bo viability jest
  O12-niezależne (UNDERDETERMINED czekało na O12; tu O12 nie jest potrzebne).
- **Sektor grawitacyjny TGP_v1:** konforemny FALSIFIED (PR-001/004/025) + disformalny BROKEN-via-viability
  ⇒ **cały sektor radiacyjny/dalekozasięgowy sfalsyfikowany, z czystym dowodem strukturalnym.**

## §6 — Anti-Lakatos (Phase 1)
✓ 5/5 PASS; werdykt z flag; 0 hardcoded. ✓ induced-TT NIE użyte jako dowód (F-VIA-B formalnie void) —
naprawiono błąd op-disformal-stability Phase 2. ✓ Geometria EXACT; element heurystyczny (skaling
screeningu) jawnie oznaczony jako jedyna nie-twarda przesłanka (F-VIA-D, §4). ✓ O12-niezależność jawna.
✓ Liczby poprzedników nietknięte. ✓ Budżet stałych 0. ✓ Werdykt nie wymuszony — NOT-BROKEN był
pre-deklarowany (Phase 0 §3) i odrzucony rachunkiem (trylemat EMPTY).

**Następny krok (user „działaj"): Phase FINAL — LOCK F-VIA-E; domknięcie op-disformal-stability
poprawnym argumentem; propagacja (FOUNDATIONS CL-6: sektor radiacyjny UNDERDETERMINED→FALSIFIED-via-viability;
REALITY_CONTACT_AUDIT; STATE; dyspozycja PR). To terminalny werdykt sektora grawitacyjnego — wymaga
świadomej autoryzacji.**
