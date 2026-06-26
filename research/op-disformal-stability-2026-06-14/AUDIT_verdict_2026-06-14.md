---
title: "AUDYT WERDYKTU (reviewer) — op-disformal-stability: argument Phase 2 (BROKEN via induced-TT) jest NIEROBUSTNY, ale KONKLUZJA (BROKEN) jest prawdopodobnie POPRAWNA via inny, solidny mechanizm: disformal viability (g_eff traci Lorentzowskość przy |u|=1=r_V). Werdyktu NIE lockować na obecnym argumencie; re-derywacja w cyklu Hamiltonowskim."
date: 2026-06-14
type: reviewer-audit
status: "AUDIT (zewnętrzny przegląd; NIE faza cyklu; NIE rewiduje liczb poprzedników; rekomenduje wstrzymanie Phase FINAL na obecnym argumencie)"
cycle_audited: op-disformal-stability-2026-06-14
auditor: "Claudian, sesja #26 (user: 'sprawdź czy wyniki poprawne / czy coś umknęło' → 'C')"
sympy: "[[./AUDIT_verdict_sympy.py]] / [[./AUDIT_verdict_sympy.txt]] — 3 rachunki EXACT"
tags:
  - reviewer-audit
  - verdict-dispute
  - disformal-viability
  - signature-flip
  - anti-Lakatos
---

# AUDYT WERDYKTU — op-disformal-stability (BROKEN)

## §0 — Konkluzja audytu w jednym akapicie

Phase 1 (B<0 = jedyny zdrowy+ekranujący znak dla **fizycznego skalara**) jest **poprawna**.
Phase 2 (BROKEN przez SIGN-CONFLICT na $c_T$ z `prop:cT`) jest **nierobustna**: opiera się na
trybie, który rdzeń **sam** oznacza jako niefizyczny (`rem:GW-scope-2026`), policzonym **naiwnym
operatorem** (≠ właściwy operator fluktuacji k-essence), ekstrapolowanym poza WKB. **Jednak
końcowa KONKLUZJA (BROKEN) jest prawdopodobnie POPRAWNA** — z innego, solidnego powodu, który
Phase 2 przeoczyła: **disformal viability**. Efektywna metryka $g_{\rm eff}$ traci sygnaturę
Lorentzowską dokładnie przy $|u|=1$ (=$r_V$) dla $B<0$, a reżim silnego ekranowania ($|u|\gg1$)
jest albo niefizyczny geometrycznie ($B<0$), albo ghostowy ($B>0$). **Werdyktu NIE wolno
zalockować na argumencie Phase 2; należy go re-wyprowadzić via viability (cykl Hamiltonowski).**

## §1 — Co jest poprawne (Phase 1) ✅

Fizyczny skalar δΦ na statycznym tle disformalnym (operator k-essence $Z^{\mu\nu}$, EXACT):
dla **B<0** mamy no-ghost ($L'>0$), $c_s^2>0$, $Z^{rr},Z^\perp>0$, ekranowanie $S<1$ — zdrowy
i ekranujący dla całego $u<0$. **B<0 jako jedyny zdrowy znak skalara — potwierdzam.**

## §2 — Dlaczego argument Phase 2 jest nierobustny ❌

Phase 2 §5 mówi wprost: *„Konflikt znaku pochodzi WYŁĄCZNIE z $c_T$ (`prop:cT`)."* A z tym $c_T$:

1. **`rem:GW-scope-2026` (rdzeń, tuż pod `prop:cT`, l.6606–6613):** δΦ^TT to *„efektywna część
   bezśladowa indukowanej perturbacji metryki, a NIE niezależne tensorowe stopnie swobody — to
   nie jest sektor spin-2 w sensie OTW."* ⇒ jego $c_T^2<0$ **nie jest** niestabilnością gradientową
   propagującego modu. To obiekt **slaved** do δΦ (cała δg jest funkcją δΦ; metryka emergentna).
2. **Naiwny operator:** `prop:cT` używa $A\Box + b\,\partial\bar\Phi\partial\bar\Phi\,\partial\partial$
   (d'Alembertian metryki na perturbacji), a **właściwy** operator fizycznego DOF to k-essence
   $Z^{\mu\nu}=2(A-bX)\eta^{\mu\nu}-4b\,\partial\bar\Phi\partial\bar\Phi$. **Rachunek (AUDIT_sympy §1):**
   fizyczna dyspersja skalara $c_{\rm skalar}^2/c_0^2=(1-3u)/(1-u)$ — dla B<0 (u=−W): $(1+3W)/(1+W)>0$
   **zawsze (zdrowy)**. Czyli „niestabilność", na której stoi BROKEN, **nie istnieje** w fizycznym modzie.
3. **WKB poza zakresem:** `prop:cT` jest WKB strefy falowej ($k\gg1/r_{\rm src}$); Phase 2 ekstrapoluje
   ją do strefy bliskiej $|u|>1$ (sama to przyznaje).
4. **Niespójność rdzenia:** boxed `eq:cT` ma projekcję $\perp$ w mianowniku, dowód (l.6586) projekcję
   $\parallel$ w liczniku — Phase 2 wybrała jedną bez flagowania sprzeczności.

**Wniosek §2:** SIGN-CONFLICT „zdrowy skalar B<0 vs zdrowy tensor B≥0" jest **pozorny** — „tensor"
to niefizyczny slaved-TT, a fizyczny skalar jest zdrowy dla B<0.

## §3 — Dlaczego KONKLUZJA (BROKEN) jest jednak prawdopodobnie POPRAWNA ✅ (właściwy argument)

To, co Phase 2 przypisała „niestabilności tensora przy $|u|=1$", jest **naprawdę** degeneracją
efektywnej metryki. **Rachunek (AUDIT_sympy §3, EXACT):**
$$g_{\rm eff}=A\eta+b\,\partial\bar\Phi\partial\bar\Phi=\mathrm{diag}(-A,\;A+bG^2,\;A,\;A),\qquad
\det g_{\rm eff}=-A^4(1+u).$$
Wartość własna radialna $A(1+u)$ zmienia znak przy $u=-1$:
- **B>0** (u>0): $1+u>0$ zawsze ⇒ $g_{\rm eff}$ Lorentzowska dla każdego u.
- **B<0** (u=−W): $A(1-W)$ — przy $W=1$ **degeneracja** ($\det=0$), przy $W>1$ **flip sygnatury**
  (kierunek radialny staje się czasowy ⇒ dwie osie czasowe ⇒ utrata Lorentzowskości/przyczynowości).

To jest **standardowy warunek disformal viability**: $1+(B/A)X>0$. Próg $|u|=1$ **pokrywa się
dokładnie** z $r_V$ — to nie przypadek, to ta sama degeneracja, którą Phase 2 zinterpretowała błędnie.

**Trylemat (AUDIT_sympy WNIOSEK):**
| reżim | $g_{\rm eff}$ | skalar | screening $S=1/|1-u|$ |
|---|---|---|---|
| B<0, $\|u\|>1$ | **flip sygnatury** ✗ | zdrowy | silny ✓ |
| B>0, $\|u\|>2$ | Lorentz ✓ | **ghost** ✗ | silny ✓ |
| B<0, $\|u\|<1$ | Lorentz ✓ | zdrowy | **$S\!\to\!1$ brak** ✗ → PR-025 konforemny stoi |

$$\boxed{(g_{\rm eff}\ \text{Lorentz})\ \wedge\ (\text{skalar zdrowy})\ \wedge\ (\text{silne ekranowanie})=\varnothing.}$$

Żaden znak B nie daje jednocześnie zdrowej geometrii, zdrowego skalara i nietrywialnego ekranowania.
**Sektor radiacyjny jest sfalsyfikowany** — ale przez **niekompatybilność viability↔ekranowanie**, nie
przez „niestabilność tensora". Agent **trafił w odpowiedź, chybił w dowodzie.**

## §4 — Korekta mojej własnej poprzedniej oceny (uczciwość)

W poprzedniej turze oceniłem „BROKEN over-claimed, prawdopodobnie SIGN-PINNED B<0". To było
**przedwczesne**: oparte na obserwacji zdrowia skalara (poprawnej), ale **przed** sprawdzeniem
viability $g_{\rm eff}$. Pełny rachunek determinantu pokazuje, że reżim ratujący (B<0, $|u|\gg1$)
jest geometrycznie niedopuszczalny. Korekta: **nie SIGN-PINNED — najpewniej BROKEN, lecz via
inny mechanizm niż podał agent.** Rachunek nadrzędny nad pierwszą intuicją.

## §5 — Status i residualna subtelność (do cyklu Hamiltonowskiego)

- **NIE lockować Phase FINAL na argumencie Phase 2** (induced-TT) — jest nierobustny.
- Werdykt **prawdopodobnie BROKEN via viability**, ale wymaga formalnego potwierdzenia:
  1. **Hyperboliczność/sygnatura $g_{\rm eff}$** w reżimie ekranowania (czy flip jest nieusuwalny przez
     dopuszczalne $B(\Phi)$, znak, O12) — rachunek (AUDIT §3) wskazuje TAK, ale na pełnej dynamicznej analizie.
  2. **Czy silne ekranowanie realnie wymaga $|u|>1$** (skaling $S=1/|1-u|$ to heurystyka z
     op-disformal-radiation-resolution — jeśli pełny DRW da inny skaling, trylemat trzeba przeliczyć).
  3. **Magnituda B/M_* (O12):** czy $|u|>1$ jest osiągane fizycznie — ale **kluczowe:** trylemat jest
     **niezależny od O12** (dla DOWOLNEGO B: albo |u|<1 → brak screeningu → PR-025, albo |u|>1 →
     patologia). To czyni BROKEN robustnym względem nieznanego B — silniejsze niż UNDERDETERMINED rodzica.
- Residual (gdyby viability dało się ocalić): skalar B<0 jest **nadświetlny** ($c^2\to3c_0^2$) —
  osobna kwestia przyczynowości w teorii z preferowaną klatką (drugorzędna wobec sygnatury).

## §6 — Anti-Lakatos (audyt)

✓ Liczby poprzedników nietknięte (PR-025, survival, parent — LOCKED). ✓ Werdykt agenta NIE odrzucony
pochopnie: Phase 1 potwierdzona, błąd Phase 2 wskazany rachunkiem, konkluzja **wzmocniona** lepszym
argumentem (nie osłabiona). ✓ Własna poprzednia ocena skorygowana jawnie (§4). ✓ Wszystkie 3 rachunki
EXACT, reprodukowalne (AUDIT_verdict_sympy). ✓ Residualne subtelności (skaling screeningu, O12,
nadświetlność) wypisane jako otwarte, nie wygładzone. ✓ 0 danych; 0 nowych stałych.

## §7 — Rekomendacja

Cykl **`op-disformal-hamiltonian-viability`** ([[../op-disformal-hamiltonian-viability-2026-06-14/]],
scoping [[../../meta/SCOPING_op-disformal-hamiltonian-viability_2026-06-14.md]]): formalnie
rozstrzygnąć BROKEN via viability (sygnatura $g_{\rm eff}$ + DOF count potwierdzający slaved-TT +
niezależność od O12 + skaling screeningu). Phase FINAL op-disformal-stability **dopiero po** tym —
z poprawnym argumentem, NIE z induced-TT.
