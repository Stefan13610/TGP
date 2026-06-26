---
title: "Phase 4 — moduł B: KB4 NEGATIVE_FOR_REAL_WALL (H-CP wykluczone w LIVE machinery — akcja dokładnie parzysta w χ) · SIG-2 BOUNDED (f_leak ≲ 1.2×10⁻³; A-i strukturalnie domyka kanał; antymateria trwała → 0; obserwabla radiacyjna → kandydat PR-023)"
type: phase_result
status: PHASE4_COMPLETE
phase: 4
cycle: op-frontier-bridge-and-asymmetry-2026-06-12
created_date: 2026-06-12
authorization: "User 2026-06-12: 'Phase 4'"
sympy_script: "[[./Phase4_sympy.py]]"
sympy_output: "[[./Phase4_sympy.txt]]"
sympy: "9/9 PASS; 0 hardcoded T_pass; circularity guard FP9 (anchor wyłącznie w linii porównawczej FP7)"
falsifier_status: "KB4 = NEGATIVE_FOR_REAL_WALL (DERIVED, wszystkie rzędy) + GAP deklarowany (textured wall / RP² holonomia) · SIG-2 = BOUNDED (górna granica wyprowadzona pod zadeklarowanym modelem wag; porównanie literal-anchor PASS) · F-BA-6 klasyfikacja → Phase FINAL"
anti_lakatos_lock: PRESERVED
---

# Phase 4 — moduł B: amplituda H-CP + leakage SIG-2

## §0 — Verdict at a glance

| Element | Wynik | Esencja |
|---|---|---|
| **KB4 (H-CP)** | **NEGATIVE_FOR_REAL_WALL (DERIVED, wszystkie rzędy)** | Akcja TGP rozwinięta wokół REALNEGO profilu ściany jest **dokładnie parzysta w χ** (C: Φ→Φ* ⟺ χ→−χ; \|Φ\|² zawiera tylko χ²) ⇒ amplituda kreacji NIE rozróżnia Φ/Φ* w żadnym rzędzie perturbacyjnym ⇒ **H-CP wykluczone w LIVE machinery**. GAP deklarowany: ściana z teksturą fazową / sektory holonomii RP² |
| **SIG-2 (leakage)** | **BOUNDED** | Kanał leakage otwiera się dopiero na głębokości ξ > ln 72 ≈ 4.28 ≫ locus kreacji 1.32; **A-i (LOCKED) ogranicza kreację do warstwy ⇒ strukturalne domknięcie kanału**; górna granica f_leak ≲ 1.2×10⁻³ (pas konwencji 1.3×10⁻⁴–1.2×10⁻³); wyciekła antymateria = TRANSIENT (anihilacja par w bulku) ⇒ frakcja trwała → 0 — **comparison-only: PASS vs f̄_max = 10⁻⁶**; obserwabla wtryskowa f_rad ≈ 2f_leak → kandydat PR-023 |
| **F-BA-6** | OPEN → Phase FINAL | Zbiór hipotez zawęża się: H-CP wykluczone (real wall); HONest_NEGATIVE wykluczone (KB1 ≠ null); pozostaje H-SORT z warunkiem KB3 — klasyfikacja + decyzje = FINAL |

## §1 — KB4: czy amplituda kreacji rozróżnia Φ od Φ*? (H-CP)

### §1.1 Setup (LIVE machinery, deklaracje)

Fluktuacja zespolona wokół realnego tła ściany: Φ = (w(x) + h) + iχ, konwencja kanoniczna
przy L = ½|∂Φ|² (concept §3.2). Operacja C: Φ → Φ* ⟺ **χ → −χ** (h niezmienione).

### §1.2 Wynik (EXACT, wszystkie rzędy)

- **FP1 (spektrum):** m_h² = λ(3w²−Φ₀²) — odtwarza LOCKED m_eff² (FM P1 FP6) ✓;
  m_χ² = λ(w²−Φ₀²) — mod fazowy: Goldstone w bulku (m_χ²(Φ₀) = 0 ✓), zdestabilizowany
  w całej warstwie (w < Φ₀ ⇒ m_χ² < 0) — spójne z paskiem pre-geometrycznym (SCOPING Q3).
- **FP2 (parzystość):** V(h, −χ) − V(h, χ) = 0 EXACT; rozkład wielomianowy: **zero monomianów
  nieparzystych w χ** — bo \|Φ\|² = (w+h)² + χ² zawiera χ wyłącznie kwadratowo, a człon
  kinetyczny jest parzysty trywialnie. To nie jest własność rzędu — to własność całej akcji.
- **FP3 (audyt wierzchołków):** wszystkie wierzchołki oddziaływań mają PARZYSTĄ liczbę nóg χ
  ({(0,4), (1,2), (2,2), (3,0), (4,0)}) ⇒ każda funkcja korelacji z nieparzystą liczbą
  nóg χ znika ⇒ **żaden proces kreacji z realnej ściany nie może wytworzyć asymetrii Φ/Φ*
  w amplitudzie — w żadnym rzędzie**.

### §1.3 Werdykt i rezyduał

**KB4 = NEGATIVE_FOR_REAL_WALL (DERIVED).** Sakharov-analog: trzeci warunek (C/CP w amplitudzie)
jest w LIVE TGP **strukturalnie niedostępny z realnej ściany** — asymetria, jeśli istnieje,
musi być topologiczno-sortująca (H-SORT), nie amplitudowa. Rezyduał deklarowany (NIE cichy):
tło ściany z teksturą fazową θ(x) ≠ const lub sektory holonomii RP² mogłyby złamać parzystość —
poza LIVE machinery (profil ściany realny per FM D-Q4, deklaracja (vi) FM Phase 0 §8(e));
ewentualny przyszły cykl. Zgodnie z Phase 0 §1.4: gałąź textured/RP² = **GAP deklarowany**.

## §2 — SIG-2: granica leakage i los wyciekłej antymaterii

### §2.1 Struktura kanałów (z P1, moving-wall picture)

Ściana recesuje od materii z prędkością względną c − v_c = c/3 (LOCKED kinematyka). Trzy kanały
losu pary: **(1)** ciasne (ξ_s < ξ_s*): anihilacja w warstwie (energia wraca do ściany);
**(2)** posortowane + antykink wchłonięty przez ścianę podczas tranzytu — **kanał H-SORT**
(+1 materia netto, antymateria w sektorze frontowym); **(3)** posortowane, ale ściana ucieka
antykinkowi — para K–A̅ pozostaje w bulku ⇒ **opóźniona anihilacja w bulku = LEAKAGE**
(źródło rozproszonego tła promienistego).

### §2.2 Domknięcie strukturalne kanału 3 (FP4-FP6)

- **Głębokie związanie (FP4):** \|V_wA(ξ_d)\|/M_K = **12(2−√3) EXACT ≈ 3.21** (Manton;
  pas konwencji amplitudy: 2.48 przy C z fitu P1) — antykink przy kreacji jest związany ze
  ścianą energią PRZEKRACZAJĄCĄ jego masę spoczynkową.
- **Szybka absorpcja (FP5):** τ_acc/τ_tr = **(2+√3)/72 EXACT ≈ 0.052 ≪ 1** (pas kryterium
  Δv ∈ {c, c/3}: 0.017–0.052) — przyspieszenie ku ścianie jest ~20-60× szybsze niż tranzyt warstwy.
- **Próg otwarcia kanału (FP6):** absorpcja zawodzi dopiero dla kreacji na głębokości
  **ξ > ξ_leak = ln 72 ≈ 4.28** (pas: do ln 216 ≈ 5.38). Ale **A-i (LOCKED) ogranicza kreację
  do warstwy konwersji** (waga ∝ V(w(x)) ∝ sech⁴(ξ/2), skupiona przy ξ ~ 1.3) ⇒ kanał leakage
  jest karmiony wyłącznie wykładniczym ogonem wagi: **f_leak = [2/3 − tanh(ln72/2) +
  tanh³(ln72/2)/3]/(2/3)**, z tanh(ln72/2) = 71/73 EXACT ⇒ **f_leak ≈ 1.1×10⁻³**
  (konserwatywnie; pas 1.3×10⁻⁴–1.2×10⁻³).

### §2.3 Los wyciekłej antymaterii + porównania (FP7)

Wyciekła para K–A̅ w bulku jest związana przyciąganiem (LOCKED V_int) i w 1D nie ma kanału
ominięcia ⇒ anihiluje w skończonym czasie ⇒ **frakcja TRWAŁEJ antymaterii w bulku → 0**
(transient). Porównania (mechaniczne, comparison-only — anchor wyłącznie w tej linii):

- **Literal anchor (pre-rejestrowany):** trwała frakcja antymaterii 0 < f̄_max = 10⁻⁶ — **PASS**.
  H-SORT przewiduje brak trwałych domen antymaterii w bulku „za darmo" (spójność klasy F6/F9).
- **Nowa obserwabla (flagowana, NIE porównywana — anti-goalpost):** wtrysk radiacyjny
  z opóźnionych anihilacji: **f_rad ≈ 2f_leak ≈ 2.2×10⁻³** energii kreowanej → rozproszone
  tło promieniste w bulku. Brak pre-rejestrowanego anchora dla tej obserwabli ⇒ porównanie
  ZAKAZANE w tym cyklu; **to jest naturalna treść falsyfikatora PR-023** (przyszła
  pre-rejestracja z własnym anchorem). Uczciwie: rząd 10⁻³ może być za duży lub w sam raz —
  rozstrzygnie przyszły cykl z uczciwym anchorem.

### §2.4 SIG-1 refinement (FP8) + deklaracje

Mnożnik popytu przy kreacji parowej ≥ 2 (kanały anihilacyjne 1 i 3 nie dają materii netto;
kanał 1 recyklinguje energię do ściany) ⇒ **t_*^(B) ≥ √2·t_*** — wynik P1 staje się dolną
granicą (monotoniczność √mult). Deklaracje per-use: model wagi kreacji (∝ sech⁴ — MODEL,
nie wyprowadzenie z kinematyki kreacji); kryterium pościgu zgrubne (pas konwencji raportowany);
poprawka od przyciągania partnera subdominująca (e^{−ξ_s} < e^{−ξ_d} dla ξ_s > ξ_s*);
**1D-proxy ≠ 3D claim** (w 3D możliwe ominięcia — pas f_leak to stwierdzenie proxy).

## §3 — Stan modułu B po Phase 4 (przed FINAL)

| Kryterium | Status |
|---|---|
| KB1 | PASS (selekcja topologiczna + energetyka) |
| KB2 | PASS-DERIVED (kierunek sortowania) |
| KB3 = SIG-3 | CONDITIONAL (ξ_s > ln(3+2√3) EXACT; nie pełny zakres) |
| KB4 | **NEGATIVE_FOR_REAL_WALL** (H-CP wykluczone w LIVE) + GAP (textured/RP²) |
| SIG-1 | PASS EXACT (≥ √2·t_* po refinement) |
| SIG-2 | BOUNDED (f_leak ≲ 1.2×10⁻³; literal anchor PASS; f_rad → PR-023 kandydat) |

**Zawężenie zbioru hipotez (mechaniczne):** H-CP — wykluczone dla realnej ściany;
HONEST_NEGATIVE — wykluczone (wymaga KB1 null AND KB4 null; KB1 jest pozytywne);
pozostaje **H-SORT z dokładnym warunkiem istnienia (KB3)**. Klasyfikacja F-BA-6 w zbiorze
CLOSED + decyzja progowa KB3 = **Phase FINAL (wyłącznie user)**. **NO PR-023** (forbidden
move #8). Moduł A: bez zmian (strict BRIDGE_PARTIAL — decyzja GAP-3 w FINAL).

## §4 — Anti-Lakatos (Phase 4): COMPLIANT ✓

0/9 hardcoded; anchor f̄_max wyłącznie w linii porównawczej FP7 (audit FP9: nieobecny we
wszystkich wyrażeniach derywacyjnych); nowa obserwabla radiacyjna flagowana BEZ porównania
(zakaz wymyślania anchorów mid-cycle); KB4 NEGATIVE raportowane wprost (zamyka H-CP — wbrew
pokusie utrzymania dwóch żywych hipotez); pasy konwencji (amplituda C, kryterium pościgu)
raportowane zamiast wyboru korzystnego punktu; naprawy plumbingu (konwencja normalizacji
kanonicznej §3.2 + funkcja pierwotna sech⁴ z weryfikacją) udokumentowane bez zmiany progów;
LOCKED read-only; 0 predecessor modified; 0 nowych stałych.
