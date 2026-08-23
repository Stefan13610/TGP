---
title: "Phase0_balance — LOCK planu op-spectral-analysis-Phi (CP-7 / L03)"
date: 2026-07-03
type: phase0-lock
tgp_owner: research/op-spectral-analysis-Phi-2026-07-03
status: LOCKED (przed jakimkolwiek obliczeniem Phase 1/2)
anti_lakatos_lock: PRESERVED
tags: [L03, CP-7, spectral, BVP, ghost-wall, lock]
related:
  - "[[../../audyt/L03_K_phi_stability/README.md]]"
  - "[[../op-L03-spectral-stability-2026-05-06/spectral_synthesis.md]]"
  - "[[../../meta/AUDYT_GLEBOKI_2026-06-28.md]]"
---

# Phase 0 — balance & LOCK (przed obliczeniami)

## Cel cyklu (z AUDYT_GLEBOKI CP-7, werdykt L03 = 🟡 PARTIAL)

> „Diagonalizacja K(φ)∂²+V''(φ) na tłach próżnia/Yukawa/soliton (sympy +
> BVP numeryczny), m²>0 dla φ>0, rozstrzygnięcie ghost-wall przy φ→0."

Poprzedni cykl (`op-L03-spectral-stability-2026-05-06`) dostarczył
**analityczną** syntezę S-L (pozytywność punktowa p, q + tw. Coddington-
Levinson). Audyt głęboki: *„tylko punkt próżni domknięty; brak analizy
spektralnej na tłach niejednorodnych, ghost-wall przy φ→0 nietknięte"*.
Ten cykl = faktyczna diagonalizacja numeryczna.

## Stan faktyczny formulacji (inwentaryzacja PRZED obliczeniem — bez relitygacji L04)

Współistnieją w repo DWIE samospójne pary (funkcjonał ↔ tło), które
NIE są równoważne polowo (znany dualizm L04; sek08b `rem:formulation-dictionary`):

**F-A (grawitacyjna kanoniczna, α=2; sek08_formalizm `thm:field-eq`):**
```
E_A[ψ] = ∫ [ ½·K_geo·ψ⁴·|∇ψ|² + U_A(ψ) ] d³x,
U_A'(ψ) = K_geo·γ·(ψ⁷ − ψ⁶)          (β=γ, próżnia ψ=1)
```
EL(E_A) ⇔ eq:full-field: ∇²ψ + (2/ψ)(∇ψ)² + γψ² − γψ³ = −qΦ₀ρ.
Linearyzacja rdzeniowa: m_sp² = 3γ−2β = γ (prop:vacuum-stability).

**F-S (solitonowa log-forma α=2 — TA, KTÓREJ UŻYWA KORONA:
`a3d_soliton_brannen_r.py`, `ls10_third_generation_selection.py`):**
```
E_S[g] = ∫ [ ½·f(g)·|∇g|² + W(g) ] d³x,
f(g) = 1 + 4·ln g   (ghost wall g* = e^{−1/4} ≈ 0.7788, f(g*)=0),
W'(g) = g²(1−g)     (W''(1) = −1)
```
EL(E_S) ⇔ ODE a3d: f·(g''+2g'/r) + (2/g)g'² = g²(1−g).
UWAGA zarejestrowana PRZED obliczeniem: (i) sek08b preferuje substratową
K=g² (α=1, `cor:alpha1-preferred`), ale skrypty korony mają ALPHA=2.0
(log-forma) i mechanizm selekcji N=3 UŻYWA odbić od ściany g*;
(ii) W''(1)=−1<0 ⇒ w F-S kontinuum wokół g=1 startuje od wartości
ujemnej (ogon oscylacyjny) — znany element dualizmu L04, będzie
ZMIERZONY i zaraportowany, nie „naprawiany".

**F-S′ (substratowa α=1, preferowana przez sek08b):**
```
f_sub(g) = K_geo·g², W'(g) = g²(1−g)
```
Diagnostyka równoległa do F-S (porównanie dyktowane słownikiem 12/12).

## Operator fluktuacji (do wyprowadzenia sympy w Phase 1 — forma oczekiwana)

Dla E = ∫ ½F(u)|∇u|² + W(u), tło u₀(r), sektor kątowy l:
```
L̂[v] = −(1/r²)·d/dr[ r²·F(u₀(r))·dv/dr ] + [ Q(r) + F(u₀)·l(l+1)/r² ]·v = λ·v
Q(r) = W''(u₀) − ½F''(u₀)·(u₀')² − F'(u₀)·[u₀'' + (2/r)u₀']
```
Phase 1 MUSI wyprowadzić Q niezależnie (druga wariacja + całkowanie
przez części) i potwierdzić formę; rozbieżność = STOP i korekta LOCK-a
z adnotacją (dozwolona TYLKO korekta algebry operatora, nie kryteriów).

## Tła (LOCKED)

| ID | Formulacja | Tło | Parametry |
|----|-----------|-----|-----------|
| B1 | F-A | próżnia ψ≡1 | — |
| B2 | F-A | Yukawa ψ=1+A·e^(−√γ·r)/r | A ∈ {0.01, 0.1, 0.3, 1.0} |
| B3e | F-S | soliton g(r), g(0)=g₀^e=1.24915, g'(0)=0 | ODE a3d, r_max=60 |
| B3μ | F-S | soliton g(0)=φ·g₀^e≈2.02100 | jw. |
| B3τ | F-S | soliton g(0)=3.18912 (Koide) | jw. (z odbiciami od ściany jak a3d) |
| B4 | F-S | próżnia g≡1 (kontrola dualizmu) | — |
| B5 | F-S′ | soliton α=1 substrat, te same g₀ | ODE substratowe |

Jednostki: K_geo = γ = β = 1. Yukawa: profil liniowej odpowiedzi
(δψ>0, rem:sign-convention); dla A=1.0 sprawdzić |δψ| poza reżimem
liniowym — wtedy relaksacja nieliniowa pełnego EOM (jeśli osiągalna,
inaczej odnotować ograniczenie).

## Kryteria PASS/FAIL (LOCKED — żadnych zmian po uruchomieniu)

- **C1 (sympy, F-A):** druga wariacja E_A na ψ=1 daje dokładnie
  m_sp² = γ (sympy simplify → 0). PASS = tożsamość exact.
- **C2 (BVP, B1):** całe σ(L̂) ≥ −tol_num; krawędź kontinuum = 1 ± 2%
  (jedn. γ=1). tol_num = 10⁻⁸·max(1, |λ_max|).
- **C3 (BVP, B2):** dla KAŻDEGO A ze skanu: N_neg = 0
  (żadnego λ < −tol_phys, tol_phys = 10⁻⁶). PASS = 4/4.
- **C4 (BVP, B3e/μ/τ, l=0):** REPORT N_neg (liczba modów λ < −10⁻⁶
  potwierdzona zbieżnością siatki). Kryterium opisowe (bez wishful):
  N_neg=0 ⇒ stabilność solitonu potwierdzona spektralnie (wspiera koronę);
  N_neg≥1 ⇒ **wynik NEGATYWNY zgłaszany wprost** (konfrontacja z
  `prop:nonlinear-stability`; wpis do Limitations korony; NIE unieważnia
  automatycznie dopasowań mas — unieważnia claim „pełna stabilność
  spektralna").
- **C5 (BVP, B4):** oczekiwane analitycznie kontinuum od W''(1)/f(1) = −1;
  numeryka MA to potwierdzić (to dokumentacja dualizmu L04 u źródła,
  wynik „ujemny" tu NIE jest zaskoczeniem ani FAIL-em cyklu — jest
  zmierzonym faktem o F-S).
- **C6 (ghost wall):** (a) min_r g(r) dla B3e/μ/τ vs g*=0.7788; czy
  f(g(r)) zmienia znak NA tle (utrata eliptyczności operatora —
  konsekwencje spektralne udokumentować); (b) analitycznie: klasyfikacja
  osobliwego końca ψ→0 dla F-A (K=ψ⁴; zmienna kanoniczna χ=ψ³/3,
  potencjał efektywny, limit-point vs limit-circle). PASS = rozstrzygnięte
  (dowolny kierunek).

## Forbidden moves (anti-Lakatos)

1. Zmiana kryteriów / tolerancji / listy teł PO uruchomieniu obliczeń.
2. Tuning g₀, r_max, siatki pod wynik (wartości wyżej = zalockowane;
   testy wrażliwości TYLKO jako dodatkowe: N ∈ {2000,4000,8000},
   r_max ∈ {40,60,80} — wszystkie raportowane).
3. Odrzucenie modu ujemnego jako „artefakt numeryczny" bez testu
   zbieżności na 3 siatkach.
4. Relitygacja L04/α=2 (zamknięte #49/#53) i mechanizmu korony (#56).
5. Edycje core/.tex w tej sesji (propozycje → NEEDS.md, user-gated).

## Warunki brzegowe (LOCKED)

l=0: regularność v'(0)=0 (staggered grid r_j=(j+½)h); Dirichlet v(R)=0.
Dyskretyzacja: FD 2. rzędu formy samosprzężonej
−(1/r²)[r²F v']' (zachowuje symetrię macierzy po symetryzacji wagą);
solver: scipy eigh_tridiagonal / eigsh.

## Deklaracja

Ten LOCK zapisano PRZED napisaniem i uruchomieniem jakiegokolwiek kodu
Phase 1/2. Wszystkie wyniki (również negatywne dla korony/F-S) zostaną
zaraportowane w README + STATE.md.
