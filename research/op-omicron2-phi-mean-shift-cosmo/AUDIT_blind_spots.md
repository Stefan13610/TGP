---
title: "AUDYT: gdzie tkwi blind spot? — 13 kandydatów + 1 STRUCTURAL CONTRADICTION"
date: 2026-05-03
parent: "[[README.md]]"
type: audit
status: in_progress
trigger: "User intuition: Phi_0 should track global matter — gap of 8 orders cannot be right"
---

# AUDYT: blind spots w analizie Φ_0 / Λ_eff cosmologicznej ewolucji

> **Trigger:** User pisze 2026-05-03: "intuicja podpowiada mi, że czegoś tutaj
> nie dostrzegam, ale nie wiem czego". Quick check (Opcja A) potwierdził
> M10.5 verdict (gap 8 dekad), ale gap jest TAK ogromny vs sek01 ontology
> claim, że coś musi być pominięte.

## 🔴 CENTRAL FINDING: Strukturalna sprzeczność sek01 ↔ sek08a

To jest **najgłębsze odkrycie tego audytu** — nie jest "blind spot" calc, to jest
**inkonsystencja TGP-foundationowa** między dwiema warstwami:

### sek01_ontologia.tex (linia 681-683):

> "Φ_0 **JEST generowane przez materię** — jest po prostu **średnim wynikiem
> tej generacji**, przyjętym za punkt odniesienia."

### sek01_ontologia.tex (Tablica stanów, linia 700-703):

> "Tło referencyjne Φ_0: **Średni stan generowany przez kosmologiczną zawartość
> materii**. 'Próżnia' TGP: brak lokalnych źródeł, ale przestrzeń istnieje
> (**generowana globalnie**)."

### sek08a_akcja_zunifikowana.tex (formalism):

```
V(Φ) = (β/3)Φ³ - (γ/4)Φ⁴
V'(Φ_0) = 0  →  Φ_0 = β/γ = 1·Φ_0  (gdy β=γ vacuum cond.)
```

**Φ_0 fixed jako vacuum minimum, NIEZALEŻNIE od matter content.**

### Sprzeczność:

| Warstwa | Co mówi o Φ_0 | Implikacja |
|---|---|---|
| sek01 ontology | Φ_0 = **<Φ>_globalnej kosmicznej materii** | Φ_0(t) **zależne** od ρ̄(t) |
| sek08a formalism | Φ_0 = vacuum minimum V'(Φ_0)=0 | Φ_0 **niezależne** od ρ |
| Closure 2026-04-26 T-Λ | Φ_eq = ℏH_0 (Hubble cutoff), today's H_0 | Φ_0 ~ H_0² (TODAY) |

Te trzy są w sprzeczności. T-Λ wybiera "today's H_0" arbitralnie — w obu pozostałych
ramach są inne natural choices.

---

## 13 KANDYDATÓW na blind spot — pełna lista

### A. Direct ρ-Φ coupling beyond ax:metric-coupling

**Hipoteza:** Może istnieje DODATKOWY direct coupling ρ↔Φ poza minimal `q·φ·ρ`.
- ax:metric-coupling dotyczy fields ψ_m (level 3a)
- Density ρ (level 3b) ma osobne sprzężenie minimalne
- Mógłby być czlon ~ρ²/Φ_0², ρ·∇Φ, etc.

**Magnitude potential:** **Średnia** — ograniczone aksjomatami, ale możliwe.
**Modyfikacja TGP:** TAK (poszerza coupling)
**Test:** sympy zmodyfikowane S, oblicz Φ-EOM, porównaj z PPN.

### B. Sek05 §de-struktury vs Buchert — same magnitude, but different mechanism

**Hipoteza:** Sek05 mówi o LOCAL "Φ_strukt > Φ_homog", Buchert obliczona Q_D.
Czy to ten sam mechanism kwantytatywnie?

**Quick check (Opcja A) result:** TAK, oba ~10⁻⁸. Buchert formula captures sek05.
**Magnitude potential:** **Mała** — confirmed.
**Status:** **CLOSED** — to NIE jest blind spot.

### C. c(Φ) cosmologicznie zmienne (user's hypothesis)

**Hipoteza:** sek04 mówi że c, ℏ, G stają się Φ-zależne. Może c(z) ≠ c(today).
**Quantitative:** δc/c ~ δΦ/(2Φ) — bounded by Φ perturbations ~ 10⁻⁵.
**Magnitude:** **Bardzo mała** (10⁻⁵).
**Status:** **CLOSED** — niewystarczające.

### D. β/γ matter-density-dependent

**Hipoteza:** Parametry V mogą być effective, zależne od `<ρ>`. Wtedy V minimum
przesuwa się z `<ρ>(t)`.

**Magnitude potential:** **Duża** — może dać dowolny shift Φ_0(t).
**Modyfikacja TGP:** **TAK głęboka** (β, γ tracą status fundamental constants).
**Test:** Sympy: V(Φ, ρ) = V_0(Φ) + α·ρ·δV(Φ), oblicz Φ_0(ρ), porównaj.
**Risk:** może kolidować z particle physics predictions (mass spectrum z β, γ).

### E. H(z) modification w TGP

**Hipoteza:** TGP modifies expansion history przy z konkretnym, zmienia H(z).
**M10.1 already check:** w(z) ≈ -1 numerically (10⁻⁴). H(z) ~ ΛCDM.
**Status:** **CLOSED** — tym kanałem nie ma drogi.

### F. UV.3 cosmological flow

**Hipoteza:** Z_Φ = 14/3 ma cosmological dynamics, nie tylko UV→IR static.
**Previous discussion:** Z_Φ jest static algebraic, η = 0.044 daje 4% per dekada,
ale nie temporal flow.
**Magnitude:** **Mała** (<5%).
**Status:** **WEAK** — niewystarczające bez dodatkowej hipotezy.

### G. Local-frame vs cosmological-frame measurement

**Hipoteza:** SH0ES mierzy w lokalnej overdensity, Planck w globalnym CMB frame.
**Local Group / Virgo Supercluster Ψ_local ~ 10⁻⁵ - 10⁻⁷.**
**Magnitude:** **Mała** (10⁻⁵).
**Status:** **CLOSED** — niewystarczające.

### H. VSL (Variable Speed of Light)

**Hipoteza:** c(t) varies cosmologically (Albrecht-Magueijo, Moffat).
**Magnitude w TGP:** δc/c ~ δΦ/Φ ~ 10⁻⁵.
**Status:** **CLOSED** — niewystarczające.

### I. Foreground/background separation

**Hipoteza:** "Λ_observed" = funkcja Φ_0 (vacuum) AND <δΦ²> (perturbations).
Może ich kompozycja jest źle rachowana.
**Quick check:** powtórzy sek05 lemma + M10.5. Same answer.
**Status:** **CLOSED**.

### J. Nonlinear amplification beyond perturbation theory

**Hipoteza:** TGP field eqn ma nonlinearity `(∇Φ)²/Φ` która może amplify NL.
**Standard PT:** lowest correction is δ³, daje variance³ ~ 10⁻²⁴ (smaller).
**Non-perturbative:** soliton configurations z innym vacuum branch.
**Magnitude:** **Niewiadoma** — non-PT calc needed.
**Status:** **OPEN, low priority** — speculative.

### K. 🔴 **Φ_0(t) ∝ <ρ>(t) — sek01 ontology poza-formalism implication**

**Hipoteza:** Sek01 ontology mówi Φ_0 = `<Φ>_globalnej materii`. Jeśli formalism
sek08a tego NIE captures, to brakuje **fundamental piece** TGP.

**Quantitative test (NOWY rachunek):**

Z sek08a Φ-EOM w jednorodnym tle ρ̄(t):
```
∇²Φ + ... + βΦ²/Φ_0 - γΦ³/Φ_0² = -q Φ_0 ρ
```

W FRW homogeneous: ∇² → 0 (translacyjnie niezmienne tło). Φ jest tylko czasem zależne.
Statyczna granica: 0 + βΦ²/Φ_0 - γΦ³/Φ_0² = -q·Φ_0·ρ̄

Z β=γ i przy Φ ≈ Φ_0(1+δ_bg):
```
βΦ_0(1+2δ_bg) - γΦ_0(1+3δ_bg) = -q·Φ_0·ρ̄
0 + (2β-3γ)Φ_0·δ_bg = -q·Φ_0·ρ̄
-γ Φ_0 δ_bg = -q Φ_0 ρ̄
δ_bg = q ρ̄ / γ
```

Z `q·Φ_0 = 4πG/c²` and `γ ~ H_0²/c_0²`:
```
δ_bg = (4πG/c²)·ρ̄/(H_0²/c_0²) = 4πG·ρ̄/H_0² = (3/2)·Ω_m(z)
```

🔴 **TO JEST OGROMNY EFEKT!**
- Dziś: δ_bg ≈ (3/2)·0.315 = 0.47 (~50% odchylenie!)
- Recomb: Ω_m(z=1100) ≈ 1, δ_bg ≈ 1.5
- **ΔΦ_0/Φ_0 między recomb i dziś ≈ 1.5 - 0.47 = 1.03 (factor of 2 change!)**

**WAIT — to jest 50% NA POZIOMIE TŁA, w current sek08a equations.**

Ale to jest skoro V z β=γ jest unstable maximum! Linearyzacja wokół ψ=1
(maximum) nie daje stabilnego rozwiązania, tylko slow-roll. M10.5 pokazało że
TEMPORAL evolution daje slow-roll z `ψ → 1` attractor.

W tym ujęciu rachunek powyżej (δ_bg = (3/2)Ω_m) jest **statycznym fixed-point**
przy obecności ρ̄, nie cosmological evolution.

**Open question:** Czy TGP cosmological FRW solution rzeczywiście tracks ρ̄ ten sposób,
czy też ψ → 1 attractor dominuje i Φ_0 staje się niezależne od ρ?

**To jest KLUCZOWE pytanie — wymaga dedicated rachunek.**

**Magnitude potential:** **OGROMNA** (do 100% Φ_0 evolution)
**Modyfikacja TGP:** Może NIE — może to tylko inny INTERPRETATION istniejących równań
**Test:** Eksplicit FRW Φ-EOM solver z ρ̄(t) source, śledzenie Φ_0(t)

### L. γ ~ H(z)² instead of H_0² (full tracker quintessence)

**Hipoteza:** T-Λ closure picks γ ~ M_Pl² · H_0² (today's). Why not γ ~ H(z)²?

**Magnitude:** Λ_eff ~ M_Pl²·H(z)² ~ ρ_total(z) — full tracker.
- z=1100: Λ_eff ~ 10⁹·Λ_today (way too much)
- Jeśli partial: Λ_eff = α·M_Pl²·H_0² + (1-α)·M_Pl²·H(z)², dla α ~ 0.99 daje 1% effect
- Dla 25% effect: α ~ 0.75 — moderately tuned

**Modyfikacja TGP:** TAK — fundamental change to T-Λ closure.
**Test:** explicit construction; check if any TGP principle picks H_0 vs H(z).

### M. Hidden NEW coupling — completely outside current axioms

**Hipoteza:** Maybe ax:metric-coupling itself is too restrictive. Allow material
fields to couple to ∇Φ or Φ̇ directly.

**Magnitude:** Unknown.
**Modyfikacja TGP:** **TAK głęboka** — modifies TGP_FOUNDATIONS §4.
**Test:** speculative.

---

## RANKING kandydatów (priority for investigation)

| # | Kandydat | Magnitude | Plausibility | Test cost | Priority |
|---|---|:---:|:---:|:---:|:---:|
| **K** | **Φ_0(t) ∝ <ρ>(t) — sek01/sek08a contradiction** | 🔴 OGROMNA | 🟢 sek01 explicit | 1-2h | **🔴 HIGHEST** |
| L | γ ~ H(z)² instead of H_0² | 🟡 Średnia | 🟡 możliwa | 1h | **🟡 MEDIUM** |
| D | β/γ matter-density-dependent | 🟡 Duża | 🔴 modifies fund | 4-8h | 🟡 MEDIUM |
| A | Direct ρ-Φ coupling (extra) | 🟡 Średnia | 🟡 plausible | 2-4h | 🟡 MEDIUM |
| F | UV.3 cosmological RG flow | 🔴 Mała | 🟡 needs RG eqs | 4-8h | 🟢 LOW |
| J | Non-PT amplification | ❓ Unknown | 🔴 speculative | weeks | 🟢 LOW |
| M | Outside-axiom new coupling | ❓ Unknown | 🔴 violates fund | unknown | 🟢 LOW |
| B,C,E,G,H,I | Already tested ~10⁻⁵ | 🔴 Mała | — | — | **CLOSED** |

---

## KLUCZOWE PYTANIE — punkt zaczepienia dla użytkownika

**Twoja intuicja może wskazywać dokładnie KANDYDATA K:**

Sek01 ontology mówi explicit: **"Φ_0 jest generowane przez kosmologiczną zawartość
materii"**. Sek08a formalism traktuje Φ_0 jako fixed parameter. Te dwa ŚCIERAJĄ
SIĘ **w niemniej fundamental way**.

Pytanie: **Czy sek08a Φ-EOM solved cosmologicznie (homogeneous ρ̄(t) source)
faktycznie daje Φ_0(t) ≠ const?**

Quick algebraic estimate **w tym audicie** (kandydat K) sugeruje że:
```
δ_bg = (3/2)·Ω_m(z)
```
co byłoby **OGROMNYM** efektem (50%-100% Φ_0 evolution).

ALE — to jest fixed-point statyczny, NIE cosmological history. Trzeba explicit
solver FRW Φ-EOM z time-dependent ρ̄(t), z proper initial conditions, by
zobaczyć czy:
- (a) ψ → 1 attractor dominuje → Φ_0 = const (current TGP picture)
- (b) Φ_0 tracks ρ̄(t) → ogromny Hubble tension shift
- (c) Mid-ground: częściowe tracking, daje moderate effect

**To jest kandydat K test. Sugeruje go zbadać NEXT.**

---

## Następny krok — rekomendacja

**Opcja K-test (~1-2h):** explicit FRW Φ-EOM solver z homogeneous ρ̄(t):

```python
"""
omicron2_K_test.py — czy sek01 ontology realnie wpływa na cosmological Φ_0(t)?

Solve FRW Φ-EOM:
  Φ̈ + 3H Φ̇ + V'(Φ)/K(Φ) = -q ρ̄(t)·Φ_0/K(Φ)
  H² = (8πG/3)·(ρ̄ + ρ_Φ)
  ρ_Φ = (1/2)K(Φ)Φ̇² + V(Φ)
  ρ̄(t) = ρ_m(z) + ρ_r(z)
  
Question: does Φ(t) track ρ̄(t), giving cosmological Φ_0 evolution?
Or does ψ → 1 attractor lock Φ at vacuum, making Φ_0 = const?

If (a) tracking → Hubble tension solution candidate.
If (b) attractor → M10.5 verdict reinforced.
"""
```

Albo:

**Opcja przerwania:** zostawić to jako research note z 13 candidates documented,
flag K jako TOP priority, czekać na user decision.

Co wolisz?
