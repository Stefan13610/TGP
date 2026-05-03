---
title: "G.0 Phase 1 — Structural setup: 4 EOM diagnostics + 3 reduction strategies"
date: 2026-05-02
phase: 1
parent: "[[README.md]]"
status: ACTIVE
tags:
  - TGP
  - G0
  - phase1
  - structural-setup
  - sek08-audit
  - R3
  - M9.1pp
---

# G.0 Phase 1 — Structural setup

## Cel Phase 1

Rozstrzygnąć **strukturalnie** (na poziomie sympy/numerycznym, ale bez
pełnych dowodów) która z trzech hipotez (H1/H2/H3 z [[README.md]]) najlepiej
opisuje relację między R3 ODE a TGP-canonical na metryce M9.1''.

**Score gate:** ≥1/3 PASS w sub-tasks G0a/G0b/G0c → Phase 2 forward
(formalizacja zwycięskiej hipotezy w sympy LOCK + testy aniologiczne).
0/3 PASS → Phase 2 redirect do testowania H3 (TGP-canonical jest błędne).

---

## 0. Punkt wyjścia — diagnoza 4 EOM (rerun na dysku 2026-05-02)

Krytyczne dla Phase 1: musimy wiedzieć jakie EOM dokładnie porównujemy.
Cztery konkurencyjne równania na statyczny soliton sferycznie symetryczny:

### 0.1. EOM (a) — TGP_FOUNDATIONS *claim*

**Źródło:** `TGP_FOUNDATIONS.md` lin. 80–86.

```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀² = -qΦ₀ρ
```

W bezwymiarowych ψ = Φ/Φ₀, β=γ, ρ=0:

```
ψ'' + (2/r)ψ' + 2(ψ')²/ψ + γψ²(1-ψ) = 0
```

**Status:** *claimed* canonical Φ-EOM TGP. Brak explicit derivation z
S_TGP w plikach core (potwierdzone bezpośrednim odczytem 2026-05-02).

### 0.2. EOM (b) — sek08a wariacja z M9.1 (FALSIFIED metric)

**Źródło:** `sek08a` lin. 244–298 (`prop:field-eq-from-action`),
z `√(-g_eff) = c₀ψ` wyprowadzonym z M9.1 (FALSIFIED 2026-04-25).

Po pełnym rachunku statyczny EOM:

```
ψ'' + (2/r)ψ' + 2(ψ')²/ψ = γ(15ψ-16)/(12ψ)        [β=γ]
```

**Status:** **derived correctly z S_TGP, ale na FALSIFIED metryce M9.1.**
Co istotne: **różni się od (a)** — sek08a *claim* o reprodukcji canonical
EOM jest błędny.

### 0.3. EOM (c) — wariacja S_TGP z poprawnym √(-g) M9.1''

**Źródło:** ten dokument, rerun 2026-05-02. Z M9.1'':
- `g_tt = -c₀²(4-3ψ)/ψ`, `g_rr = ψ/(4-3ψ)`
- `g_eff^rr = (4-3ψ)/ψ`
- `√(-g) = c₀ψ/(4-3ψ)`

Z `K(ψ)=ψ⁴` i `V(ψ)=γψ³(4-3ψ)/12` (z β=γ identity):

```
ψV(ψ)/(4-3ψ) = γψ⁴/12        ← (4-3ψ) idealne skraca!
S_static = 4πc₀ ∫r²dr × [½ψ⁴(ψ')² − γψ⁴/12]
EL → ψ'' + (2/r)ψ' + 2(ψ')²/ψ = −γ/(3ψ)
```

**Status:** **co naprawdę wynika z S_TGP po update sek08 do M9.1''.**
Ma jeden ekstrem w ψ=0 (singular), brak vacuum w ψ=1 — **niefizyczne!**
To wskazuje, że albo K(ψ)=ψ⁴ jest błędne dla M9.1'', albo V(ψ) musi być
inne, albo całe S_TGP wymaga reformulacji.

### 0.4. EOM (d) — R3 ODE (mass spectrum 0.001%)

**Źródło:** `research/why_n3/r3_phase1_psi_g0_identification.py` (α=2):

```
g'' + (2/r)g' + 2(g')²/g + (g-1)/g² = 0
```

**Wariacyjne pochodzenie:** `L = ½g⁴(g')² - V(g)` z `V = g³/3 - g⁴/4`,
**na płaskim 3D background, BEZ √(-g) = c₀ψ wegetującego sprzężenia z metryką**.

```
EL z L = ½g⁴(g')² − V(g):
  g⁴g'' + 2g³(g')² + 2g⁴/r·g' = V'(g)
  g'' + 2(g')²/g + (2/r)g' = V'/g⁴
  
Dla V = g³/3 - g⁴/4:
  V' = g² - g³ = g²(1-g)
  V'/g⁴ = (1-g)/g²        ← dokładnie R3 RHS!
```

**Status:** R3 ODE **JEST** wariacyjnie spójne, ale na **flat 3D background
(η_ij)** — nie na M9.1''. To jest centralny insight Phase 1.

### 0.5. Tabela porównawcza

| EOM | RHS (β=γ=1) | Vacuum punkt | Bariera | Wariacyjnie spójne? | Empirycznie? |
|---|---|---|---|---|---|
| (a) claim | `-ψ²(1-ψ)` | ψ=1 ✓ | ψ→0 | ❓ (no proof) | — |
| (b) M9.1 | `γ(15ψ-16)/(12ψ)` | ψ=16/15 ✗ | ψ→0 | ✓ | sprzeczne z (a) |
| (c) M9.1'' | `−γ/(3ψ)` | brak (singular) | ψ→0 | ✓ | NIEFIZYCZNE |
| (d) R3 | `−(1-g)/g²` | g=1 ✓ | g→0, g₀_crit=1.874 | ✓ (flat 3D) | mass 0.001% PDG |

### 0.6. Strategiczny wniosek

R3 ODE (d) jest **jedyną** wariacyjnie spójną EOM dającą prawidłowe vacuum
i empirię. Jednocześnie wynika z Lagrangianu **bez** √(-g) M9.1''. To
sugeruje 3 możliwe interpretacje:

1. **(H1)** R3 jest *frame-redefined* version EOM (c) — jest pewna
   transformacja `(ψ,r) → (g,r̃)` która eliminuje (4-3ψ)-related anomalie
   w (c) i daje (d).
2. **(H2)** R3 jest *Einstein-frame projection* — po Weyl rescaling
   `g_eff → Ω²·η`, EOM dla soliton excitation reduce się do (d).
3. **(H3)** R3 jest *niezależnym formalizmem* — żadna transformacja nie
   reduce EOM (a), (b) lub (c) do (d). TGP wymaga reformulacji warstwy 1
   (akcja) lub akceptacji DWÓCH formalizmów (gravity ≠ matter).

Phase 1 testuje wszystkie 3 hipotezy systematicznie.

---

## 1. Sub-task G0a — Volume integration test (H1, primary)

### Cel

Sprawdzić czy wariacja S_TGP z poprawnym M9.1'' volume element redukuje się
do R3 ODE w pewnym sensowanym limicie (np. dominant balance dla solitonu,
lub po wyłączeniu jednego z członów potencjału).

### Strategia

**Etap 1:** Symbolicznie wyprowadzić EOM (c) z `S_TGP[Φ, g_eff^M9.1'']` z
różnymi formami `K(ψ)` i `V(ψ)`:

| Wariant | K(ψ) | V(ψ) | EOM |
|---|---|---|---|
| K1 | ψ⁴ (oryginalne) | γψ³(4-3ψ)/12 (z β=γ) | (c): `−γ/(3ψ)` |
| K2 | ψ⁴ × (4-3ψ)/ψ | γψ³(4-3ψ)/12 | TBD (rerun) |
| K3 | [ψ/(4-3ψ)]² | γψ³(4-3ψ)/12 | TBD |
| K4 | ψ⁴ | inne V (np. log) | TBD |

**Etap 2:** Dla każdego wariantu sprawdzić:
- Czy istnieje stabilne vacuum (V'(ψ_v)=0, V''(ψ_v)>0)?
- Czy istnieje topologiczna bariera (singularity solitonu)?
- Czy RHS EOM ma postać `±(1-ψ)/ψ²` (R3-like) lub redukuje się do tego
  w limicie ψ→1?

**Etap 3:** Linearize wokół ψ=1 i porównaj soliton tail decay z R3:
- R3: `g(r) - 1 ~ A_tail/r` z `m_eff² = 1` (tail decay)
- TGP-canonical (a) claim: `m_eff² = γ` (z lin. of `−γψ²(1-ψ)`)
- (c) M9.1'' z K1: `m_eff²` z lin. `−γ/(3ψ)` wokół ψ=1: `δψ/3` → `m² = 1/3`?

### PASS criterion

**Znalezienie kombinacji `(K, V, √(-g))` dającej:**
1. Stabilne vacuum w ψ=1 ✓
2. Topologiczna bariera w pewnym ψ_crit ✓
3. RHS EOM = `(1-ψ)/ψ²` (lub ekwiwalentne po reparametryzacji liniowej)
4. Mass spectrum α=2 zachowany (m_μ/m_e = 206.77)

≥ 3/4 z tych = G0a PASS.

### Plik

`phase1_G0a_volume_integration.py` — sympy + numerical, przewidywany
~250–350 linii.

---

## 2. Sub-task G0b — Field redefinition test (H1 alt + H2)

### Cel

Symbolicznie szukać transformacji `ψ = T(g)` (z `r' = ρ(r,g)` opcjonalnie)
która przekształca jedną z EOM (a)/(b)/(c) w R3 ODE (d).

### Strategia

**Etap 1:** Substytucja **liniowa** `ψ = a·g + b` (znana z PHASE1 why_n3:
a=0.3814, b=0.6186 mapuje g₀_crit=1.874 ↔ ψ=4/3):

Dla każdej z EOM (a)/(b)/(c):
- Wstawić ψ = a·g + b w EOM
- Rozwijać `ψ' = a·g'`, `ψ'' = a·g''`, `1/ψ = 1/(a·g+b)` etc.
- Wyrazić wszystko jako funkcję g i g'
- Porównać z R3 RHS `(1-g)/g²`
- Zbadać residuum (jeśli różnica jest mała, ID linearne pasuje;
  jeśli różnica jest istotna, trzeba bardziej ogólnej transformacji)

**Etap 2:** Substytucja **potęgowa** `ψ = g^k · h(g)` z `h(1)=1`:
- Sympy LOCK na k i h(g) takich, że EOM po podstawieniu = R3
- Jeśli k jednoznaczne (LOCK), H1/H2 aktywne
- Jeśli no LOCK, H3 wzmacniana

**Etap 3:** Substytucja **ogólna** `ψ = T(g)` z 3 constraints:
- `T(1) = 1` (vacuum match)
- `T(g₀_crit=1.874) = 4/3` (bariera↔horyzont)
- `T'(1) = 1` (slope unity dla małych perturbacji, optional)

Sympy: `solve(EOM_TGP[T(g)] − EOM_R3[g] = 0, T)`.

### PASS criterion

**Znalezienie konkretnej `T(g)` (linear, power, lub general) dającej:**
- EOM po podstawieniu redukuje się do R3 ODE z błędem RMS < 1% na zakresie
  `g ∈ [0.5, 1.7]` (cały fizyczny zakres lepton g₀ + ε margin)
- Plus zachowanie vacuum (`T(1)=1`)
- Plus zachowanie bariery (`T(1.874)=4/3` lub odwrotny mapping)

= G0b PASS.

### Plik

`phase1_G0b_field_redefinition.py` — sympy heavy, ~300 linii.

---

## 3. Sub-task G0c — Einstein-frame projection (H2, primary)

### Cel

Wykonać Weyl rescaling `g_eff^M9.1'' = Ω²(ψ)·η_μν` + scalar redefinition
`φ̃ = ∫√(K_eff(ψ))dψ`, i sprawdzić, czy w Einstein-frame lokalna EOM dla
soliton excitation reduce się do R3 ODE.

### Strategia

**Etap 1:** Wybrać Weyl factor `Ω(ψ)`:

Dla M9.1''` ds² = -c²(4-3ψ)/ψ dt² + ψ/(4-3ψ) δ_ij dx^i dx^j`:
- Konformalna płaska forma wymaga `Ω²·g_flat = g_M9.1''`
- Próby: `Ω² = ψ/(4-3ψ)` (ze spatial part) lub `Ω² = 1/√[(4-3ψ)/ψ]`
  (geometric mean), lub inne

**Etap 2:** Po Weyl rescaling akcja S_TGP staje się:

```
S_E = ∫d⁴x √(-η) [½K_E(φ̃)·η^μν ∂_μφ̃ ∂_νφ̃ - V_E(φ̃)]
```

z `K_E(φ̃)` i `V_E(φ̃)` zależnym od konkretnego Ω(ψ) i field redef.

**Etap 3:** Sympy: znaleźć takie Ω(ψ) i field redef φ̃(ψ), żeby:
- `K_E(φ̃) = φ̃⁴` (R3 kinetic prefactor)
- `V_E(φ̃) = φ̃³/3 - φ̃⁴/4` (R3 potential)

Jeśli istnieje rozwiązanie — H2 PASS. Jeśli sprzeczne równania — H2 FAIL.

### PASS criterion

**Sympy LOCK na konformalnym frame `Ω(ψ)` i field redef `φ̃(ψ)`:**
- Wybór jednoznaczny (LOCK lub nieskończenie wiele rozwiązań tej samej rodziny)
- Reprodukcja R3 Lagrangianu z błędem analitycznym = 0 (czysty sympy)
- Plus weryfikacja PPN i FRW preserved (Einstein-frame zachowuje observables)

= G0c PASS.

### Plik

`phase1_G0c_einstein_frame_projection.py` — sympy + analytical, ~300 linii.

---

## 4. Score gate Phase 1

```
Score = G0a + G0b + G0c (każdy 0 lub 1)

≥ 2/3 PASS  →  Strong evidence dla H1/H2 → Phase 2 forward
1/3 PASS    →  Weak evidence, Phase 2 z fokusem na zwycięską hipotezę
0/3 PASS    →  H3 favored → Phase 2 redirect: enumerate alternative actions
```

### Phase 2 paths

**Jeśli G0a PASS** → Phase 2: pełen sympy LOCK na S_TGP[poprawne K, V] i
weryfikacja całej why_n3 mass formuły (m_μ/m_e, m_τ/m_e, Koide).

**Jeśli G0b PASS** → Phase 2: pełen sympy LOCK na transformacji T(g),
plus interpretacja fizyczna co T(g) reprezentuje (Berry phase? frame
transformation? gauge fixing?).

**Jeśli G0c PASS** → Phase 2: pełen sympy LOCK na Einstein-frame redef
plus weryfikacja PPN γ=β=1 + cosmological κ=3/(4Φ₀) preserved.

**Jeśli 0/3 PASS** → Phase 2: enumerate alternative actions S' z
różnymi V'(ψ) lub K'(ψ) form i sprawdzić czy któraś daje (a)+(d) jednocześnie.

---

## 5. Sub-task ordering + timing

Order: **G0a → G0b → G0c** (prosty → trudny w sympy)

| Sub-task | Trudność | Czas | Wymaga |
|---|---|---|---|
| G0a | Średnia | ~3 dni | sympy + R3 numerical solver |
| G0b | Wysoka | ~4 dni | sympy LOCK + chain rule |
| G0c | Bardzo wysoka | ~5 dni | conformal transformations |

**Total Phase 1: ~12 dni** (2 tygodnie kalendarzowe).

---

## 6. Hard anchors do reprodukcji

Te liczby muszą być reprodukowane przez DOWOLNĄ wybraną redukcję:

| Wielkość | Wartość | Test wariant |
|---|---|---|
| g₀_crit (R3, α=2) | 1.874 | bariera solitonu |
| ψ_horizon (M9.1'') | 4/3 = 1.333 | Lorentzian boundary |
| g₀^e | 0.86941 | input dla R3 |
| m_μ/m_e | 206.77 | mass formuła z α=2 |
| m_τ/m_e | 3477 | Koide K=2/3 |
| PPN γ_PPN | 1.000 | M9.1'' lin. weak field |
| PPN β_PPN | 1.000 | M9.1'' "master formula" |
| κ (FRW) | 3/(4Φ₀) | po update do √(-g)=ψ/(4-3ψ) |

Każdy PASS w G0a/b/c musi być sprawdzony przeciwko subset tych anchors.

---

## 7. Plik scaffolding

```
research/op-g0-r3-from-canonical-projection/
├── README.md                          ✓ utworzony
├── Phase1_setup.md                    ✓ ten dokument
├── phase1_G0a_volume_integration.py   ← następne
├── phase1_G0b_field_redefinition.py   ← następne
├── phase1_G0c_einstein_frame_projection.py  ← następne
├── phase1_G0a_volume_integration.txt
├── phase1_G0b_field_redefinition.txt
├── phase1_G0c_einstein_frame_projection.txt
└── Phase1_results.md                  ← po wszystkich 3 sub-tasks
```

---

## 8. Status meta + ryzyka

**Status:** ACTIVE (start 2026-05-02).

### Ryzyka Phase 1

1. **Sympy może utknąć** w niektórych transformacjach (K_eff(ψ) z chain rule
   na M9.1''). Mitigation: numerical fallback dla każdej testu.

2. **R3 ODE może mieć 'drugą naturę'** — np. być EOM dla średniego pola
   po quantum integration, nie classical EOM. Wtedy żadna z H1/H2 nie zachodzi.
   Mitigation: explicit notatka w Phase1_results.md.

3. **Wszystkie 3 PASS** (overdetermined) — jeśli np. G0a i G0c oba PASS,
   musi być spójny mapping między ich rozwiązaniami. Sprawdzić w Phase 2.

### Open questions z [[README.md]] §10

Powtórzone tutaj jako tracker:

1. K(ψ) dla M9.1''? — testowane w G0a etap 1
2. A8 (`G(Φ) = G₀·Φ₀/Φ`) spójne z M9.1''? — out-of-scope Phase 1, do Phase 3
3. R3 logarithmic potential `V_R3(g) = -1/g - ln(g)` topologiczne pochodzenie?
   — Phase 2 jeśli H2 PASS

---

**Order startu sub-tasks:**
1. G0a (etap 1+2 sympy: 4 warianty K) — najpierw szybko enumerate
2. G0b (etap 1 linear T(g), potem power) — sprawdzić czy linear wystarczy
3. G0c (Einstein-frame, najtrudniejszy) — ostatni; jeśli G0a lub G0b PASS,
   może być skip

**Next step:** uruchomienie `phase1_G0a_volume_integration.py`.
