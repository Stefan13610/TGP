# K-like Universals Scan — bridge §9.2 closure

**Data:** 2026-05-02
**Status:** ✅ NEGATIVE RESULT (strukturalny lock-in K-action interpretation)
**Skrypt:** [[r5_phase2_k_like_universals_scan.py]]
**Output:** [[r5_phase2_k_like_universals_scan.txt]]
**Powiązane:** [[R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02]] §9.2

---

## 0. TL;DR

> **Odkrycie:** Dla wszystkich `I(p, q, s) = ∫ g^p · |g'|^q · r^s · dr` z dobrym konwergentnym scalingiem zachodzi:
>
> $$\boxed{I(p, q, s) \sim A^q}$$
>
> **niezależnie od p, s, oraz α (substrate vs canonical).**

**Konsekwencja:**
1. **K = ∫(1/2)·g^(2α)·(g')²·r²·dr ~ A² NIE jest "topological invariant"** — to po prostu **tail derivative scaling** (g')² → A² generic.
2. **Brak szerszej klasy K-like universals.** Hipoteza bridge §9.2 (że K~A² to deeper structure) **negative closure**.
3. **Phase 2 fundamentality** leży w `g₀^[e²(1-α/4)]` factor (non-trivial g₀-dependence), NIE w A² kinematic prefactor.
4. **R5 K² = [(g')²]² to specyficzna konstrukcja** ("outer square"), nie generic (g')^4.

**Pattern fit (50 tight universals):**
$$\text{slope}_{\text{emp}} \approx 0.00046 \cdot p + 0.99966 \cdot q - 0.00112 \cdot s + 0.00184$$

Coefficient na q: 0.9997 (≈1.0). Coefficient na p, s: ~0. Intercept ~0.

---

## 1. Setup

### 1.1 Pytanie z bridge §9.2

> Czy `∫g^p·(g')^q·r^s dr ~ A²` dla innych (p, q, s)? Czy K~A² to izolowany przypadek czy klasa?

### 1.2 Test parametryczny

- **Substraty:** α=1 (R5 substrate), α=2 (TGP-canonical)
- **Lepton g₀:** g₀_e = 0.86941, g₀_μ = φ·g₀_e = 1.40674 (golden ratio)
- **Siatka:** (p, q, s) ∈ {0,1,2,3,4} × {0,1,2,3,4} × {0,1,2,3} → 100 punktów
- **Slope:** `slope_emp = log(I(g₀_μ)/I(g₀_e)) / log(A_μ/A_e)`

### 1.3 ODE setup

Soliton z R3 substrate ODE:
$$g'' = (1-g)\, g^{2-2\alpha} - \frac{\alpha}{g}(g')^2 - \frac{2}{r}\, g'$$

z BC: g(0) = g₀, g'(0) = 0. Rozwiązanie ma core (g ≈ g₀ dla r ≈ 0) i linear tail (g(r) → 1 + A·sin(r-φ)/r dla r → ∞).

---

## 2. Wyniki (CASE 1, α=1)

50/100 (p, q, s) ma **near-integer slope** (within 1%). Pełny pattern:

| q | slope_emp range | Comment |
|---:|:---|:---|
| 0 | 0.000 – 0.012 | Bez (g')²: ∫g^p·r^s ≈ const(g₀) |
| 1 | 0.999 – 1.004 | (g') daje A^1 |
| 2 | 1.991 – 2.002 | (g')² daje A^2 (R5 K case!) |
| 3 | 2.996 – 3.001 | (g')³ daje A^3 |
| 4 | 3.999 | (g')⁴ daje A^4 (sparse, low convergence) |

**Pattern jest klarny:** slope = q.

---

## 3. Wyniki (CASE 2, α=2)

Identyczny pattern:

| q | slope_emp range | Comment |
|---:|:---|:---|
| 0 | 0.000 – 0.008 | Stałe |
| 1 | 0.998 – 1.004 | A^1 |
| 2 | 1.998 – 2.000 | A^2 |
| 3 | 2.998 – 3.000 | A^3 |

**Universal:** Wszystkie tight slope_α=1 = slope_α=2 within 0.5% (50/50).

---

## 4. Linear pattern fit (CASE 4)

Dopasowanie `slope = a·p + b·q + c·s + d` na n=50 tight universals:

```
slope ~= 0.00046*p + 0.99966*q + (-0.00112)*s + 0.00184
```

Residua **|res| ≤ 0.008** dla wszystkich (max ~0.008 dla outlier p=4, q=0).

**Wniosek:** Wzór `slope = q` jest exact (z dokładnością do ODE precision ~0.005).

---

## 5. Strukturalna interpretacja

### 5.1 Skąd `slope = q`?

W tail regime soliton ma:
$$g(r) \xrightarrow{r \to \infty} 1 + A \cdot \frac{\cos(r - \phi)}{r}$$

(vacuum value 1, oscillating decay z amplitudą A_tail).

Dla potęg:
- $g^p \approx 1^p = 1$ (vacuum dominates) — **NIE skaluje się z A**
- $g'^q \approx A^q \cdot O(1/r^q)$ — **skaluje się z A^q**
- $r^s$ — pure r factor, niezmienny przy g₀ → g₀'

Więc:
$$I(p, q, s) \approx A^q \cdot \int [\text{O}(1/r^q)]^? \cdot r^s \, dr = A^q \cdot \text{const}(g_0)$$

Where const(g₀) zawiera contribution z core (gdzie g ≈ g₀) — daje numerical różnicę między g₀_e i g₀_μ ale NIE skaluje się z A.

### 5.2 Co wyklucza ten wynik?

| Pre-scan hipoteza | Status post-scan |
|---|---|
| K~A² to "topological invariant" | ✗ FALSE — to tail derivative scaling |
| Generic K-like universals (różne k) | ✗ FALSE — wszystkie q-classes są (g')^q TYPE |
| Bridge α=1 closure ma deeper structure | ✗ FALSE — Bridge α=1 to algebra, nie struktura |
| K~A² requires nonlinear core | ✗ FALSE — linear tail wystarczy (g→1, g'~A/r) |
| K~A² uniwersalne dla każdego α | ✓ TRUE — generic (g')² scaling, α-independent |

### 5.3 Co potwierdza?

| Hipoteza | Status |
|---|---|
| Phase 2 fundamentality w `g₀^n(α)` | ✓ STRENGTHENED (A² to trywialny prefactor) |
| R5 K² = c·A⁴ to algebra, nie nowa fizyka | ✓ CONFIRMED (K² = [(g')²]² = generic) |
| g₀-dependence to "real" content | ✓ STRENGTHENED |

---

## 6. Implikacje dla TGP architecture

### 6.1 Hierarchia "fundamental" → "derivative"

**Pre-scan:**
```
Phase 2 m_obs(g₀, α) = c·A²·g₀^n(α)
                       ↑    ↑
                    A²?    g₀^n(α)?
                    Co jest fundamental?
```

**Post-scan:**
```
A² = (g')² tail derivative scaling
   = GENERIC kinematic prefactor (uniwersalny dla any q=2 integral)

g₀^n(α) = NON-TRIVIAL g₀-dependence
        = TGP-FUNDAMENTAL content (encoded n(α) = e²(1-α/4))
```

### 6.2 Falsyfikacja R5 "topological" interpretation

**R5 README claim (pre-scan):**
> "K = ∫(1/2)·g^(2α)·(g')²·r²·dr scales as A² universally — topological invariant"

**Post-scan correction:**
> "K is one specific (g')² integral; ALL ∫g^p·|g'|²·r^s·dr scale as A². 'Universal' in the trivial sense (kinematic), nie w 'topological' sensie. Strukturalny content R5 leży gdzie indziej (np. w c_M calibration, jeśli istnieje)."

**Re-framing R5:** R5 K-action jest jedną konkretną reprezentacją Lagrangianu (∝ kinetic energy density), nie "topological structure". Slope=2 nie wyróżnia K spośród inne (g')² integrals.

### 6.3 Phase 2 fundamentality wzmocniona

n(α) = e²(1-α/4) jest jedyną nietryvialną częścią mass formula.

- **e² = exp(2) = Euler²** — verified via μ/e exact fit do 0.0007% (sub-tension closure 2026-05-01)
- **(1-α/4) factor** — connects do dimensional analysis (Hobart-Derrick, X=e²/4)
- **g₀^n(α) factor** — to GDZIE jest "fizyka"; A² to tylko kinematic prefactor

Co potwierdza: **fundamental question to derivacja n(α), NIE A² scaling.** X=e²/4 RG-derivation (osobny agent) jest właściwym kierunkiem.

---

## 7. Otwarte (po tym scan)

Pytania ZAMKNIĘTE przez ten scan:
- ✅ Czy istnieje klasa K-like universals? **NIE** (wszystkie redukują do `slope=q`).
- ✅ Czy K~A² to topological? **NIE** (kinematic).
- ✅ Czy Phase 2 A² ma deeper structure? **NIE** (trywialne).

Pytania OTWARTE:
- 🔓 **Derywacja n(α) = e²(1-α/4) z RG/dimensional analysis** (osobny agent, X=e²/4)
- 🔓 **C_core ≈ 1.09 zamknięta forma** (R5 closure)
- 🔓 **m_obs vs M_full distinction** (Hobart-Derrick / ADM-Komar)
- 🔓 **ψ(g, r) full nonlinear profile** (Phase 1 follow-up)

---

## 8. Pliki

- **[[r5_phase2_k_like_universals_scan.py]]** — skrypt z 5 cases (full scan + pattern fit + R5 verify)
- **[[r5_phase2_k_like_universals_scan.txt]]** — output run
- **[[R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02]]** — bridge document (§9.2 source)

---

## 9. Status

| Pytanie | Pre-scan | Post-scan |
|---|---|---|
| Czy istnieje klasa K-like universals? | open | **NEGATIVE** — wszystkie kinematic |
| Czy K~A² to "topological invariant"? | claim | **FALSE** — generic tail (g')² |
| Czy Phase 2 A² to fundamental? | unclear | **NIE** — kinematic prefactor |
| Czy g₀^n(α) to fundamental? | claim | **YES + STRENGTHENED** |

✅ **NEGATIVE RESULT — strukturalny lock-in: K-action to NIE topological structure, lecz tail derivative scaling.**

To **zacieśnia** fundamental question TGP do n(α) derivation (X=e²/4) — Phase 2 fundamentality jest w g₀-dependence, nie w A² kinematics.
