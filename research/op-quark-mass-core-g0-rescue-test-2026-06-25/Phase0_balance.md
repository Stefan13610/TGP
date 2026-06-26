---
title: "Phase 0 LOCKED — op-quark-mass-core-g0-rescue-test: rekoncyliacja HALT-B(0/5) ⟷ dodatekX(≤2,5%); tolerancje RESCUE-CONFIRMED/PARTIAL/FAILED LOCKED przed rachunkiem"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-25
cycle: op-quark-mass-core-g0-rescue-test-2026-06-25
authorization: "User 2026-06-25: 'op-quark-mass-core-g0-rescue-test działaj' — fraza aktywacyjna pokrywa Phase 0 LOCK (precedens #15-#19, #40); Phase 1+ = osobne 'działaj'. Licencja z #40 (NORM-OVERLOAD, if_recovery: 'sufit VOIDED ⇒ HALT-B reopened + licencja rescue-test')."
methodology_binding: "CALIBRATION_PROTOCOL §3.6 BINDING; CYCLE_KICKOFF §1-3; dodatekX LIVE read-only (zakaz re-litygacji R2); HALT-B verdict IMMUTABLE; PRE_REGISTERED_FALSIFIERS §0.3 append-only"
anti_lakatos_lock: PRESERVED
PR_candidate: "PR-025 (tylko przy RESCUE-CONFIRMED; LOCK = user w Phase FINAL)"
---

# Phase 0 — pre-registration LOCKED (przed jakimkolwiek rachunkiem)

## §0 — Tożsamość cyklu (nienegocjowalne)

1. **Rekoncyliacja, nie nowa derywacja.** `dodatekX` (1216 linii, v45-v47b, LIVE,
   LP-5 13/13 PASS) jest **read-only**. Re-weryfikujemy NIEZALEŻNIE **jedną** centralną
   liczbę (samozgodne m_3) + diagnozujemy strawmana HALT-B — **nie** przepisujemy dodatekX.
2. **HALT-B IMMUTABLE.** Werdykt „0/5, sufit 2,68×" pozostaje — re-scoping mówi *co* testował
   (strawman: jeden wzór z błędną domeną g₀ #40, bez m_0/φ-drabiny/Koide), nie że rachunek był zły.
3. **Pokusa ratunku = R1 META.** Tolerancje 5%/15% LOCKED PONIŻEJ; m_t pod oboma schematami;
   „zero parametrów" zakazane bez uczciwego licznika. RESCUE-PARTIAL/FAILED = pełnoprawne.
4. Każda faza = osobne „działaj".

### §0.4 — Pre-flight read confirmation
- [x] #40 (NORM-OVERLOAD; [0,817;0,891]=anchors; D1: A_tail² g₀^(e²/2) vs A_tail⁴)
- [x] HALT-B README+FINAL (strawman formula)
- [x] dodatekX CAŁY (φ-FP universal; 𝒜=a_Γ/φ=C_F²α_s²; shifted Koide; eq:X-self-consistent; LP-5 13/13; l.1207 anchors)
- [x] dodatekJ (M∝A_tail⁴; eq:J-ode)
- [x] PRE_REGISTERED_FALSIFIERS §0-2

**Sign-off:** Claudian @ 2026-06-25.

## §1 — Kontrakt L1 (observable-first)
```yaml
L1_native:
  output_observable: "werdykt rekoncyliacji + niezależne m_3 (m_b,m_t) z samozgodnego domknięcia; uczciwy licznik inputów"
  measurement_instrument: "PDG 2024 quark masses; stałe LIVE leptonowe a_Γ=0,040, φ, 𝒜=a_Γ/φ=0,02472; maszyneria dodatekX"
  native_coefs_constrained: ["𝒜=a_Γ/φ (universal, shared)", "m_0=𝒜·m_3/m_1", "shifted-Koide closure"]
  falsification_rule: "§0.2 README — RESCUE-CONFIRMED/PARTIAL/FAILED, tolerancje LOCKED"
  pre_registration_date: "2026-06-25"
```

## §2 — Pytania (CLOSED; kolejność wymuszona)

**Q1 (strawman? Phase 1):** czy formuła HALT-B `m=c_M·A_tail²·g₀^(e²/2)` z wszystkimi
g₀∈[0,817;0,891], BEZ m_0, BEZ φ-drabiny per-sektor, BEZ Koide — różni się od maszynerii
dodatekX? Zbiór składników dodatekX (CLOSED): {(i) addytywne m_0=𝒜·m_3/m_1, (ii) φ-FP
per-sektor g₀^(2)=φ·g₀^(1), (iii) shifted-Koide closure}. HALT-B ma 0/3.

**Q2 (re-weryfikacja, Phase 1):** niezależne rozwiązanie samozgodne dla down i up:
dane (m_1,m_2)+𝒜 → m_3. Porównać z PDG. Tolerancje §0.2.

**Q3 (licznik, Phase 1):** ile faktycznych inputów? per-sektor {m_1, m_2 (lub r_21 via
g₀^(1))} + shared {a_Γ, φ}. Output: m_3. Czy r_31 to predykcja (tak, jeśli m_3 nie był
inputem) czy fit?

## §3 — Falsyfikatory (CLOSED; tolerancje LOCKED 2026-06-25)

| Werdykt | Warunek (wszystkie) |
|---|---|
| **RESCUE-CONFIRMED** | strawman (≥2/3 składników brak) ∧ δ_b≤5% ∧ δ_t≤5% (≥1 schemat) ∧ δ_𝒜≤2% |
| **RESCUE-PARTIAL** | strawman ∧ domknięcie działa, ale 5%<max(δ_b,δ_t)≤15% LUB tylko 1 schemat m_t LUB r_31 słabsza predykcja (m_2 jako input) |
| **RESCUE-FAILED** | max(δ_b,δ_t)>15% w obu schematach LUB domknięcie cyrkularne (m_3=input) |

**Value-blind:** werdykt WYLICZONY z flag T2-T5; tolerancje IMMUTABLE; m_t pod oboma
schematami zaraportowany ZANIM wybrany lepszy (forbidden #5).

## §4 — Analytical pre-derivation (PRZED rachunkiem)

### §4.1 — Maszyneria dodatekX (zinwentaryzowana, LIVE)
- **φ-FP** (res:X-phiFP-universal): g₀^(2)=φ·g₀^(1) → r_21 wszystkie sektory; anchors
  {down 0,8171; lepton 0,8695; up 0,8905} = [0,817;0,891] (l.1207). [Domyka #40 NORM-OVERLOAD.]
- **𝒜=a_Γ/φ=0,02472** (thm:X-A-golden) vs 𝒜_emp 0,02464 (0,3%); 𝒜=C_F²α_s² (prop:X-A-from-tube-tension)
  → α_s=√𝒜/C_F=0,11792 vs PDG 0,1179 (0,03σ) [warunkowe na K_geo·m_sp²=π·Φ₀², CG nie domknięte — R5].
- **shifted Koide** (prop:X-shifted-koide): Q_K(m_i+m_0)=2/3; m_0(d)=21,9 MeV, m_0(u)=1981,5 MeV.
- **self-consistent** (eq:X-self-consistent): m_0=(a_Γ/φ)m_3/m_1 ∧ Q_K=2/3 → m_3 z m_1. dodatekX: m_b 0,24%, m_t 2,5%.

### §4.2 — Oczekiwanie pre-rachunkowe (INFORMATIONAL)
Niezależna re-weryfikacja POWINNA odtworzyć dodatekX (m_b ~0,2%, m_t ~2,5%) skoro maszyneria
LIVE. Jeśli NIE — sygnał, że dodatekX przeszacowany (RESCUE-FAILED), co byłoby ważnym
negatywem. Strawman HALT-B: oczekiwany 0/3 składników (formuła jawnie inna).

### §4.3 — Uczciwy licznik (anti-„zero parametrów")
Per sektor kwarkowy: m_1 (skala), r_21 (via g₀^(1), 1 param) = **2 wejścia**; m_3 = **predykcja**
(via shifted-Koide + 𝒜 shared). Shared z leptonami: a_Γ, φ. Czyli NIE „zero parametrów" —
to „r_31 predykowane z r_21 + uniwersalne 𝒜". Headline musi to oddać.

### §4.4 — Klasyfikacja stałych (§3.6.13)
| Wielkość | Klasa | Nota |
|---|---|---|
| a_Γ=0,040, φ | (α) SHARED z leptonów | nie per-sektor |
| 𝒜=a_Γ/φ=0,02472 | (α) pochodna | universal |
| m_0 | (α) pochodna | =𝒜·m_3/m_1; 0 dla leptonów |
| α_s(M_Z) | (β) SM input / lub predykcja z 𝒜 | most warunkowy R5 |
| PDG quark masses | comparison target | DOZWOLONE (to rescue-test); ale m_3 NIE jako input do własnej predykcji |
| Nowe stałe | budżet 0 | luka=GAP |

## §5 — Plan faz (re-weryfikacja first)
| Faza | Zakres | Werdykt | Gate |
|---|---|---|---|
| **Phase 1** | T1 anchors → T2 samozgodne m_3 (CENTRALNY) → T3 𝒜+α_s → T4 strawman → T5 licznik → T6 werdykt → T7 scheme-caveat | RESCUE-* | „działaj" |
| Phase FINAL | re-scope HALT-B+R12; korekta STATE WIP; D1; PR-025 (user); propagacja | — | „działaj" |

## §6 — Forbidden moves (LOCKED; 12)
1. Re-derywacja/modyfikacja dodatekX (LIVE read-only).
2. Modyfikacja werdyktu HALT-B (IMMUTABLE).
3. Dostrajanie 𝒜/m_0/schematu m_t do zgodności.
4. „Zero parametrów" bez uczciwego licznika (§4.3).
5. Wybór schematu m_t po obejrzeniu wyniku (raport pod oboma ZANIM ocena).
6. Post-hoc zmiana tolerancji 5%/15%.
7. Użycie m_3 jako inputu do własnej predykcji m_3 (cyrkularność).
8. Hardcoded T_pass.
9. Miękkie domknięcie po negatywie (RESCUE-FAILED = sukces metodologiczny).
10. Nadinterpretacja: RESCUE-CONFIRMED ≠ „cały SM z 3 inputów" — to częściowa predykcja kwarkowa.
11. Scope-creep (parameter-counting M03 = osobny cykl).
12. Append PR-025 przed Phase FINAL / bez user.

## §7 — Risk register
| ID | Ryzyko | Sev | Mitygacja |
|---|---|---|---|
| R1 | Pokusa ratunku (najsilniejsza) | META/HIGH | tolerancje LOCKED; werdykt wyliczony; licznik jawny |
| R2 | Re-litygacja dodatekX | HIGH | tylko re-weryfikacja 1 liczby; LIVE read-only |
| R3 | Cyrkularność φ-FP (r_21 input) / m_3 input | HIGH | T5 jawny licznik; forbidden #7; r_31 = prawdziwa predykcja |
| R4 | m_t scheme (pole vs MS-bar ~6%) maskuje zgodność | MED-HIGH | raport pod oboma (forbidden #5) |
| R5 | 𝒜=C_F²α_s² warunkowe (CG nie domknięte) | MED | most, NIE domknięcie; status [AN+NUM warunkowy] jawny |
| R6 | RESCUE-CONFIRMED odebrane jako triumf | MED | #10: częściowa predykcja, caveaty m_t/scheme/CG |

## §8 — Anticipated outcomes (PRZED rachunkiem)
1. **Najbardziej prawdopodobny: RESCUE-PARTIAL** — maszyneria dodatekX odtworzy m_b (~0,2%)
   i m_t (~2,5%, w MS-bar lepiej), strawman HALT-B potwierdzony (0/3 składników), ALE uczciwy
   licznik pokaże, że to „r_31 z r_21+𝒜", nie „zero parametrów", a m_t scheme-zależne. To
   uczciwie częściowa predykcja, nie pełny triumf.
2. **Możliwy RESCUE-CONFIRMED** jeśli δ_t≤5% w MS-bar i 𝒜 universal ≤2%.
3. **Pokusa #1:** ogłosić „kwarki rozwiązane, 0 parametrów". Antidotum: #10, T5 licznik.
4. **Pokusa #2:** wybrać schemat m_t dający lepszą zgodność. Antidotum: forbidden #5.
5. **RESCUE-FAILED** (mniej prawdopodobny): jeśli niezależny solver NIE odtworzy dodatekX —
   ważny negatyw (dodatekX przeszacowany), zgłoszony wprost.

## §9 — Anti-Lakatos compliance (Phase 0)
Zbiory CLOSED ✓ · tolerancje LOCKED przed rachunkiem ✓ · dodatekX LIVE read-only (zakaz
re-litygacji) ✓ · HALT-B IMMUTABLE ✓ · m_t pod oboma schematami ✓ · uczciwy licznik wymuszony
(zakaz „zero param") ✓ · 𝒜=C_F²α_s² jako most warunkowy (nie domknięcie) ✓ · nadinterpretacja
zakazana (#10) ✓ · 0 nowych stałych ✓ · anticipated outcomes z pokusami ✓ · RESCUE-FAILED =
pełnoprawny ✓.

**Phase 0 LOCKED 2026-06-25. Następny krok: Phase 1 — wymaga user „działaj".**
