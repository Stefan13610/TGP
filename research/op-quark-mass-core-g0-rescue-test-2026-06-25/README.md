---
title: "op-quark-mass-core-g0-rescue-test — REKONCYLIACJA: HALT-B (2026-05-16, '0/5 structural insufficiency, sufit 2,68× vs 80000×') ⟷ dodatekX (2026-04-05, LIVE: masy kwarków odtworzone ≤2,5% przez φ-FP + addytywne m_0=a_Γ/φ·r_31 + shifted Koide). Czy sektor kwarkowy jest 'uratowany' (HALT-B testował strawmana), czy claimy dodatekX nie wytrzymują niezależnej weryfikacji?"
date: 2026-06-25
type: research-cycle
folder_status: parking
parent: "[[../../audyt/L08_kink_fermion_closure/README.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  L1_native:
    output_observable: "Werdykt rekoncyliacji sprzeczności HALT-B vs dodatekX + niezależna re-weryfikacja centralnej predykcji kwarkowej: m_3 (m_b, m_t) z samozgodnego domknięcia {m_0=𝒜·m_3/m_1 ; shifted-Koide(m_i+m_0)=2/3}, 𝒜=a_Γ/φ uniwersalne (z sektora leptonowego). Wielkość pochodna: uczciwy licznik inputów sektora kwarkowego (czy r_31 jest PREDYKCJĄ czy fitem)."
    measurement_instrument: "PDG 2024 (MS-bar 2 GeV / m(m) / pole): m_u=2,16 m_d=4,67 m_s=93,4 m_c=1270 m_b=4180 m_t=172760 MeV. Stałe LIVE z sektora leptonowego: a_Γ=0,040, φ=(1+√5)/2, 𝒜=a_Γ/φ=0,02472. Maszyneria LIVE: dodatekX (φ-FP res:X-phiFP-universal; thm:X-A-golden 𝒜=a_Γ/φ; prop:X-shifted-koide; LP-5 13/13 PASS)."
    native_coefs_constrained:
      - "𝒜 = a_Γ/φ (uniwersalna stała konfinementu; SHARED z leptonami; NIE per-sektor wolna)"
      - "m_0^(sektor) = 𝒜·m_3/m_1 (addytywna masa rury kolorowej; m_0=0 dla leptonów)"
      - "samozgodne domknięcie: {m_0-relacja ∧ shifted-Koide=2/3} → m_3 z (m_1,m_2)"
      - "φ-FP base anchors {down 0,817; lepton 0,870; up 0,891} = sek08b:529 [0,817;0,891] (potwierdza NORM-OVERLOAD #40)"
    falsification_rule: "§0.2 — RESCUE-CONFIRMED / RESCUE-PARTIAL / RESCUE-FAILED, tolerancje LOCKED PRZED rachunkiem."
    pre_registration_date: "2026-06-25"
  L2_framework_reduction:
    target_frameworks:
      - "Rekoncyliacja dwóch wewnętrznych warstw TGP: cykl HALT-B (audyt) ⟷ dodatekX (manuskrypt LIVE)"
      - "𝒜 = C_F²·α_s²(M_Z) (most do α_s; prop:X-A-from-tube-tension)"
    reduction_type: "internal-reconciliation + independent re-verification (NIE nowa derywacja — zakaz re-litygacji dodatekX)"
    failure_disposition: "L1-stands"
  L3_falsification_map:
    - { bound: "m_b PDG = 4180 MeV", constrains: "samozgodne domknięcie down z (m_d,m_s,𝒜)", window: "≤5% RESCUE / 5-15% PARTIAL / >15% FAILED", status: "pending Phase 1" }
    - { bound: "m_t PDG = 172760 MeV (pole)", constrains: "samozgodne domknięcie up z (m_u,m_c,𝒜)", window: "≤5% / 5-15% / >15% (uwaga: pole vs MS-bar ~6% scheme)", status: "pending Phase 1" }
    - { bound: "𝒜_emp (down,up) = 0,0245/0,0248", constrains: "𝒜=a_Γ/φ=0,02472 universal", window: "≤2%", status: "pending Phase 1" }
    - { bound: "α_s(M_Z) PDG = 0,1179±0,0009", constrains: "𝒜=C_F²α_s² → α_s=√𝒜/C_F", window: "≤1σ", status: "pending Phase 1" }
    - { bound: "HALT-B verdict IMMUTABLE", constrains: "werdykt NIETYKALNY; re-scoping dotyczy ZAKRESU (co testował), nie poprawności dla testowanej formuły", window: "anti-Lakatos invariant", status: "LOCKED" }
    - { bound: "dodatekX status LIVE [AN+NUM]/[POST+NUM+UNIV]", constrains: "treść manuskryptu read-only; re-weryfikacja niezależna, NIE re-derywacja", window: "anti-Lakatos invariant", status: "LOCKED" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: L1
  kind: reconciliation-audit + independent-verification
  output_type: observable (m_3 prediction) + epistemic (input accounting + HALT-B re-scope)
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: true
  depends_on:
    - "op-L08-quark-g0-tail-vs-core-audit-2026-06-25 (#40 NORM-OVERLOAD: [0,817;0,891]≠domena core g₀; T11 VOIDED)"
    - "op-L08-Phase6-quark-sector-mass-formula-2026-05-16 (HALT-B; strawman do zdiagnozowania)"
    - "partial_proofs/quark_sector/dodatekX (φ-FP universal; 𝒜=a_Γ/φ=C_F²α_s²; shifted Koide; LP-5 13/13 PASS)"
    - "partial_proofs/hierarchia_mas/dodatekJ (M∝A_tail⁴; eq:J-ode; φ-FP lepton chain)"
  impacts:
    - "audyt/L08 problem #3 (quarks) — status: INDETERMINATE-PENDING-RESCUE → rozstrzygnięcie"
    - "STATE.md WIP 'quark-mass HALT-B (sufit 2,68× vs 80000×)' — potencjalnie fałszywie pesymistyczny status do korekty"
    - "core/sek07_predykcje R12 (m_b/m_t recovery) — status update"
    - "PREDICTIONS_REGISTRY — sektor kwarkowy: re-klasyfikacja (jeśli RESCUE-CONFIRMED)"
    - "DOUBT D1 (#40): m=c_M·A_tail²·g₀^(e²/2) [HALT-B] vs M=c_M·A_tail⁴ [dodatekJ] — rozstrzygnięcie który wzór kanoniczny"

predecessors:
  - "[[../op-L08-quark-g0-tail-vs-core-audit-2026-06-25/]] (#40 NORM-OVERLOAD; licencja tego cyklu)"
  - "[[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/]] (HALT-B strawman)"
  - "[[../op-FFS-quark-object-2026-05-20/]] (kwark = FFS z ogonem; ładunki ułamkowe)"

related:
  - "[[../../partial_proofs/quark_sector/dodatekX_quark_sector.tex]] (maszyneria kwarkowa LIVE)"
  - "[[../../partial_proofs/hierarchia_mas/dodatekJ_ogon_masy.tex]] (M∝A_tail⁴; φ-FP)"
  - "[[../../core/sek07_predykcje/sek07_predykcje.tex]] R12 (m_b/m_t recovery)"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] (PR-025 candidate)"
  - "[[../op-theta-quark-koide/]] (wcześniejsze prace quark-Koide)"

classification: RECONCILIATION-AUDIT — rozstrzygnięcie sprzeczności HALT-B(0/5) vs dodatekX(≤2,5%) + niezależna re-weryfikacja predykcji kwarkowej + uczciwy licznik inputów
priority: high (rozstrzyga, czy STATE WIP 'quark-mass HALT-B' to fałszywie pesymistyczny status; licencja z #40)
goal: "Value-blind: (1) zdiagnozować, że HALT-B testował STRAWMANA (formuła bez m_0/φ-drabiny/Koide, z błędną domeną g₀ #40); (2) NIEZALEŻNIE re-zweryfikować centralną predykcję dodatekX — samozgodne domknięcie {m_0=𝒜·m_3/m_1; shifted-Koide=2/3} z 𝒜=a_Γ/φ → m_b,m_t vs PDG; (3) uczciwie policzyć inputy (czy r_31 to PREDYKCJA czy fit). Werdykt: RESCUE-CONFIRMED / RESCUE-PARTIAL / RESCUE-FAILED."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 re-weryfikacja sympy + Phase FINAL)"

risk_flags:
  - "R1 (META/HIGH): pokusa ratunku — najsilniejsza. Antidotum: tolerancje LOCKED §0.2 przed rachunkiem; werdykt WYLICZONY; m_t scheme-caveat jawny; uczciwy licznik inputów (NIE 'zero parametrów' jeśli per-sektor potrzeba m_1,m_2)."
  - "R2 (anti-re-litygacja): dodatekX jest LIVE — NIE re-derywujemy, tylko NIEZALEŻNIE re-weryfikujemy 1 centralną liczbę (m_3) + diagnozujemy strawmana HALT-B."
  - "R3 (input honesty): φ-FP 'reprodukuje r_21' jest częściowo cyrkularne (1 param g₀^(1) per sektor ustawia r_21). Prawdziwa predykcja = r_31 z (r_21 + 𝒜). Policzyć jawnie."
  - "R4: m_t convention (pole vs MS-bar ~6%) może maskować/poprawiać zgodność — pre-rejestrować, który schemat; raport pod oboma."
  - "R5: 𝒜=C_F²α_s² wymaga identyfikacji K_geo·m_sp²=π·Φ₀² (CG-1/CG-3 nie domknięte) — to most, nie domknięcie; status [AN+NUM warunkowy]. NIE liczyć jako pełnej derywacji."

phase_plan:
  Phase_0: "Balance + 8/8 gate + tolerancje LOCKED + forbidden moves + pre-derywacja (maszyneria dodatekX zinwentaryzowana)"
  Phase_1: "Sympy: T1 φ-FP anchors=[0,817;0,891] (zamknięcie z #40); T2 NIEZALEŻNE samozgodne m_3 (m_b,m_t) vs PDG; T3 𝒜=a_Γ/φ universal + 𝒜=C_F²α_s²→α_s; T4 diagnoza strawmana HALT-B (3 brakujące składniki); T5 uczciwy licznik inputów (r_31 predykcja vs fit); T6 werdykt wyliczony; T7 scheme-caveat"
  Phase_FINAL: "Werdykt RESCUE-*; re-scope HALT-B + R12; korekta STATE WIP; rozstrzygnięcie D1; PR-025 (jeśli CONFIRMED, user-gated); propagacja"

tags:
  - L08
  - quark-sector
  - rescue-test
  - reconciliation-HALT-B-vs-dodatekX
  - phi-FP
  - shifted-koide
  - confinement-m0
  - A-equals-aGamma-over-phi
  - input-accounting
  - value-blind
  - anti-Lakatos-LOCKED
  - cycle-scaffold-2026-06-25
---

# op-quark-mass-core-g0-rescue-test-2026-06-25

> **Cel (licencja z #40 NORM-OVERLOAD):** rozstrzygnąć **sprzeczność wewnętrzną**
> między cyklem HALT-B (2026-05-16: „sektor kwarkowy structural insufficiency, 0/5,
> sufit 2,68× vs 80000×") a **żywym manuskryptem `dodatekX`** (2026-04-05: masy kwarków
> odtworzone do ≤2,5% przez **φ-FP + addytywne m_0 = a_Γ/φ·(m_3/m_1) + shifted Koide**).
> Czy sektor jest **„uratowany"** (HALT-B testował strawmana), czy claimy dodatekX nie
> wytrzymują **niezależnej** re-weryfikacji? Werdykt **RESCUE-CONFIRMED / PARTIAL / FAILED**.

## §0 — Cel + native-first contract

### §0.0 — Tożsamość cyklu (nienegocjowalne)

1. **AUDYT REKONCYLIACYJNY, NIE NOWA DERYWACJA.** `dodatekX` jest **żywą treścią
   manuskryptu** (status [AN+NUM], LP-5 „[POST+NUM+UNIV]" 13/13 PASS) — **read-only,
   zakaz re-litygacji** (R2). Werdykt **HALT-B IMMUTABLE** — re-scoping dotyczy *co* testował
   (formuła), nie poprawności rachunku dla tej formuły.
2. **Rola cyklu = niezależna weryfikacja + uczciwy licznik.** (a) Zdiagnozować, że HALT-B
   testował **strawmana** (formuła bez m_0/φ-drabiny/Koide + błędna domena g₀ z #40);
   (b) **niezależnie** odtworzyć 1 centralną predykcję dodatekX (samozgodne m_3); (c) policzyć
   uczciwie inputy (czy r_31 to predykcja czy fit).
3. **Pokusa ratunku = najsilniejsze ryzyko META (R1).** Tolerancje LOCKED PRZED rachunkiem;
   m_t scheme-caveat jawny; „zero parametrów" NIE wolno głosić, jeśli per-sektor potrzeba
   (m_1, m_2). RESCUE-PARTIAL/FAILED to pełnoprawne wyjścia.
4. Każda faza = osobne „działaj".

### §0.1 — Native observable target

**Centralna re-weryfikacja (niezależna od dodatekX):** dla każdego sektora kwarkowego,
mając (m_1, m_2) [PDG] + uniwersalne 𝒜 = a_Γ/φ [z leptonów], rozwiązać samozgodnie:
```
m_0 = 𝒜 · m_3 / m_1            (addytywna masa rury kolorowej)
Q_K(m_1+m_0, m_2+m_0, m_3+m_0) = 2/3   (shifted Koide)
```
→ predykcja **m_3** (m_b dla down, m_t dla up) → porównać z PDG. To genuine predykcja
(r_31 z r_21 + 𝒜), NIE fit m_3.

### §0.2 — Pre-registered falsification rule (LOCKED 2026-06-25, PRZED rachunkiem)

> Niech `δ_b = |m_b^pred − m_b^PDG|/m_b^PDG`, `δ_t = |m_t^pred − m_t^PDG|/m_t^PDG`
> (m_t: raport pod OBOMA schematami — pole 172760 i MS-bar ~162500). Niech `δ_𝒜` =
> odchylenie 𝒜=a_Γ/φ od 𝒜_emp; `δ_αs` = odchylenie α_s=√𝒜/C_F od PDG (w σ).
>
> - **RESCUE-CONFIRMED** ⟺ (a) **strawman zdiagnozowany** (formuła HALT-B udowodniona jako
>   ≠ maszyneria dodatekX: brak ≥2 z 3 składników {addytywne m_0, φ-drabina per-sektor,
>   shifted-Koide}), **ORAZ** (b) niezależne samozgodne domknięcie daje `δ_b ≤ 5%` **I**
>   `δ_t ≤ 5%` (w co najmniej jednym schemacie m_t), **ORAZ** (c) 𝒜=a_Γ/φ uniwersalne
>   `δ_𝒜 ≤ 2%`. ⇒ HALT-B re-scoped: testował misformułowanie; sektor kwarkowy reprodukuje
>   spektrum. R12 zawężony; PR-025 candidate.
> - **RESCUE-PARTIAL** ⟺ strawman zdiagnozowany ORAZ domknięcie działa, ale z istotnym
>   zastrzeżeniem: `5% < max(δ_b,δ_t) ≤ 15%` LUB zgodność tylko w jednym schemacie m_t
>   LUB licznik inputów ujawnia, że r_31 jest słabszą predykcją niż leptonowe (np. m_2
>   potrzebne jako input). ⇒ sektor kwarkowy CZĘŚCIOWO predyktywny; HALT-B re-scoped z caveatami.
> - **RESCUE-FAILED** ⟺ niezależna re-weryfikacja NIE odtwarza dodatekX (`max(δ_b,δ_t) > 15%`
>   w obu schematach) LUB domknięcie okazuje się cyrkularne (m_3 faktycznie input). ⇒ claimy
>   dodatekX przeszacowane; HALT (status quo) utrzymany; rozbieżność dodatekX↔HALT-B = otwarty problem.

```yaml
pre_registration_date: 2026-06-25
recovery_scope:
  allowed_directions:
    - "Niezależna re-weryfikacja 1 centralnej liczby (m_3) z LIVE maszynerii dodatekX"
    - "Diagnoza strawmana HALT-B (porównanie formuł, bez modyfikacji werdyktu)"
    - "Uczciwy licznik inputów (per-sektor m_1,m_2 + shared a_Γ,φ); raport m_t pod oboma schematami"
  forbidden_directions:
    - "Re-derywacja/modyfikacja dodatekX (LIVE; read-only)"
    - "Modyfikacja werdyktu HALT-B (IMMUTABLE)"
    - "Dostrajanie 𝒜, m_0, schematu m_t do zgodności (value-blind violation)"
    - "Głoszenie 'zero parametrów' bez uczciwego licznika"
    - "Post-hoc zmiana tolerancji 5%/15%"
    - "Nowe pola/stałe (S05; budżet 0)"
  if_recovery_exhausted:
    - "RESCUE-FAILED: HALT status quo + dodatekX↔HALT-B rozbieżność = osobny otwarty problem"
```

### §0.3 — TGP-native check (mandatory)

- [x] **Q1:** Pattern 2.7 (A_tail asymptotic); maszyneria LIVE dodatekX.
- [x] **Q2:** NONE — brak BD-form; klasyczne masy solitonu + addytywna energia rury.
- [x] **Q3 (Inherited LOCKs):** #40 NORM-OVERLOAD ([0,817;0,891]=anchors); dodatekX LP-5 LIVE; HALT-B IMMUTABLE; lepton φ-FP+Koide LIVE.
- [x] **Q4:** sympy nsolve samozgodnego układu; Q_K algebra; universal.
- [x] **Q5:** N/A (klasyczne masy).
- [x] **Q6:** N/A.
- [x] **Q7 (ASK-RULE):** „masa kwarka" = m_substrate(ogon solitonu) + m_0(energia rury kolorowej, reżim III) — declarative; m_0=0 dla leptonów (brak koloru).
- [x] **Q8:** self-audit Phase FINAL.

### §0.4 — Pre-flight read confirmation

- [x] #40 cały cykl (NORM-OVERLOAD; [0,817;0,891]=anchors; D1)
- [x] HALT-B README + Phase_FINAL (strawman: m=c_M·A_tail²·g₀^(e²/2), wszystkie g₀∈I, bez m_0/Koide)
- [x] **dodatekX CAŁY** (1216 l.): φ-FP universal (res:X-phiFP-universal), Koide leptonowy-specyficzny (res:X-koide-failure), 𝒜=a_Γ/φ=C_F²α_s² (thm:X-A-golden, prop:X-A-from-tube-tension), shifted Koide (prop:X-shifted-koide), self-consistent (eq:X-self-consistent), LP-5 13/13 PASS (rem:X-m0-numerical); linia 1207 potwierdza [0,817;0,891]=anchors
- [x] dodatekJ (M∝A_tail⁴; eq:J-ode g₀=warunek brzegowy)
- [x] PRE_REGISTERED_FALSIFIERS §0-§2 (PR-025 candidate)

**Sign-off:** Claudian @ 2026-06-25.

### §0.5 — Sympy plan (Phase 1)

| Test | Klasa | Pytanie | PASS |
|---|---|---|---|
| T1 | FP | φ-FP base anchors {down,lepton,up}={0,817;0,870;0,891}=sek08b:529 I (zamknięcie #40) | odtworzenie zbioru |
| T2 | FP (CENTRALNY) | NIEZALEŻNE samozgodne {m_0=𝒜·m_3/m_1; Q_K=2/3} → m_b,m_t vs PDG | δ_b,δ_t vs tolerancje |
| T3 | FP | 𝒜=a_Γ/φ vs 𝒜_emp(down,up); 𝒜=C_F²α_s² → α_s=√𝒜/C_F vs PDG | δ_𝒜≤2%, δ_αs≤1σ |
| T4 | FP | diagnoza strawmana: HALT-B formuła vs dodatekX (3 brakujące składniki) | ≥2/3 brakuje |
| T5 | FP | uczciwy licznik inputów: per-sektor (m_1,m_2) + shared (a_Γ,φ) → m_3 predykcja | jawny bilans |
| T6 | FP | werdykt RESCUE-* wyliczony z T2-T5 | determinowany |
| T7 | LIT/DEC | scheme-caveat m_t (pole vs MS-bar); S05; HALT-B/dodatekX IMMUTABLE | raport |

**Substance:** 6 FP + 1 LIT/DEC. 0 hardcoded. Werdykt wyliczony.

---

## §1 — Phase 0
[Patrz `Phase0_balance.md` — LOCKED]

## §2 — Phase 1 [wymaga „działaj"]
## §FINAL — Werdykt RESCUE-*

---

## Status
🟡 **PARKING — scaffold + Phase 0 LOCK 2026-06-25.** Autoryzacja: user „op-quark-mass-core-g0-rescue-test działaj". **Phase 1 wymaga osobnego „działaj".**

**Cross-references:**
- [[../op-L08-quark-g0-tail-vs-core-audit-2026-06-25/Phase_FINAL_close.md]] (#40 NORM-OVERLOAD)
- [[../../partial_proofs/quark_sector/dodatekX_quark_sector.tex]] (maszyneria LIVE)
- [[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/Phase_FINAL_close.md]] (HALT-B strawman)
