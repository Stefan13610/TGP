---
title: "Phase 0 — pre-registration LOCK: pulsar orbital-decay energy balance (PR-025 candidate)"
date: 2026-06-13
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-13 (before any computation)
cycle: op-PSR-Pdot-energy-balance-2026-06-13
authorization: "User 2026-06-13: 'ok, sprawdźmy to ;)' (po dyskusji trylemat energetyczny)"
origin_discussion: "sesja 2026-06-13 — ekspert external review; user wyprowadził: dipol znika (C_i=q·m_i), monopol znika (Q̇=0), skalar kwadrupolowy P_φ=(G/30)⟨I⃛²⟩=(1/6)P_GR, brak cutoffu m_sp~H₀, σ-channel h_TT^σ=h_TT^GR przy c₀ξ_eff=16πGΦ₀²"
tags:
  - pre-registration
  - PR-025-candidate
  - pulsar-Pdot
  - energy-balance
  - anti-Lakatos-LOCKED
---

# Phase 0 — pre-registration LOCK

## §1 — Pytanie badawcze

**Czy sektor radiacyjny TGP (skalar Φ + kanał σ_ab) reprodukuje obserwowany zanik
orbity podwójnego pulsara J0737−3039 na poziomie bilansu ENERGII (nie amplitudy)?**

Kluczowe rozróżnienie ustalone w dyskusji 2026-06-13: zamknięty wynik
`op-sigma-3PN-radiative` (h_TT^σ/h_TT^GR = 1) jest twierdzeniem o **amplitudzie
odpowiedzi detektora**. Ṗ_b pulsara mierzy **bilans energii** (dE/dt przez sferę
w nieskończoności). W repo nie istnieje żaden rachunek strumienia energii
(grep: Isaacson / energy flux / dE/dt / Ṗ_b → 0 trafień w cyklach radiacyjnych).

## §2 — Wejścia LOCKED (wszystkie sprzed 2026-06-13, źródła w repo)

| Input | Wartość | Źródło LOCK |
|---|---|---|
| Sprzężenie skalarne | C_i = q·m_i uniwersalne; q² kalibrowane limitem Newtona (do wyprowadzenia w T2a, oczekiwane q²=4πG) | S05 + Phase 5 `op-emergent-metric-from-interaction` (G_eff = q²/4πΦ₀²K₁) |
| Masa skalara kosmologicznego | m_sp ~ H₀/c (bezmasowy na skalach orbitalnych) | Appendix E / op-G-substrate-derivation-2026-05-24 |
| Masa kanału σ | m_σ ≈ 0.71 meV | `op-sigma-3PN-radiative` Phase 3 needs_blocker + `op-sigma-yukawa-audit-2026-05-09` |
| Lagrangian σ (Path A) | L_σ = −¼(∂σ_ab)² − ½m_σ²σ² − (ξ_eff/2)σ_ab T^{ab,TT} | `op-sigma-3PN-radiative` Phase 3 §1.1 (OP-7 T3.1 LOCK) |
| Normalizacja amplitudy | c₀·ξ_eff = 16πGΦ₀² (⟹ h_TT^σ = h_TT^GR) | T3.4 amendment, FOUNDATIONS §3.6.10.6 |
| P_GR (kwadrupol OTW) | P_GR = (G/5c⁵)⟨I⃛_ij I⃛_ij⟩ | textbook anchor (Maggiore §3) — LOCKED external |
| Status LIVE mechanizmu σ-mediacji | P6 CONDITIONAL; mechanizm (iii) FAILS przy m_Φ ~ M_Pl (LIGO scales) | FOUNDATIONS §3.6.9 CL-1 annotation 2026-06-01 |

**Dane obserwacyjne LOCKED (literatura, nie repo):**

| System | Wielkość | Wartość | Źródło |
|---|---|---|---|
| PSR J0737−3039A/B | R_obs ≡ Ṗ_b^int/Ṗ_b^GR | 0.999963 ± 0.000063 | Kramer et al., PRX 11, 041050 (2021) |
| PSR J0737−3039A/B | P_b | 0.10225156248 d | jw. |
| PSR J0737−3039A/B | e | 0.0877775 | jw. |
| PSR B1913+16 (cross-check) | R_obs | 0.9983 ± 0.0016 | Weisberg & Huang, ApJ 829, 55 (2016) |

## §3 — Gałęzie pre-deklarowane (WSZYSTKIE przed rachunkiem; anti-Lakatos)

- **Gałąź A (LIVE default):** m_σ = 0.71 meV jak w audycie. Pole masywne nie
  propaguje przy ω < m_σ/ℏ ⟹ kanał σ martwy radiacyjnie (ℏω_GW ~ 10⁻¹⁸ eV).
  Energia ucieka wyłącznie skalarem. Oczekiwane R = P_φ/P_GR.
- **Gałąź B (mechanizm-v hipotetyczny):** σ efektywnie bezmasowy u detektora.
  Strumień energii κ_E·P_GR, gdzie κ_E wynika z normalizacji kinetycznej σ
  — UWAGA pre-deklarowana: warunek amplitudowy (λξ_eff = 16πG) i warunek
  strumienia (κ_E = 1) są NIEZALEŻNE; amplitude-lock NIE pinuje energii.
  R = P_φ/P_GR + κ_E.
- **Gałąź C (exact match):** wymaga skasowania kanału skalarnego przy pełnym
  sprzężeniu newtonowskim — brak mechanizmu w LOCKED frameworku.
- **Gałąź D (screening skalara przy źródle):** samo-destrukcyjna strukturalnie:
  skalar niesie CAŁĄ statykę Newtona (σ przy 0.71 meV ma zasięg ℏc/m_σ ≈ 0.28 mm);
  ekranowanie radiacji ekranuje też statykę ⟹ brak wiązania orbity ⟹ sprzeczność.

## §4 — Falsyfikatory LOCKED

| ID | Test | Reguła decyzyjna (IMMUTABLE) |
|---|---|---|
| **F-PDOT-A** | Współczynnik skalarny: wyprowadzić P_φ dla binarki kołowej ab initio (całki kątowe + uśrednienie po okresie, sympy-exact, zero fitów). | Wynik = wyrażenie zamknięte; porównać z claimem usera 1/6·P_GR. Jakikolwiek hardcode = HALT. |
| **F-PDOT-B** | Bramka propagacji σ: ℏω_GW vs m_σ c² dla J0737−3039. | ℏω_GW < m_σ c² ⟹ Gałąź A obowiązuje (kanał σ nie unosi energii). |
| **F-PDOT-C** | Bilans: R_TGP = Ṗ_b^TGP/Ṗ_b^GR per gałąź vs R_obs. | **FALSIFIED-BRANCH jeśli \|R_TGP − R_obs\| > 5σ_obs.** TGP-sektor PASS tylko jeśli JAKAŚ gałąź oparta WYŁĄCZNIE na wartościach LOCKED (§2) daje \|R−R_obs\| ≤ 5σ. Wartość κ_E dobrana post-hoc do danych = FAIL (anti-Lakatos). |
| **F-PDOT-D** | Robustność e: poprawki ekscentryczności (e = 0.0878) do STOSUNKU kanałów. | Jeśli poprawka O(e²) może zmienić werdykt F-PDOT-C ⟹ DOWNGRADE do INCONCLUSIVE. |

## §5 — Pre-deklarowane oczekiwanie (honest, wzór PR-004 Phase 0 §4)

Oczekuję (ekspert zewnętrzny, przed rachunkiem):
- F-PDOT-A: P_φ = (16/15)Gμ²d⁴ω⁶/c⁵ ⟹ P_φ/P_GR = 1/6 EXACT (zgodnie z claimem usera),
- F-PDOT-B: NON_PROPAGATING (ω/m_σ ~ 10⁻¹⁵),
- F-PDOT-C Gałąź A: R = 1/6 ⟹ odchylenie ~1.3×10⁴σ ⟹ FALSIFIED-BRANCH,
- F-PDOT-C Gałąź B: κ_E nieokreślone przez LOCK; wariant kanoniczny ξ=λ ⟹ κ_E=1 ⟹ R=7/6 ⟹ ~2.6×10³σ ⟹ FALSIFIED-BRANCH; PASS wymagałby κ_E = 0.8333 ± 0.0003 jako NOWEJ kotwicy strojonej daną ⟹ forbidden.

Jeśli rachunek da co innego niż oczekiwanie — wynik rachunku nadrzędny.

## §6 — Forbidden moves

1. Zakaz modyfikacji m_σ, q², c₀, ξ_eff po LOCK.
2. Zakaz nowych stałych (budżet 0) i fitowania czegokolwiek do R_obs.
3. Zakaz post-hoc reinterpretacji, który kanał "naprawdę" unosi energię.
4. Zakaz zmiany progów 5σ ex post.
5. Eccentricity: wolno tylko OSZACOWAĆ wpływ na stosunek (bound), nie fitować.
6. Wynik FALSIFIED nie uruchamia recovery w tym cyklu (osobny scoping, wzór PR-004 → galactic-substrate-tail).
