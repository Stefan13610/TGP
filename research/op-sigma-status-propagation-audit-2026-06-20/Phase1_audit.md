---
title: "Phase 1 AUDIT — tabela trafień + klasyfikacja"
date: 2026-06-22
type: phase_audit
status: 🟢 DONE (klasyfikacja wyliczona z reguły Phase0 §3)
phase: 1
cycle: op-sigma-status-propagation-audit-2026-06-20
---

# Phase 1 — AUDYT (tabela trafień, klasyfikacja wyliczona)

Metoda: grep wszystkich celów (§4 LOCK) na {sigma, σ_ab, κ_E, C_σ, tensor, GW, breathing,
polaryzacja, radiacyjn, „wolny parametr", „2 parametry", derived/wyprowadz/predict/zamknięty} +
odczyt kontekstu. Każdy werdykt = reguła §3 LOCK.

## §A — SPÓJNE (bez zmian)

| # | Plik / kotwica | Treść | Podstawa SPÓJNE |
|---|---|---|---|
| S1 | `PREDICTIONS_REGISTRY.md:1266` GW7 | „C_σ (≡κ_E) = nieredukowalny WOLNY PARAMETR", 2 dowody, param +1, M911-* warunkowe | = F wprost (zaktualizowane #34) |
| S2 | `PREDICTIONS_REGISTRY.md:1260–1265` GW1–GW6 | c_T=c_s, 3 DOF, h_b=h_L, m_σ²/m_s²=2 LOCKED, no-vector, masslessness | R1+R2+R4 — niezależne od κ_E |
| S3 | `core/sek08_formalizm.tex:6989` rem:sigma-Csigma-free | 2 dowody UV+brak bieguna; M²=2m_s²=OPE; M911-* warunkowe | = F wprost (cel docelowy spójności) |
| S4 | `core/sek08_formalizm.tex:6959` rem:sigma-params | „C_σ dowiedzenie nieredukowalny param swobodny; bilans 2+1=3" | = F wprost (zaktualizowane #34) |
| S5 | `core/sek08_formalizm.tex:6853` thm:amplitude-matching | ξ_eff=4πG₀σ₀Φ₀ jako **warunek dopasowania** | R3 — matching, nie predykcja |
| S6 | `core/sek08_formalizm.tex:7030` rem:sigma-logic (tabela) | σ_ab amplituda = „Propozycja", masa/modes = „Hipoteza" | status pokorny, brak overclaim |
| S7 | `core/formalizm/dodatekQ_…tex:219–231` Q.5 | „residual C_σ… nie jest problemem substratu: C_σ jest [wolny]" (Aktualizacja 2026-06-20) | = F wprost (zaktualizowane #34) |
| S8 | `core/_meta_latex/status_map.tex:219` | σ_ab = „Hipoteza"; M_*=„Propozycja, brak mikro-derywacji"; c_GW=„Twierdzenie" | status pokorny; R4 |
| S9 | `papers_external/tgp_english_summary/sec07_tests.tex:316` | „tensor sector (σ_ab) has status ``Sketch''" | status pokorny, brak overclaim |
| S10 | `papers_external/tgp_english_summary/sec07_tests.tex:138–140` | „amplitude matched to GR via ξ_eff" | R3 — matching |
| S11 | `tgp_letter.tex:345` | mod oddechowy = 3. polaryzacja; c_GW exact | R2+R4 — niezależne od κ_E |
| S12 | `papers_external/tgp_english_summary/sec07_tests.tex:349` | „TGP (3 free parameters): Ω_m, H_0, Φ₀" | inny kontekst (fit kosmologiczny DE), nie param lagranżjanu |

## §B — DO-POPRAWY (twarde; konflikt z F lub rem:sigma-params)

| # | Plik / kotwica | Treść (cytat) | Dlaczego DO-POPRAWY | Proponowana poprawka (addytywna) |
|---|---|---|---|---|
| **P1** | `TGP_FOUNDATIONS.md:778–779` CL-7 | „Status LIVE sektora radiacyjnego: **UNDERDETERMINED** (… domknięcie = **pinowanie κ_E z substratu**)" | #33/#34 dowiodły, że pinowanie κ_E jest NIEMOŻLIWE (UV-czuły + brak bieguna). „Domknięcie = pinowanie" już rozstrzygnięte negatywnie | Dodać annotację **CL-8 (2026-06-22)**: κ_E = nieredukowalny wolny parametr (dowiedzione #33/#34); domknięcie = opcja (b), NIE pinowanie. Status LIVE → RESOLVED-via-free-parameter |
| **P2** | `core/sek07_predykcje/sek07_predykcje.tex:246` | „TGP ma **2 wolne parametry**: Φ₀ … i a_Γ" | Globalne (niezakresowane) „2" pomija nieredukowalny C_σ; sprzeczne z rem:sigma-params (3). Sama sekcja w wierszu GW (l.150) wymienia σ_0 | Dopisać: sektor tensorowy dodaje 1 nieredukowalny parametr (C_σ≡κ_E, rem:sigma-params/rem:sigma-Csigma-free) ⇒ łącznie 3; parsimony ratio cytowany dla sektora skalarno-flavorowego |
| **P3** | `core/_meta_latex/tabela_epistemiczna.tex:22` | „TGP posiada **2 wolne parametry**:" | Jak P2 — globalne „2" pomija C_σ | Dopisać cross-ref: + 1 nieredukowalny C_σ sektora tensorowego (rem:sigma-Csigma-free); 2 (skalar) / 3 (z tensorem) |
| **P4** | `README.md:97–102, 139` narracja OP-7 | „OP-7 CLOSED — … observationally complete … h_TGP_peak ≈ h_GR (6% fit)"; „σ_ab tensor (T2–T5 closed: … GW150914/GW170817 numerical fit)" | Amplituda GW (∝ κ_E=C_σσ₀²) jest dziś dowiedzenie wolnym parametrem; „fit" osiągalny przez **wybór** κ_E, nie predykcja bezparametrowa. „CLOSED/complete" przecenia | Dodać notę statusową (2026-06-20 #33/#34): normalizacja amplitudy tensorowej κ_E=C_σσ₀² = nieredukowalny wolny parametr; dopasowanie GW150914/170817 warunkowe na κ_E; treść modowa (2TT+breathing, c_GW=c) pozostaje ważna |

## §C — NIEJASNE (poprawne w wąskim zakresie; lekki cross-ref, niski priorytet)

| # | Plik / kotwica | Treść | Uwaga | Rekomendacja |
|---|---|---|---|---|
| N1 | `core/sek08_formalizm.tex:197` rem:parsimony, `eq:parsimony-N2` | „Warstwa II … dwa wolne parametry … N_param=2" | Zakresowane do Warstwy II (skalar) — broni się; ale „N_param=2" w równaniu brzmi globalnie | Dopisać cross-ref do rem:sigma-params (tensor +1 → 3 globalnie); ratio dotyczy sektora skalarnego |
| N2 | `papers_external/tgp_english_summary/sec06_formalism.tex:376` | „two free parameters (Φ₀, a_Γ) generate ≥12 predictions" | Jak N1 (zewnętrzny paper); wewn. niespójność z sec07_tests:349 (kontekstowo różne „3") | Opcjonalny przypis: sektor grawitacyjny tensorowy dodaje κ_E (osobny) |
| N3 | `core/sek07_predykcje.tex:151` wiersz GW | status „**Zamknięty** (tw. amplitude-matching)" | Amplituda matched (R3), nie predykcja; l.250 już rozróżnia matching/W, ale label „Zamknięty" myli | Dopisać krótki caveat: amplituda matched, normalizacja κ_E wolna (rem:sigma-Csigma-free) |
| N4 | `tgp_companion.tex:914, 1186`; `README.md:126` | „Free params 7 / ratio 5.7"; „40 predictions from 3 inputs" | Inne księgowanie (redukcja SM+ΛCDM); κ_E sektora grawitacyjnego poza tym zestawem | Opcjonalny przypis o κ_E jako parametrze sektora grawitacyjnego (nie wśród 7) |

## §D — Potwierdzenia ochronne (anti-Lakatos, §2 LOCK)
- ✅ GW4 (m_σ²/m_s²=2, LOCKED) — NIE downgradowany (R1). Występuje jako OPE/masa, odróżnione od C_σ.
- ✅ GW2/GW3/GW5/GW6 — NIE ruszane (R2). Audyt ich nie dotyka.
- ✅ thm:amplitude-matching — pozostaje warunkiem dopasowania (R3), bez przekwalifikowania.
- ✅ c_GW=c₀ (GW1) — niezależne od κ_E (R4).

## §E — Werdykt Phase 1
**DO-POPRAWY**: 4 trafienia twarde (P1–P4) + 4 niejasne/miękkie (N1–N4).
Propagacja statusu z #34 była **niepełna**: zaktualizowano sek08+dodatekQ+registry,
ale pominięto sek07_predykcje, tabela_epistemiczna, FOUNDATIONS CL-7 oraz README.
Wynik **nie jest negatywny** — istnieją realne niespójności do skorygowania.
Poprawki: addytywne/korekcyjne, zero nowych stałych; wymagają autoryzacji + build verification (P1–P4 w core; N1 core).
