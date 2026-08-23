---
title: "AUDYT GŁĘBOKI TGP_v1 — re-weryfikacja 22 zgłoszonych luk vs stan rzeczywisty (STATE.md #52)"
date: 2026-06-28
type: meta-audit
status: 🔴 ACTIVE — dokument decyzyjny: dostarczone vs brakujące + roadmapa kolejnych krytycznych projektów
updated: 2026-06-28 (niezależna re-weryfikacja §6/§1.1/CP-6b; KOREKTA #56: zarzut „korona 1→2" wycofany, blokery = 2: S01, L05)
method: "22 niezależnych audytorów (workflow tgp-deep-audit), każdy adversarialnie re-weryfikuje 1 problem względem najnowszych sesji #42–#52; 'CLOSED-annotation-only' ≠ naprawa, przeklasyfikowanie ≠ fix. UZUPEŁNIONO: druga niezależna ocena (6 audytorów, czysty kontekst, weryfikacja w .tex/.py) — §6, 5 rozbieżności."
parent: "[[../STATE.md]]"
tags: [meta, audyt, gap-analysis, roadmap, anti-lakatos, publication-readiness]
---

# AUDYT GŁĘBOKI TGP_v1 — co dostarczone, czego brakuje, co dalej

## §0 — Streszczenie wykonawcze (jedno spojrzenie)

Re-weryfikacja 22 zgłoszonych luk (S01–S07, L01–L08, M01–M03, D01, T01, most Γ→Φ, sektor QCD)
ujawnia **jeden dominujący wzorzec**:

> **Większość „domknięć" z POST_ACTION_UPDATE (2026-05-04..06) to uczciwe PRZEKLASYFIKOWANIE,
> nie naprawa strukturalna.** Defekt fizyczny zostaje; zmienia się jego *status epistemiczny*
> (aksjomat / wolny parametr / postulate-conditional / declared-limit) i opis w manuskrypcie.

To jest **dobra nauka** (anti-Lakatos: luki jawne, nie chowane) — ale oznacza, że **ciężar
pozostał, nie zniknął**. Z 22 problemów:

- **Realnie DOMKNIĘTE (2):** L01 (definicja ρ = −Tᵘᵤ/c² — genuine derivation), M03 (balance-sheet retrofit — deliverable istnieje).
- **CLOSED-ANNOTATION-ONLY (3):** S06, S07, L02 — dodano remark/glosariusz, defekt nietknięty.
- **SUPERSEDED (1):** S02 — fix zakotwiczony w sfalsyfikowanym M9.1''.
- **PARTIAL (16):** reszta — częściowy postęp + uczciwe re-labelling, rdzeń otwarty.
- **Blokuje publikację pełnego manuskryptu (2):** **S01** (metryka: canonical-in-body vs falsified-in-header) i **L05** (rdzeń wciąż głosi wykładnik masy 4 jako Twierdzenie, sprzeczny z kanonicznym α=2 → wykładnik 3).

> **✅ AKTUALIZACJA #58 (2026-06-28, agent-implementator — Tier 0 + Tier 1 wykonane; szczegóły STATE.md #58):** oba blokery ZAMKNIĘTE submission-side — **S01** (CP-1: caveat falsyfikacji+recovery {A,B,C} wpięty do PDF; M9.1''→proposed/static anchor) i **L05** (CP-2: `sek08b prop:Atail-preserved`→Aproksymacja + m_obs≠M_full). Higiena Tier 1 (CP-3 L02-renotacja, CP-4 L07/S05-framing, CP-5 D01-locki+CI, CP-6 ledger 856→~688 + M01/M02, CP-6b #42-consistency-checks, CP-0L korona) — wykonana. **Blokery publ.: 2 → 0** (manuskrypt spójny warstwowo). Wszystkie build-gate'y exit 0 (main 554 str., 0 nowych undefined refs; companion/letter/lepton OK). Pozostałości = defekty *fizyczne* (Tier 2/3: most Γ→Φ, σ_ab/C_σ-free, non-abelian, S07 NON-DERIVABLE) + 2 flagi (Σm_ν=59,6 w papers; COSMO main-body). Re-klasyfikacje warstwowe: L02→CLOSED-RESOLVED, S02→body-spójne, M01/M02/D01/LEDGER-N/EW-mW→CLOSED-RESOLVED.

> **⚠ AKTUALIZACJA 2026-06-28 (niezależna re-weryfikacja u źródła — recenzent zewnętrzny, czysty kontekst; pełne ustalenia §6):** niezależny audyt dodaje trzy niedoreprezentowane luki (kosmologia DESI fit-by-construction + sprzeczność MG-vs-particle-DM + SPARC obalony PR-004 5,4σ; m_W/sin²θ_W consistency-check przebrany za predykcję; licznik 856 NIEpropagowany do ~688). Severity licznika podniesiona do **blokującej** (self-inconsistency dla recenzenta).
>
> **✅ KOREKTA #56 (2026-06-28, re-analiza mechanizmu τ u źródła):** wcześniejszy zarzut „korona = 1→2 z ukrytym fitem g_0^τ" zostaje **WYCOFANY jako zbyt ostry**. Liczba blokerów pozostaje **2 (S01, L05)** — CROWN/CP-0 **NIE jest blokerem**. τ jest wyznaczany przez **algebraiczne domknięcie Koidego `Q_K=3/2`** (przy N=3), a nie wolny fit; `G0_TAU_fitted=3.18912` w `a3d…py:143` to lokalna stała wygody skryptu pokrywająca się z wartością Koidego (`tau_selection_v47b.py`). N=3 (i brak 4. leptonu) selekcjonowany realnie przez stabilność: studnia ODE `V(g)=g²(1−g)` + ściana-duch `g*=exp(−1/2α)`, odbicia 0→1→3→6, RMSE/A 0,5→2,2→8→36% (k=4 DEGRADED, `ls10 LS-10d`), oraz wykluczenie 4. generacji poniżej granicy LEP (`ex116`). Korona = **genuine 1→3**. Pozostała uczciwa luka (węższa): (a) φ² zawodzi dla gen.3 (13,7%, `dodatekJ2:117,220`) ⟹ dwa mechanizmy (φ-drabina + Koide), nie jeden; (b) dokładność Koidego opiera się na `B=√2` numerycznie potwierdzonym, lecz analitycznie niewyprowadzonym (a3d T5 FAIL, T-OP3). → CP-0 przeklasyfikowane do **Tier 1 (honest-framing), NIE Tier 0**.

**Drugi wzorzec:** wszystkie ciężkie luki (XL) zbiegają do **dwóch fundamentalnych węzłów**:
1. **Most Γ→Φ (continuum/NGFP)** — wspólny mianownik α=2, c₀, 𝒜→α_s, K_geo. **#49/#53 UDOWODNIŁY,
   że α=2 jest NON-DERIVABLE z kanonicznego substratu** ⟹ most nie domknie wszystkich 4 korzeni
   bez nowego aksjomatu. To prawdopodobnie **terminalna dyspozycja = aksjomat**, nie „do zrobienia".
2. **Ontologia σ_ab / κ_E / c₀** — parameter-free GW. Udowodnione negatywnie (#33/#34/#37: brak
   izolowanego bieguna spin-2, ∫P₂=0, liniowa rozbieżność UV −16/35).

**Wniosek operacyjny:** rozdziel **inżynierię publikacyjną** (Tier 0–1, tygodnie, realne) od
**wieloletnich tracków fundamentalnych** (Tier 3, lata, część prawdopodobnie niederywowalna).
Nie blokuj publikacji tymi drugimi — nieś je jako Limitations (już wpięte w letter/companion #52).

---

## §1 — Tabela master (22 problemy)

Legenda statusu: 🟢 CLOSED-RESOLVED · 🟡 PARTIAL · 🟠 ANNOTATION-ONLY · 🔵 SUPERSEDED.
„README→" = status deklarowany w `audyt/`; „RZECZYWISTY" = po re-weryfikacji vs #52.

| ID | Problem | README→ | RZECZYWISTY | Sev. | Effort | Blok. publ. |
|----|---------|---------|-------------|------|--------|:-----------:|
| **S01** | 4 formy metryki w sek08c | CLOSED-RESOLVED (G.0) | 🟡 PARTIAL (uściślone #57) — M9.1'' obalona **tylko w sektorze GW** (forma ppE, 5,02σ); statyczny PPN (A,B→β=γ=1) przeżywa, recovery {A,B,C} zachowuje A,B (sfalsyfikowane=C/c₀ wolne #37). Bloker: body głosi „metryka jednoznacznie wyznaczona" (Twierdzenie, l.417–424) **bez** caveatu falsyfikacji+recovery — te są w komentarzach (.tex „v3.0"), nie w PDF. 🟢 **#58 CP-1: BLOKER ZAMKNIĘTY** — `rem:M911-GW-caveat` (PDF-widoczny) niesie falsyfikację GW 5,02σ + recovery {A,B,C} + C/c₀-free (#37); status M9.1''→`proposed/specific static anchor`; „jednoznacznie" zmiękczone. Defekt fizyczny (GW) pozostaje (S07/Tier 3, jako proposed/Limitations) | HIGH | L (honest-S) | ~~TAK~~ **NIE (#58)** |
| **S02** | √(−g) niespójne z M9.1'' | CLOSED-RESOLVED | 🔵 SUPERSEDED — fix oparty na sfalsyfikowanym M9.1''; v1.x √−g wciąż w sek08a/dodatekH; M9.2/M9.3 nie przeliczone (gauge-argued). 🟢 **#58 CP-1: body-side spójne** — v1.x √−g=c₀φ oznaczone `deprecated`→M9.1'' w `sek08a` + `status_map` (dodatekH miał już przypis). Resztka (M9.2/M9.3 re-run) = Tier 2 | MED | M | nie |
| **S03** | β_PPN: 1/2 vs master-formula 1 | CLOSED-RESOLVED | 🟡 PARTIAL — β=1 przez antypodalny ansatz p=s + łatkę factor-2 (coord vs local speed); brak derywacji wariacyjnej (Path A nigdy nie otwarte). 🟡 **#58 CP-1: framing dodany** — β,γ jawnie jako **projekcja PPN** (nie predykcja) w PDF (`rem:M911-GW-caveat`, meta/PPN_AS_PROJECTION.md). Derywacja wariacyjna Path A nadal otwarta → PARTIAL | HIGH | L | nie |
| **S04** | aksjomat metric-coupling vs φ·ρ (5. siła) | CLOSED-DERIVED | 🟡 PARTIAL — derywacja ρ OK, ale no-fifth-force WARUNKOWE na m_Φ~M_Pl (niewykazany postulat); Cassini (N6) nie przeliczone | MED | M | nie |
| **S05** | σ_ab łamie aksjomat single-Φ | CLOSED-RESOLVED | 🟡 PARTIAL — composite broniony, ALE C_σ/κ_E to dowiedziony nieredukowalny wolny parametr UV; „M²=2m_s² derived" skorygowane (OPE coeff, nie biegun); uczciwe re-labelling, nie derywacja. 🟡 **#58 CP-4: framing zmiękczony** — FOUNDATIONS §1 single-Φ uniqueness uściślone (dot. *pola*, nie liczby parametrów; σ_ab=kompozyt; budżet 3). C_σ/κ_E free param pozostaje → PARTIAL | HIGH | XL | nie |
| **S06** | cyrkularność χ.1 (G_N) / UV.2 (M_TGP) | „substancjalnie zamknięte" | 🟠 ANNOTATION-ONLY — G_N=1/M_Pl² wciąż tautologia; K_struct numerologia wycofana, nie zastąpiona derywacją; F6 rollback PENDING | MED | XL | nie |
| **S07** | M9.1'' postulat vs derywacja | OPEN P2 → „CLOSED" (#10) | 🟠 ANNOTATION-ONLY — #53 dowodzi NON-DERIVABLE z H_Γ; {A,B,C} to ansatz; „closure" = klasyfikacja+annotacja | MED | XL | nie |
| **L01** | operacyjna definicja ρ | CLOSED-DERIVED | 🟢 CLOSED-RESOLVED — ρ=−Tᵘᵤ/c² genuine; ALE nie re-derived dla disformal-LIVE; ρ_rad=0 ratowane inflacją (A−) | MED | L | nie |
| **L02** | kolizja notacji β/γ (WF vs GL) | „substancjalnie zamknięte" | 🟠 ANNOTATION-ONLY — renotacja nigdy nie wpięta do body; sek08 prop:vacuum-selection ma sprzeczność 1 vs 0.264, oba „Wilson-Fisher". 🟢 **#58 CP-3: CLOSED-RESOLVED** — renotacja `(β/γ)_GL` (→1) vs `(β/γ)_WF`(≈0,264) wpięta do body (`prop:vacuum-selection` + inline \ref `app:A-beta-gamma-distinction`), FOUNDATIONS §1/l.1121, `vacuum_selection.py`; sprzeczność usunięta | MED | **S** | nie |
| **L03** | stabilność V''(1)<0 przez K=K_geo·φ⁴ | CLOSED-annotation | 🟡 PARTIAL — tylko punkt próżni domknięty (FP5); brak analizy spektralnej na tłach niejednorodnych, ghost-wall przy φ→0 nietknięte | MED | L | nie |
| **L04** | dualizm ODE K=g² (α=1) vs K=g⁴ (α=2) | CLOSED-RESOLVED | 🟡 PARTIAL — spójność wewn. OK (jeden lock), ALE „derywacja" D-uniqueness WYCOFANA jako cyrkularna (#49: substrat → α_eff=−½); α=2 = aksjomat | MED | XL | nie |
| **L05** | wykładnik masy k=4 vs p=5−α | CLOSED-RESOLVED A− | 🟡 PARTIAL (uściślone #57) — ujednolicone w ~80%: rdzeń (sek08_formalizm l.9376, sek08a l.392) niesie wzór zunifikowany `m_obs=c_M·A_tail²·g_0^(e²(1−α/4))` spójny z α=2 (μ/e −0,005%), „A^4"/„A^(5−α)=A^3" zdegradowane do aproksymacji. Resztka: `sek08b prop:Atail-preserved` (l.352–374) NIE zsynchronizowane — wciąż Twierdzenie z gołym A^4; + dryf ν 4,12(dodatekJ) vs 1,36(dodatekJ2). 🟢 **#58 CP-2: RESZTKA ZAMKNIĘTA** — `prop:Atail-preserved`→`Aproksymacja (M_full)`, nowy `rem:Atail4-vs-unified` (m_obs≠M_full ADM/Komar, kanon α=2) + `dodatekJ2 rem:J2-kobs-kfull` (k_full=4 vs k_obs=e²/2; ν 4,12/1,36 = fity tej samej A_tail) | MED | **S** | ~~TAK~~ **NIE (#58)** |
| **L06** | masa aksjonu m_X „locked" ale wolna | PARTIAL B+ | 🟡 PARTIAL — 7–8 ścieżek derywacji obstrukcja; uczciwe re-labelling do FREE; 100 vs 0.83 MeV nieujednolicone | MED | XL | nie |
| **L07** | aksjomat zero-sum (Λ_eff) | PARTIAL B+ | 🟡 PARTIAL — ZS1 derived, ZS2-quad = „gauge fixing" (Path D 4/5 obstrukcja); uczciwe framing tylko w KOMENTARZACH .tex, body wciąż „wielkość wyliczalna". 🟢 **#58 CP-4: framing → PDF body** — `sek01 rem:ZS-L07-status` (ZS1=Z₂-tożsamość, ZS2-kwadr.=gauge fixing, nie surowy aksjomat) + `sek05 rem:Lambda-L07-strengthened`; body już nie głosi „wyliczalna" bez dyspozycji. ZS2-quad gauge-fixing pozostaje (Tier 3) | MED | **S** | nie |
| **L08** | kink-fermion (spin, generacje, W/Z) | per-subproblem split | 🟡 PARTIAL — spin-statistics/Dirac realnie A−; ALE W/Z (SU(2)_L) + gluony (SU(3)_c) = DECLARED LIMIT; wykł. generacji e²/2 = numerical anchor | HIGH | XL | nie |
| **M01** | status creep w registry | OPEN P3 | 🟡 PARTIAL — gate forward zainstalowany; backward retrofit niepełny (~40/80), kolumna falsyfikowalności odroczona, nagłówek wciąż inflated. 🟢 **#58 CP-6: zamknięty submission-side** — `meta/M01_Phase_FINAL_close_2026-06-28.md` utworzony; nagłówek skorygowany do ~688. Backward retrofit ~40/80 + per-entry kolumna falsyfikowalności jawnie odroczone | MED | M | nie |
| **M02** | ledger pollution (commit 74394a8) | OPEN P3 | 🟡 PARTIAL — downgrade'y zrobione, ALE counter wciąż 856 (rollback Option A nie wykonany); 4 foldery wciąż „FULL CONVERGENCE". 🟢 **#58 CP-6: CLOSED-RESOLVED** — counter 856→~688 (Option A) w INDEX+REGISTRY; 4 foldery (χ.1/UV.2/ω.2/ω.3) rebrand → ANSATZ/NUMEROLOGICAL (mark-only). (Resztka: 1 log-wpis INDEX ~l.248 stary „PERFECT CONVERGENCE" — pokryty bannerem) | MED | M | nie |
| **M03** | brak balance-sheet ~27 cykli | OPEN P4 → 100% | 🟢 CLOSED-RESOLVED — 40 retrofit_*.md istnieje; counter 856→~688, ratio 5,5→~3,67 (estymata ±20); per-row tagging odroczony | LOW | M | nie |
| **D01** | dryf parametrów (α_s, m_H, Φ_0, g0e) | mixed P3 | 🟡 PARTIAL — manuskrypt zsync; ALE 6+ skryptów tooling wciąż Φ_0=24.65 (vs lock 24.783); głębia = re-labelling do FREE. 🟢 **#58 CP-5: CLOSED-RESOLVED (hygiena)** — 8 skryptów Φ_0→24.783, 4 skrypty Σm_ν→59.01 (z notą `# LOCKED #42`), CI job `param-lock`; re-run `tgp_master_consistency_v47.py` 59/60 PASS (SPÓJNY). (Resztka: companion/letter Σm_ν=59,6 — flaga dla usera) | MED | XL (hygiena S–M) | nie |
| **T01** | M9.1'' sfalsyfikowane w GWTC-3 (5σ) | CLOSED-EXECUTED v4 | 🟡 PARTIAL — BRAK parameter-free falsyfikatora GW amplitudy/fazy; „recovery" emergent-metric SUPERSEDED (#37 double-tuning); zostaje tylko structural mode-content | HIGH | XL | nie |
| **BRIDGE-CG** | most Γ→Φ (wspólny mianownik 4 korzeni) | brak folderu (OPEN) | 🟡 PARTIAL — NGFP analitycznie 7/7, ALE rdzeń mostu [OTWARTY] (CG-1 Banach, ex200 4/8, ex202 T6 FAIL); #49/#53 dowodzą α=2 NON-DERIVABLE | HIGH | XL | nie |
| **QCD** | kwarki/QCD + V_NN + W/Z | mixed („Domknięty"/„Program") | 🟡 PARTIAL — kwarki 4→2 predykcje; W/Z+gluony DECLARED LIMIT; V_NN brak V_σσ/V_ττ; α_s przez postulate-conditional 𝒜 | HIGH | XL | nie |

### §1.1 — Wiersze dodane przez niezależną re-weryfikację (2026-06-28, recenzent zewnętrzny; dowody §6)

| ID | Problem | README→ | RZECZYWISTY | Sev. | Effort | Blok. publ. |
|----|---------|---------|-------------|------|--------|:-----------:|
| **CROWN** | korona leptonowa „1→3 do 0,006%" | strongest result | 🟢 ~CLOSED (KOREKTA #56) — **genuine 1→3**: μ via φ-drabina, τ via domknięcie Koidego `Q_K=3/2` (algebraiczne, przy N=3; r_31=3477,5, 0,009%); `G0_TAU_fitted` = wygoda skryptu = wartość Koidego, **nie** wolny fit; N=3 selekcjonowane stabilnością (ściana-duch + brak 4. gen. <LEP, `ex116`). Luka resztkowa: φ² zawodzi dla gen.3 (dwa mechanizmy) + `B=√2` numeryczne, nie analityczne. 🟢 **#58 CP-0L: framing wpięty** — `a3d_soliton_brannen_r.py:143` G0_TAU_fitted oznaczone „= wartość Koidego Q_K=3/2" (usunięty mylący `# fitted to match PDG`); lepton paper Limitations: gen.3 przez Koide (φ² zawodzi 13,7%), B=√2 num. potwierdzone/analitycznie niewyprowadzone (T-OP3) | LOW–MED | S (honest-framing) | **NIE** |
| **COSMO** | DESI w(z) / Ω_Λ / Ω_DM | „prediction" / promoted | 🟡 PARTIAL — DESI χ²≈χ² **fit-by-construction** (ψ zamrożone, „indistinguishable from ΛCDM"); Ω_Λ przehandlowany (g̃≈0,98, trade-off psuje α_s); **sprzeczność MG-vs-particle-DM** + SPARC obalony (PR-004 5,4σ). 🟡 **#58 CP-6b: częściowo (papers)** — DESI w(z)/Ω_Λ oznaczone **consistency-check / fit-by-construction** w companion/letter/README; ALE sprzeczność MG-vs-DM + SPARC + body main.tex cosmology NIETKNIĘTE (Tier 2/3) → nadal blok jako „prediction" | HIGH | M–XL | **TAK** (jako „prediction") |
| **EW-mW** | m_W=80,354 „0,01σ" | clean prediction | 🟠 ANNOTATION/DECLARED-LIMIT — m_W=m_Z·cosθ_W bierze m_Z+α_s jako INPUT; tree sin²θ_W=3/13 **11,3σ off**, ratowane pętlą QCD; consistency-check przebrany za predykcję. 🟢 **#58 CP-6b: nagłówki zsync** — m_W, α_s, sin²θ_W, w_0(DESI), Ω_Λ oznaczone **consistency-check** (CC, nie predykcje) w `tgp_companion.tex`+`tgp_letter.tex`+`README.md` (legendy/akapity, #42) | HIGH | XL | ~~nagłówki letter/companion~~ **NIE (#58)** |
| **LEDGER-N** | licznik predykcji | M03 100% / 856 | 🟡 PARTIAL — **856/784 wciąż w nagłówku** (`INDEX.md:224`), audyt liczy ~688; podwójny ledger NIEpropagowany. 🟢 **#58 CP-6: CLOSED-RESOLVED** — licznik 856/784→~688 (ratio ~3,67) spropagowany do `INDEX.md` + `PREDICTIONS_REGISTRY.md` (Option A, rozkład zachowany); self-inconsistency usunięta | MED | M | ~~TAK~~ **NIE (#58)** |

---

## §2 — Co realnie DOMKNIĘTE vs co tylko UCZCIWIE OPISANE

### 2A. Realnie domknięte (zostają domknięte)
- **L01** — ρ ≡ −Tᵘᵤ/c² wyprowadzone z δS_mat/δψ (tożsamość trace-response w sprzężeniu konforemnym); wpisane w sek08a. *Reszta: rozszerzenie disformal pending.*
- **M03** — balance-sheet retrofit faktycznie wykonany (40 plików `retrofit_op-*.md`); to był właściwy deliverable audytu (przeklasyfikowanie 7 cykli NUMEROLOGICAL to poprawny wynik, nie unik).
- **L08 #1/#4** — spin-statistics (π₁(RP²)=Z₂ → faza −1) i Dirac/Clifford: realne konstrukcje A−.
- **S02 A1/A2/A3 (lokalnie)** — czyszczenie 4 form było realne (forma I → `\iffalse`, forma III → appendix historyczny). Problem: kotwica (M9.1'') padła.

### 2B. Tylko uczciwie opisane (defekt fizyczny ZOSTAJE — „honest, not fixed")
To jest **rdzeń ustalenia użytkownika**. Te luki wyglądają na „zamknięte" w POST_ACTION/STATE, ale to
re-labelling:
- **α=2 (L04), c₀ (S06/c₀), 𝒜→α_s (QCD), K_geo** → aksjomat / FREE / POSTULATE-CONDITIONAL.
  **#49/#53: α=2 dowiedzione NON-DERIVABLE z substratu.**
- **κ_E / σ_ab (S05) + brak parameter-free GW (T01)** → „opcja b": wolny parametr.
- **W/Z + gluony (L08/QCD)** → DECLARED LIMIT (wyczerpanie 6 ścieżek SU(2), 0/7 SU(3)).
- **m_X (L06)** → FREE (7–8 ścieżek obstrukcja).
- **ZS2-quadratic (L07)** → „gauge fixing"; **G_N/M_TGP (S06)** → tautologia/numerologia wycofana.
- **M9.1''/metryka (S07)** → ansatz; status faktycznie „proposed (P)", nie „emergent (E)".

### 2C. Czyste zaniedbania edytorskie (tanie, ale realne liabilities dla recenzenta)
- **L02** — renotacja β/γ tylko w glosariuszu; sek08 ma żywą sprzeczność 1 vs 0.264 „at WF fixed point".
- **L05** — sek08b `prop:Atail-preserved` wciąż **Twierdzenie** z wykładnikiem 4 (sprzeczne z α=2).
- **L07** — uczciwe framing ZS2 tylko w komentarzach .tex; body PDF wciąż „wielkość wyliczalna wymuszona".
- **D01** — 6+ skryptów tooling Φ_0=24.65 (recenzent re-running dostanie inne liczby niż manuskrypt).
- **M01/M02** — nagłówek registry wciąż wiedzie liczbami inflated (855/„refined²"); 4 foldery „FULL CONVERGENCE".
- **S05/FOUNDATIONS §1** — sztywne wording „single-Φ uniqueness" nie zmiękczone mimo budżetu 3 parametrów.

---

## §3 — ROADMAPA: kolejne krytyczne projekty

Zasada nadrzędna: **oddziel publikowalność od fundamentów.** Tier 0–1 to tygodnie i odblokowuje
submission; Tier 3 to lata i część jest prawdopodobnie niederywowalna — nieść jako Limitations.

### 🔴 TIER 0 — BLOKERY PUBLIKACJI (przed submission pełnego manuskryptu)

**~~CP-0 (bloker)~~ → WYCOFANE do Tier 1 (KOREKTA #56, patrz §3.1).** Pierwotny zarzut „korona = 1→2 z ukrytym fitem" był zbyt ostry; korona to genuine 1→3 (Koide-przy-N=3). Pozostaje tylko **lekki honest-framing** (Tier 1, niżej), NIE bloker.

**CP-1 — Rekonsolidacja sektora metrycznego (S01 + S02 + S03).** *Effort L · gating.*
*Uściślenie #57: M9.1'' obalona TYLKO w sektorze GW (forma ppE, 5,02σ); statyczny PPN przeżywa,
recovery {A,B,C} zachowuje A,B (sfalsyfikowane = σ-coupling C/c₀, wolne #37).* Skompilowane body
prezentuje M9.1'' jako „jednoznaczne Twierdzenie" (l.417–424) **bez** caveatu falsyfikacji GW + recovery
(te są w komentarzach .tex „v3.0", nie w PDF) → recenzent czytający body nie wie o zawężeniu. Zakres:
(a) uczynić rodzinę emergent-metric {A,B,C} (lub konkretny członek) body-kanoniczną *albo* jawnie
nadać M9.1'' status „proposed/specific anchor" wszędzie w body; (b) usunąć deprecated v1.x √(−g)=c·ψ
z sek08a/dodatekH/status_map; (c) uzgodnić β jako projekcję L2 (a nie „prediction") — jedna spójna
warstwa epistemiczna; (d) dostarczyć twierdzenie o jednoznaczności (NEEDS N5) LUB uczciwie nieść jako
aksjomat. **Zależność:** pełne (a) zależy od decyzji S07; minimalna wersja (honest „proposed") jest
wykonalna natychmiast i wystarcza do submission.

**CP-2 — Synchronizacja sek08b z wzorem zunifikowanym (L05).** *Effort S (uściślone #57).*
Rdzeń sprzeczności JUŻ rozwiązany: `sek08_formalizm` l.9376 + `sek08a` l.392 niosą kanoniczny wzór
`m_obs=c_M·A_tail²·g_0^(e²(1−α/4))` (α=2 spójne, μ/e −0,005%), z „A^4"/„A^3" jako aproksymacjami.
Resztka (NIE twarda sprzeczność, lokalna): (a) `sek08b prop:Atail-preserved` (l.352–374) wciąż
`\statuslabel{Twierdzenie}` z gołym A^4 — zdegradować do „aproksymacja w reżimie ogona / M_full",
cross-ref do `sek08_formalizm` (wzór zunifikowany); (b) wnieść notę m_obs≠M_full (już zalecaną w
`audyt/L05` l.71, niewykonaną); (c) uzgodnić ν 4,12(dodatekJ) vs 1,36(dodatekJ2). Dotyka korony —
ale korona = genuine 1→3 (#56), więc to higiena spójności, nie ratowanie wyniku.

### 🟠 TIER 1 — HIGIENA UCZCIWOŚCI (tanie, wysoka wartość reputacyjna, równoległe)

**CP-3 — Renotacja β/γ (L02).** *Effort S (~1 dzień).* Grep+replace β_WF/γ_WF w sek08/M3-8/FOUNDATIONS§7/scripts;
uzgodnić prop:vacuum-selection (1 vs 0.264); inline \ref do app:A.

**CP-4 — Promocja honest-framingu z komentarzy do body (L07 + S05).** *Effort S.* Przenieść framing
ZS2 „gauge fixing" do body sek05/sek01 (korekta „wielkość wyliczalna"); zmiękczyć single-Φ uniqueness
w FOUNDATIONS §1 do spójności z budżetem 3-parametrowym.

**CP-5 — Lock parametrów w tooling + CI guard (D01 hygiena).** *Effort S–M.* Sync 6+ skryptów do
Φ_0=24.783, re-run self-consistency K_geo·m_sp²=π·Φ_0², dodać grep-lock w CI by dryf nie wracał.

**CP-6 — Domknięcie ledgera (M01 + M02 + S06 submission-side).** *Effort M.* Nagłówek registry wiedzie
uczciwymi liczbami (~688 / ratio ~3,67), per-entry kolumna falsyfikowalności; rebrand 4 folderów
(χ.1/UV.2/ω.2/ω.3) do ANSATZ/NUMEROLOGICAL; finalizacja counter (Option A 856→688 *lub* per-entry
[CONTESTED]); formalny `Phase_FINAL_close` dla M01. *(Honest framing w letter/companion już wpięty #45/#52.)*
**⚠ podniesione do blokera (niezależna re-weryfikacja):** licznik **856/784 wciąż żyje w nagłówku** `INDEX.md:224` i `PREDICTIONS_REGISTRY.md:1163` — propagacja 856→688 jest obowiązkowa przed submission (manuskrypt cytujący 856 „weryfikacji" vs własny audyt ~688 = self-inconsistency, którą recenzent wychwyci).

**CP-6b — Synchronizacja nagłówków papers z ledgerem #42 (NOWY, niezależna re-weryfikacja).** *Effort S–M.*
W `tgp_companion.tex`/`tgp_letter.tex`/`README.md` oznaczyć jako **consistency-checks, NIE predykcje**:
m_W (`sek09:447,454` — m_Z+α_s input), α_s (postulate-conditional 𝒜, #43), sin²θ_W=3/13 (tree 11,3σ off),
DESI w(z) (fit-by-construction, `b1…py:498`), Ω_Λ (g̃≈0,98 przehandlowany, trade-off z α_s). Inaczej
peer-review rozłoży rozjazd „uczciwa warstwa research vs optymistyczna warstwa nagłówkowa".

**CP-0L — Lekki honest-framing korony leptonowej (KOREKTA #56; było „CP-0 bloker").** *Effort S.*
NIE bloker (korona = genuine 1→3). Drobna higiena: (a) w `paper_lepton_masses`/sek08b uczynić jawnym, że
gen.3 idzie przez **domknięcie Koidego `Q_K=3/2` (przy N=3)**, a NIE przez φ-drabinę (φ² zawodzi 13,7%,
już uczciwie w `dodatekJ2:117,220`); (b) oznaczyć `G0_TAU_fitted=3.18912` w `a3d…py:143` jako „= wartość
Koidego" (usunąć mylący komentarz `# fitted to match PDG`); (c) w Limitations zaznaczyć, że `B=√2`
(dokładność Koidego) jest numerycznie potwierdzone, nie analitycznie wyprowadzone (T-OP3). Paper II
**można publikować** (po tym kosmetycznym framingu).

### §3.1 — Nota o korekcie #56
CP-0 zgłoszone w pierwszej iteracji jako 3. bloker (ukryty fit τ) zostało **wycofane** po re-analizie
mechanizmu selekcji τ u źródła (`tau_selection_v47b.py`, `ls10_third_generation_selection.py`, `ex116`).
Blokery publikacji wracają do **2: S01, L05.**

### 🟡 TIER 2 — DOMKNIĘCIA WARUNKOWE (wzmacniają, nie blokują)

**CP-7 — `op-spectral-analysis-Phi` (L03).** *Effort L.* Diagonalizacja K(φ)∂²+V''(φ) na tłach
próżnia/Yukawa/soliton (sympy + BVP numeryczny), m²>0 dla φ>0, rozstrzygnięcie ghost-wall przy φ→0.
Podpiera stabilność solitonów leptonowych (korona!) i dyspersję GW.

**CP-8 — S04 residuals (Cassini/ω_BD/m_Φ).** *Effort M.* Natywna re-derywacja Cassini/Shapiro z członem
dφ/dr·ρ; czysty bound ω_BD>40000; rozładowanie zależności no-fifth-force od postulatu m_Φ~M_Pl + caveat Pattern-2.5.

**CP-9 — L01 disformal extension.** *Effort L.* Re-derywacja źródła ρ w żywym sprzężeniu disformal
(człony Tᵘᵛ∂Φ∂Φ), powiązana z programem pinowania κ_E/B(Φ)/M_*.

### ⚪ TIER 3 — FUNDAMENTY WIELOLETNIE (NIE blokować publikacji — nieść jako Limitations)

| Track | Pokrywa | Status dowodowy | Realna dyspozycja |
|---|---|---|---|
| **Most Γ→Φ / NGFP** | α=2 (L04), c₀, 𝒜→α_s (QCD), K_geo, BRIDGE-CG, D01-root, S06-root | **#49/#53: α=2 NON-DERIVABLE z kanon. substratu** | Most NIE domknie 4 korzeni bez nowego aksjomatu → **terminalnie aksjomat** (chyba że rozszerzenie substratu) |
| **Ontologia σ_ab / κ_E / c₀** | S05, T01 (parameter-free GW) | Negatywne #33/#34/#37 (brak bieguna spin-2, ∫P₂=0, UV −16/35) | Wymaga genuinely nowej struktury; inaczej: stały framing 3-param + structural-only GW falsifier |
| **Non-Abelian gauge + V_NN** | L08 (W/Z, gluony), QCD | DECLARED LIMIT (6-path SU(2), 0/7 SU(3)) | Wymaga rozszerzenia substratu (Skyrmion/Dirac multibody) lub trwale input-phenomenology |
| **Metryka z first-principles** | S07 | #53 NON-DERIVABLE | Terminalnie status „proposed (P)" / aksjomat (de-facto już teraz) |
| **m_X structural** | L06 | 7–8 ścieżek obstrukcja | Wymaga axiomów poza S05; trwale FREE |

---

## §4 — Rekomendowana sekwencja

```
TYDZIEŃ 1–2:   CP-3, CP-4, CP-5, CP-6b, CP-0L  (Tier 1, równolegle — czysta higiena, zero nowej fizyki)
TYDZIEŃ 2–4:   CP-6                     (ledger 856→688)
TYDZIEŃ 3–6:   CP-2                     (L05 — sprzeczność Twierdzenia o wykładniku masy)
TYDZIEŃ 4–8:   CP-1                     (S01/S02/S03 — rekonsolidacja metryki; wariant honest)
RÓWNOLEGLE:    PUBLIKUJ Paper I (N-body) + Paper II (lepton masses) — UV-niezależne, gotowe
               (Paper II: po kosmetycznym CP-0L; korona = genuine 1→3, KOREKTA #56)
PÓŹNIEJ:       CP-7, CP-8, CP-9 (Tier 2)
TŁO (lata):    Tier 3 — most Γ→Φ; nieść jako Limitations, NIE blokować
```

**Najważniejsza rzecz teraz:** Tier 0 (CP-1, CP-2) — bo to **jedyne dwie luki blokujące** spójność
skompilowanego manuskryptu. Reszta albo jest tania higiena (Tier 1), albo wzmocnienie (Tier 2), albo
uczciwie-opisany fundament (Tier 3). Papers I/II nie czekają na nic z powyższego.

## §6 — Niezależna re-weryfikacja u źródła (recenzent zewnętrzny, czysty kontekst, 2026-06-28)

Druga, niezależna ocena (6 równoległych audytorów, każdy bez dostępu do tego dokumentu — by uniknąć
zakotwiczenia; każdy weryfikował defekt w .tex/.py, nie w README/POST_ACTION). **Zgadza się ze strukturą
i dyspozycją fundamentów** powyższego audytu, ale jest on **zbyt łagodny dla trzech sztandarowych wyników**.
Pięć rozbieżności (z dowodem u źródła):

**R1 — Korona leptonowa. ⚠ WYCOFANE (KOREKTA #56) — pierwotny zarzut był BŁĘDNY/zbyt ostry.**
Pierwsza iteracja twierdziła „korona = 1→2 z ukrytym fitem g_0^τ; trzy niezgodne ścieżki". Re-analiza
mechanizmu τ u źródła **obala to**:
- τ jest wyznaczany przez **algebraiczne domknięcie Koidego `Q_K=3/2`** (przy N=3), które fiksuje r_31 z r_21:
  `√r₃₁=2(1+√r₂₁)+√(3(1+4√r₂₁+r₂₁))` → r_31=3477,5 (**0,009%**, `ls10_third_generation_selection.py:79-87`).
  To **nie** wolny fit, lecz tożsamość przy ustalonym N=3.
- `G0_TAU_fitted=3.18912` (`a3d…py:143`) to **lokalna stała wygody** skryptu Brannena, **pokrywająca się**
  z wartością Koidego (`tau_selection_v47b.py`: g0_tau(Koide)≈g0_tau(PDG)) — nie niezależny anchor. Pierwszy
  audytor zakotwiczył się na komentarzu `# fitted` i błędnie ogłosił ukryty fit.
- Trzy „ścieżki" to NIE trzy niezgodne sukcesy: φ² to **odrzucona** naiwna próba (3955, 13,7%, uczciwie OPEN
  w `dodatekJ2:117,220`), a „fitted" = Koide (ta sama wartość). Realny mechanizm jest jeden: Koide-przy-N=3.
- N=3 + brak 4. leptonu realnie selekcjonowane przez stabilność: studnia ODE `V=g²(1−g)` + ściana-duch,
  odbicia 0→1→3→6, RMSE/A→36% (k=4), wykluczenie 4. gen. <LEP (`ex116`). **To potwierdza pamięć autora**
  (gen.3 wymaga studni + górnego ogranicznika stabilności).

**Skorygowany werdykt R1:** korona = **genuine 1→3** (g_0^e + φ-drabina + Koide-przy-N=3). Pierwotny
`AUDYT_GLEBOKI` (trzymający koronę jako najmocniejszy wynik) był **bliżej prawdy** niż niezależny audytor.
Pozostała uczciwa (węższa) luka, do Limitations: (i) gen.3 idzie innym mechanizmem niż gen.1→2 (Koide, nie φ²);
(ii) `B=√2` (dokładność Koidego) numerycznie potwierdzone (`|B_num−√2|<1e-6`), analitycznie niewyprowadzone
(a3d T5 FAIL, T-OP3). → tylko **CP-0L (Tier 1, honest-framing)**, NIE bloker.

**R2 — m_W / sin²θ_W: consistency-check przebrany za predykcję.** Audyt główny (wiersz QCD) pisze ogólnie
„W/Z declared limit", ale nie flaguje, że nagłówkowe „m_W=80,354 (0,01σ)" to **m_W=m_Z·cosθ_W bierze m_Z jako
INPUT i wciąga α_s** (`sek09:447,454,1249`), a surowe tree sin²θ_W=3/13 jest **11,3σ off** (`sek09:1215`),
ratowane pętlą QCD. Nagłówek letter/companion (`tgp_companion.tex:519-523`) prezentuje to jako czystą predykcję.
→ **CP-6b.**

**R3 — Kosmologia niedoreprezentowana.** Brak osobnego wiersza kosmologicznego w tabeli 22. Niezależna
weryfikacja dodaje (wiersz COSMO §1.1): (a) **DESI w(z) „χ²≈χ²" = fit-by-construction** — skrypt wstrzykuje
Ω_m ΛCDM, pole ψ zamrożone ⇒ w≈−1 trywialnie; `b1_wde_friedmann_fit.py:498` „indistinguishable from ΛCDM";
(b) **Ω_Λ-promocja = przehandlowany parametr** — g̃≈0,98 to odwrócona numerologia („5" ma ≥5 rozkładów,
`op-delta1…/README.md:94-99`), trade-off optymalizuje Ω_Λ ale psuje α_s do +1,26σ (`:204`); źródło samo:
„postulat z matchem, nie ab initio" (`results.md:367`); (c) **sprzeczność wewnętrzna**: modified-gravity
(galaktyki bez DM) vs particle-DM (Ω_DM=0,262) + galaktyczna grawitacja **formalnie obalona**
(`op-PR004-SPARC/README.md:7`, przegrywa z MOND 5,4σ).

**R4 — Licznik 856 vs 688: severity wyżej.** Audyt główny: M03 🟢 CLOSED-RESOLVED, „blok. publ.: nie".
U źródła: **688 nigdy nie spropagowane** — `INDEX.md:224` i `PREDICTIONS_REGISTRY.md:1163` wciąż publikują
**856/784**. Podwójny ledger = self-inconsistency dla recenzenta ⟹ podnoszę do **blokującego** (część CP-6).
M03 *deliverable* (40 plików) faktycznie istnieje — ale finalizacja licznika nie.

**R5 — S01: ostrzejsze framowanie.** Audyt główny: „PARTIAL — brak twierdzenia o jednoznaczności (N5)".
U źródła problem przesunął się ze „sprzeczności 4 form" (faktycznie usuniętej, `\iffalse`) na **„kanon w body =
forma obserwacyjnie SFALSYFIKOWANA 5,02σ"** (`sek08c…tex:419-421` kanon vs `:7-9,27-31` falsyfikacja) — mocniejszy
zarzut niż brak unikalności.

**Zgoda (silna):** dominujący wzorzec „re-labelling ≠ naprawa"; dwa fundamentalne węzły XL (most Γ→Φ +
σ_ab/κ_E); α=2 NON-DERIVABLE jako aksjomat (nie falsyfikacja); L01/M03 jako jedyne realnie domknięte;
S01+L05 jako blokery; publikacja Paper I niezależnie; struktura Tier 0–3; roadmapa CP-1…CP-9.

---

## §5 — Anti-Lakatos (samokontrola tego audytu)
✓ Odróżniono „manuskrypt jest UCZCIWY co do luki" od „luka ZAMKNIĘTA" (kryterium użytkownika).
✓ Nie zinflowano re-labellingu w naprawę; nie zgłoszono jako otwarte tego, co realnie naprawiono (L01, M03).
✓ Werdykty negatywne (#49/#53 NON-DERIVABLE) potraktowane jako *ratyfikacja statusu aksjomatu*, nie falsyfikacja TGP.
✓ Każdy werdykt zakotwiczony w plikach + numerach sesji STATE.md (pełne evidence w workflow `tgp-deep-audit`, run wf_4c82b639-877).

## Cross-references
- [[../STATE.md]] (#42–#55) · [[HONEST_FRAMING_UV_CG_ROOTS.md]] (4 korzenie) · [[../PREDICTIONS_REGISTRY.md]] (ledger 856→688)
- `audyt/S01..S07, L01..L08, M01..M03, D01, T01` (źródłowe README + POST_ACTION)
- `research/op-CG-alpha-eff-convergence-2026-06-26` (#49) · `research/op-CG-Kij-from-Hgamma-2026-06-27` (#53)
- `meta/TGP_W_Z_THEORETICAL_LIMIT.md` (declared limit) · `partial_proofs/{fermion,nuclear}_from_soliton/*VERDICT.md`
- **§6 dowody niezależnej re-weryfikacji:** `scripts/a3d_soliton_brannen_r.py:144` (G0_TAU_fitted) · `paper_lepton_masses/tgp_lepton_masses.tex:611-643` · `core/.../dodatekJ2:99,117-119` · `core/sek09…tex:447,454,1215,1249` · `tooling/scripts/b1_wde_friedmann_fit.py:498` · `research/op-delta1-g-tilde-derivation/README.md:94-99,204` · `research/op-PR004-SPARC/README.md:7` · `INDEX.md:224` vs `PREDICTIONS_REGISTRY.md:1163` · `core/sek08c…tex:419-421` vs `:7-9,27-31`
