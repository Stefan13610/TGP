---
title: "Phase FINAL — CLOSED-RESOLVED: (B) REFUTED-SUBSTRATE. Substrat NIE generuje α=2 — blokowo-uśredniony kompozyt Φ=⟨σ²⟩ daje K(Φ)∝Φ^{−1} (wykładnik e≈−1, analityka exact; e_inf=−0.12 MC, R²_FSS=0.73), NIE Φ^{+4} (α=2). Niespójność lematu A3 (#31) POTWIERDZONA realna. α=2 ściśle aksjomatyczne-na-gęstości (ratyfikuje status_map l.72 + #48). Ścieżka substratowa do α=2 = CLOSED-NEGATIVE."
date: 2026-06-26
cycle: op-CG-alpha-eff-convergence-2026-06-26
parent: "[[./README.md]]"
phase: FINAL
classification: SUBSTRATE_DERIVATION_STATUS_RESOLVED
verdict: REFUTED-SUBSTRATE
folder_status: closed-resolved
sympy_pass: "3/3 (analityka)"
mc_fss: "e_inf=-0.116, R2_FSS=0.729, 4 L"
hardcoded: 0
anti_lakatos_lock: PRESERVED
---

# Phase FINAL — Closure (B) REFUTED-SUBSTRATE

## §0 — VERDICT
```
████████████████████████████████████████████████████████████████████
█  op-CG-alpha-eff-convergence-2026-06-26                          █
█  SUBSTRATE_DERIVATION_STATUS_RESOLVED — (B) REFUTED-SUBSTRATE    █
█                                                                  █
█  Metryka: abar = alpha_eff(FSS) = e/2.  Manuskrypt: alpha=2 (e=4)█
█                                                                  █
█  ANALITYKA (sympy, exact, 3/3): Phi=sigma^2 (composite <s^2>) ⟹  █
█    chain-rule K(Phi)=1/(4Phi) ~ Phi^{-1}  ⟹  e=-1, alpha_eff=-1/2 █
█    alpha=2 wymaga p=1/6 (Phi=sigma^{1/3}) — bez sensu substrat.  █
█    Escape przez wymiar anomalny: Delta_e=5 (eta~O(5)) — niemozliwe█
█                                                                  █
█  NUMERYKA (phi^4 Z2, FSS L={16,24,32,40}): e_inf=-0.116,         █
█    R2_FSS=0.729 (clean), spread=0.014 ⟹ abar=-0.058              █
█                                                                  █
█  |abar-2| = 2.06 ≥ 1.0 (prog B) ∧ R2_FSS=0.73 ≥ 0.7 ⟹ (B).      █
█                                                                  █
█  ⇒ Substrat NIE generuje alpha=2. Niespójność A3 (#31)          █
█    POTWIERDZONA realna. alpha=2 ŚCIŚLE AKSJOMATYCZNE-NA-GĘSTOŚCI █
█    (ratyfikuje status_map l.72 + op-Kgeo #48). Ścieżka          █
█    substratowa do alpha=2 = CLOSED-NEGATIVE.                     █
████████████████████████████████████████████████████████████████████
```

## §1 — Ustalenia

1. **Rdzeń analityczny rozstrzyga wprost (sympy 3/3, circularity-free):**
   - **T1:** dla kompozytu `Φ = σ²` (= `⟨ŝ²⟩`, substrate composite) reguła łańcuchowa daje
     **`K(Φ) = 1/(4Φ) ∝ Φ^{−1}`** ⟹ wykładnik `e = −1`, `α_eff = −1/2`. To **dokładnie** CG34
     „K_1 ∼ 1/Φ" (#31). Manuskrypt wymaga `K(φ)=K_geo·φ⁴` ⟹ `e=+4`, `α=2`.
   - **T2:** ogólnie `Φ = σ^{2p}` ⟹ `e(p) = 1/p − 2`. Composite `p=1` ⟹ `e=−1`. `α=2` (`e=4`)
     wymagałby `p = 1/6` (`Φ = σ^{1/3}`) — mapowanie **bez jakiegokolwiek uzasadnienia
     substratowego** (`Φ=⟨ŝ²⟩` jest kwadratem, p=1, nie pierwiastkiem szóstego stopnia).
   - **T3 (zamknięcie escape route):** jedyna droga `e: −1 → +4` to wymiar anomalny `Δe = 5`,
     wymagający `η ~ O(5)`; WF 3D `η ≈ 0.036` (3 rzędy za mało). Coarse-graining unitarnego
     near-Gaussian substratu **NIE może** zmostkować `−1 → +4`.

2. **Numeryka potwierdza znak i branżę (FSS, NIEpatologiczny φ⁴ Z₂):** estymator chain-rule
   (NIE artefakt-prone log-log gęstości gradientowej z ex200), L ∈ {16,24,32,40}, `⟨Φ⟩≈0.72`
   stabilne (okno scale-separated istnieje — brak patologii runaway/frozen z CG34). Wynik:
   `e_inf = −0.116`, `R²_FSS = 0.729` (clean, per-L R² rośnie monotonicznie 0.564→0.805 z L),
   `spread = 0.014`. `ᾱ = −0.058`, `|ᾱ−2| = 2.06 ≥ 1.0` ∧ `R²_FSS ≥ 0.7` ⟹ **(B)**.

3. **Rozbieżność MC (−0.12) vs analityka (−1) — uczciwie udokumentowana:** estymator MC binuje
   `K_pt = (∇σ)²/(∇Φ)²` względem `Φ` w węźle; różnica centralna `∇(σ²) = (σ²_{x+1}−σ²_{x−1})/2`
   ≈ `4σ̄∇σ` (σ̄ = średnia sąsiadów) dekoreluje z `σ²` węzła ⟹ **systematyczny bias estymatora
   w stronę 0** (lattice locality artifact). ⟹ MC daje **dolne oszacowanie |e|** (−0.12), exact
   to −1. **Werdykt (B) jest robustny niezależnie** — oba (−1 i −0.12) są ≫ od progu `|e−4|`;
   żadne honest oszacowanie nie zbliża się do `e=+4` (α=2). MC NIE jest cytowane jako „potwierdza
   −1", lecz jako „potwierdza znak/branżę: małe/ujemne, NIE +4".

4. **Konsekwencja:** thm:D-uniqueness / thm:alpha2 ustala **FORMĘ** `K(φ)=K_geo·φ^{2α}` i (z C1-C3)
   `α=2` jako **selekcję w klasie konforemnej na gęstości** — ale substrat `⟨ŝ²⟩` *sam z siebie*
   daje `α_eff = −1/2`, nie 2. ⟹ niespójność lematu A3 (#31, dotąd „zgłoszona, nie rozstrzygnięta")
   jest **POTWIERDZONA jako realna**. α=2 jest **nieredukowalnie aksjomatyczne** (jak K_geo #48, c₀ #37).

## §2 — Spójność z manuskryptem i precedensem (KRYTYCZNE: (B) ≠ falsyfikacja TGP)

- **(B) RATYFIKUJE istniejący uczciwy status, nie obala teorii.** `status_map` l.72 już głosi:
  `K(φ)=K_geo·φ⁴` to „**selekcja w klasie konforemnej (na gęstości), NIE derywacja z substratu**".
  #48 (op-Kgeo) ustalił K_geo jako „aksjomatyczny prefaktor". Ten cykl **domyka pętlę od strony
  numeryczno-substratowej**: pokazuje, że substrat *nie tylko nie pinuje K_geo* (#48), ale *daje
  wręcz inny wykładnik* (`α_eff=−1/2`, nie 2). α=2 MUSI być (i jest) niesione jako aksjomat selekcji.
- **Wzorzec domu (czwarty korzeń):** α=2 dołącza explicite — **od strony substratu, z werdyktem
  CLOSED-NEGATIVE** — do rodziny `irreducibly conditional/axiomatic`: α=2 (#36/#38/#39 + TEN cykl),
  c₀ (#37, Ward/UV), 𝒜→α_s (#43), K_geo (#48). Wszystkie zbiegają do tego samego: makro-fenomenologia
  TGP działa, ale kluczowe stałe NIE są derywowane z substratu — są selekcjami aksjomatycznymi.
- **Wzmacnia honest-framing #42:** α=2 potwierdzone jako **aksjomat selekcji** (nie derywacja),
  zgodnie z ledgerem (N_axiom=6 zawiera α=2). Żadnej zmiany N_free; potwierdzenie, nie rewizja.

## §3 — Anti-Lakatos
✓ Werdykt **WYLICZONY** z plomby (value-blind): `|ᾱ−2|=2.06≥1.0 ∧ R²_FSS=0.73≥0.7` ⟹ (B); progi
0.3/1.0/0.7 niezmienione (anti-moving-goalposts). ✓ Wynik **NEGATYWNY zgłoszony wprost** (B realnie
osiągnięty, nie ukryty ani przemianowany na „warunkowy"). ✓ Rozbieżność MC vs analityka (−0.12 vs −1)
**udokumentowana jako bias estymatora**, nie zamieciona; werdykt zakotwiczony w exact analityce +
robustny pod oboma. ✓ Substrat NIEpatologiczny (φ⁴ Z₂, okno scale-separated) — uczciwiej niż
patologiczny `-J(φ_iφ_j)²` (CG34); pokazuje, że nawet *dobry* substrat daje `e≠+4`. ✓ Circularity guard
ENFORCED (zero α_s/mas kwarków; zakaz odwrotnego dopasowania do α=2). ✓ 0 hardcoded, 0 nowych stałych.
✓ #31/#43/#48 IMMUTABLE, zero re-litygacji. ✓ (B) jawnie odróżnione od falsyfikacji TGP.

## §4 — Dyspozycja (user-gated)
| Cel | Akcja | Status |
|---|---|---|
| dodatekQ2 lemat A3 (#31 „do dopięcia") | reframe: niespójność `α_eff=s−1` (substrat) vs α=2 (aksjomat) **POTWIERDZONA realna**; A3 NIE jest twierdzeniem derywującym α=2 z substratu | 📋 user-gated |
| dodatekQ CG-4 status | `[CZĘŚCIOWY NUM]` → odnotować: składnik α=2↔K(φ) **rozstrzygnięty NEGATYWNIE** (substrat daje α_eff=−1/2; α=2 aksjomatyczne) | 📋 user-gated |
| thm:alpha2 / status_map l.72 | wzmocnić framing: „selekcja na gęstości" potwierdzona numeryczno-substratowo; usunąć ewentualny over-claim „derywacja z substratu" gdziekolwiek występuje | 📋 user-gated |
| #42 ledger (α=2 aksjomat) | **potwierdzony** (N_axiom=6 bez zmian) | ✅ |
| STATE.md | wpis #49 (Phase 1 + FINAL) | ✅ |

## §5 — Następny krok
Cykl **rozstrzygnął** jeden z czterech wspólnych korzeni UV/CG (α=2-from-substrate = CLOSED-NEGATIVE).
Pozostałe (c₀ #37, 𝒜 #43, K_geo #48) wciąż `POSTULATE-CONDITIONAL` — ich jedyna droga to **pełne
domknięcie mostu Γ→Φ / NGFP** (op-uv-as-ngfp: AS NGFP analitycznie 7/7, ale most substrat→makro
nadal [SZKIC]). **Mapa obstrukcji jest teraz kompletna i spójna:** żaden z czterech korzeni nie jest
derywowany z substratu; wszystkie są selekcjami aksjomatycznymi pending UV/CG. To **uczciwy,
domknięty obraz** dla publikacji (analog wyniku metodologicznego G_SPA=48: negatyw też jest wynikiem).

## §6 — Sign-off
**Cykl:** `op-CG-alpha-eff-convergence-2026-06-26` · **Status:** 🟢 **CLOSED-RESOLVED — (B) REFUTED-SUBSTRATE**
**Closure:** 2026-06-26 (1 sesja: Faza A LOCK + Phase 1 FSS + FINAL). Audit trail: README + Phase0_LOCK
(LOCKED) + `Phase1_fss.py`/`.txt` (IMMUTABLE, 0 hardcoded, sympy 3/3 + MC FSS 4 L).
**Claudian** @ 2026-06-26.

## Cross-references
- [[./README.md]] · [[./Phase0_LOCK.md]] · [[./Phase1_fss.py]] · [[./Phase1_fss.txt]]
- [[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]] (#31 — niespójność A3, teraz POTWIERDZONA)
- [[../op-Kgeo-from-D-uniqueness-2026-06-26/Phase_FINAL_close.md]] (#48 — K_geo aksjomatyczny prefaktor)
- [[../../core/_meta_latex/status_map.tex]] l.72 (selekcja na gęstości, NIE derywacja z substratu)
- [[../op-uv-as-ngfp/]] (NGFP — wspólny korzeń UV/CG pozostałych 3 stałych)
