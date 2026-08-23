---
title: "Phase 0 LOCK — op-CG-Kij-from-Hgamma: balance sheet + pre-rejestracja F-CGK-A..D (value-blind) + forbidden moves"
date: 2026-06-27
type: phase0-lock
status: "LOCKED (pre-rejestracja; immutable po commit). Aktywacja Phase 1 wymaga user 'działaj'."
cycle: op-CG-Kij-from-Hgamma-2026-06-27
parent: "[[README.md]]"
pre_registration_date: "2026-06-27"
locked_inputs:
  - "#49 op-CG-alpha-eff-convergence: e=-1 (alpha_eff=-1/2), IMMUTABLE"
  - "#39 op-bond-order-RG-selection: alpha=2 RG-irrelevant, IMMUTABLE"
  - "status_map l.72: alpha=2 = selekcja na gestosci, NIE derywacja"
  - "dodatekB rem:B-v2-status: alpha=2 = bezposredni skutek aksjomatu v2 (bez MK->GL)"
anti_lakatos_note: "Wszystkie progi/reguly zapisane PRZED Phase 1. Werdykt wyliczany z flag. Wynik DERIVABLE pre-akceptowany (nie tylko NON-DERIVABLE). Endpoint +4 nie wbudowany. 0 hardcoded planowane."
---

# Phase 0 LOCK — op-CG-Kij-from-Hgamma

## §1 — Balance sheet (co wchodzi, co jest pytane, co jest deliverable)

### 1.1 — Wejscia LOCKED (NIE podlegaja rewizji w tym cyklu)
| Wejscie | Tresc | Status |
|---|---|---|
| `H_Gamma` (eq:B-H) | `sum_i[(m0^2/2)s^2+(lam0/4)s^4] - J sum_<ij> s_i s_j`, Z2 | aksjomat substratu (dodatekB) |
| Slownik | `phi=(Phi/Phi0)^{1/2}`, `Phi=<s^2>`, `g=phi^2` | dodatekB / sek08b |
| #49 | Gaussian CG: `e=-1` (`K~Phi^-1`); escape `eta~O(5)` zamkniety | LOCKED |
| #39 | `alpha=2` RG-irrelevant | LOCKED |
| T_C (predecessor) | `(phi_i phi_j)^2 ~ s^2 s^2 != s_i s_j` (bilinear) | [[../op-amplitude-density-phase-bridge-2026-06-27/Phase1_Vsub_FINDINGS.md]] |

### 1.2 — Obiekt audytu (HIPOTEZA pod testem, NIE wejscie)
Sprzezenie geometryczne v2 `K_ij=J(phi_i phi_j)^2` → `K(phi)=K_geo phi^4` → `alpha=2`.
**Pytanie:** wyprowadzalne z 1.1 czy nieredukowalny aksjomat?

### 1.3 — Deliverable
Werdykt strukturalny `F-CGK-D in {DERIVABLE, NON-DERIVABLE, UNDETERMINED}` + kompletna mapa
obstrukcji. Claim ceiling: **C (STRUCTURAL_VERIFIED)**. Output type: **structural**. Brak observable.

## §2 — Pre-derywacja (OCZEKIWANIE, NIE prog; zapisane przed Phase 1)

Standardowe wymiary 3D Ising (literatura, do cytowania w Phase 1, nie zakladania):
- `eta ≈ 0.036`, `nu ≈ 0.630`, `[s] = (1+eta)/2 ≈ 0.518`.
- Operator energii `eps ~ s^2` (= nasze `Phi`): `Delta_eps = d - 1/nu ≈ 3 - 1.587 = 1.413`.
- Operator kinetyczny kompozytu `O_kin = (grad Phi)^2`: `Delta_O = 2 Delta_eps + 2 ≈ 4.826`.

**Oczekiwanie (do testu):**
- B: `De = +4 - (-1) = 5` wymaga `eta_O` rzedu 5; bound unitarnosci 3D Z2 daje `eta` male
  (O(0.01-0.1)) ⇒ **B-REFUTED**.
- C1: `Delta_O ≈ 4.83 > d=3` ⇒ `O_kin` **irrelewantny** ⇒ nie pinowany przez FP ⇒ **C-AXIOM**.
- ⇒ oczekiwany agregat: **NON-DERIVABLE** (alpha=2 nieredukowalny aksjomat).

**To jest oczekiwanie, nie wynik.** Reguly §F sa nadrzedne; DERIVABLE pre-akceptowany jesli flagi.

## §3 — Plan Phase 1/2 (sympy substance; >=50% non-trivial)
- **Phase 1:**
  - T-A1 (F-CGK-A): rozwiniecie gradientowe bilinearnego bondu `-J s_i s_j` → `K_s=const`, `e=-1` w `Phi` (re-derywacja value-blind, anchor). [first-principles]
  - T-C1 (F-CGK-C1): wymiar RG `O_kin=Phi^n(grad Phi)^2`, n=0,1,2 przy 3D WF z `Delta_eps` (literatura); relewantny (`Delta<3`) czy irrelewantny (`Delta>3`)? [first-principles + literature-anchored]
- **Phase 2:**
  - T-B1 (F-CGK-B): wymagane `eta_O` dla `De=5` vs bound unitarnosci (cytowany); osiagalne czy nie? [literature-anchored consistency]
  - T-C2 (F-CGK-C2): dodaj czteropolowy bond `-J' sum (s_i s_j)^2` do `H_Gamma`; czy `J'` fiksowany symetria/konsystencja, czy wolny parametr? [first-principles]
- **Structural declarations** (NIE liczone w sympy PASS): S05 single-Phi preserved; scope.

## §F — FALSYFIKATORY (pre-rejestrowane, value-blind, LOCKED) {#F}

| ID | Test | Reguła decyzyjna (kandydat, LOCKED) |
|---|---|---|
| **F-CGK-A** | baseline: bilinearny bond → `e` | `e=-1` ⇒ anchor PASS (sanity #49); `e!=-1` ⇒ STOP (blad rachunku/konwencji, nie werdykt) |
| **F-CGK-B** | `eta_O` wymagane dla `De=5` vs bound unitarnosci 3D Z2 | wymagane `eta_O <= bound` ⇒ **B-DERIVABLE**; `eta_O >> bound` ⇒ **B-REFUTED** |
| **F-CGK-C1** | wymiar RG `O_kin` przy WF | `Delta_O < 3` (relewantny) ⇒ kandydat na pinowanie; `Delta_O > 3` (irrelewantny) ⇒ wklad do **C-AXIOM** |
| **F-CGK-C2** | wspolczynnik czteropolowego bondu | fiksowany (symetria/konsystencja) ⇒ **C-DERIVABLE-CONDITIONAL**; wolny ⇒ **C-AXIOM** (nowy param do #42) |
| **F-CGK-D** | **AGREGAT (z flag)** | `DERIVABLE` jesli (B-DERIVABLE ∨ C-DERIVABLE-CONDITIONAL); **`NON-DERIVABLE`** jesli (B-REFUTED ∧ C-AXIOM); `UNDETERMINED` w pozostalych |

Progi LOCKED: „relewantny" = `Delta_O < d=3` (scisle); „eta duze" = `eta_O > 1` (anti-moving-goalposts:
nawet liberalny prog `eta_O>1` wystarcza, bo wymagane ~5). Werdykt = funkcja flag, 0 hardcoded.

## §G — Forbidden moves (LOCKED)

1. **Zakaz re-litygacji #49 / #39 / status_map l.72.** F-CGK-A re-derywuje `e=-1` WYLACZNIE jako
   value-blind anchor; wynik IMMUTABLE.
2. **Zakaz przemycania `+4` przez niemotywowana redefinicje kompozytu** (`Phi=s^{1/3}`, `p=1/6`).
   Zmiana identyfikacji `Phi=<s^2>` dozwolona TYLKO z niezaleznym uzasadnieniem substratowym.
3. **Zakaz claimu „alpha=2 derived"** chyba ze `F-CGK-D = DERIVABLE` scisle z §F.
4. **Zakaz cichego dodawania parametrow.** Czteropolowy bond `J'` (jesli potrzebny) = NOWY aksjomat;
   wplyw na ledger #42 JAWNIE odnotowany.
5. **Zakaz hardcodowania** endpointu `+4`/`alpha=2`; wszystkie wykladniki liczone.
6. **Zakaz mylenia „generowany pod RG" z „relewantny pod RG".** Operator moze byc generowany,
   ale irrelewantny ⇒ nie pinuje IR ⇒ nie jest derywacja.
7. **Bound unitarnosci `eta` CYTOWANY z literatury** (3D Ising/bootstrap), nie zalozony.
8. **Budzet nowych stalych: 0** (poza jawnie zliczonym `J'` w galezi C2, jesli aktywna).
9. **S05 single-Phi axiom** zachowany bezwarunkowo (forbidden globalny).

## §H — TGP-native + methodology read (pre-Phase-1 checklist)

- [ ] Q1 Pattern coverage: `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §2 (coarse-graining patterns)
- [ ] Q2 Red flags §3 (BD-drift: nie traktowac 3D Ising CFT jako „import", to literatura-anchor wymiarow)
- [ ] Q3 Inherited LOCKs: #49, #39 mapping
- [ ] Q4 Std-physics tools (ERG, wymiar operatora, bootstrap bound) — justify: to NATYWNE narzedzia
      substratu (substrat = teoria pola 3D), nie BD-translation
- [ ] Q6 GR-limit framing: n/d (cykl strukturalny, brak sektora grawitacyjnego)
- [ ] Pre-flight read: [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]], [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]],
      [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4

## §I — Relacja do PR / ledger
- **Brak nowego PR** (output strukturalny; brak observable). Upstream do **PR-025**: jesli
  `F-CGK-D=NON-DERIVABLE`, wzmacnia HONEST_FRAMING (alpha=2 aksjomat) bez zmiany PR-025.
- **Ledger #42:** zmiana TYLKO jesli galaz C2 aktywna i `J'` wolny (nowy aksjomat) — wtedy
  `N_axiom += 1` JAWNIE. W przeciwnym razie ledger bez zmian.

## Status
🔒 **Phase 0 LOCKED 2026-06-27.** Pre-rejestracja kompletna. Phase 1 czeka na user „działaj".

## Cross-references
- [[README.md]] · [[../op-amplitude-density-phase-bridge-2026-06-27/Phase1_Vsub_FINDINGS.md]]
- [[../op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49)
- [[../op-bond-order-RG-selection-2026-06-23/Phase_FINAL_close.md]] (#39)
- [[../../meta/HONEST_FRAMING_UV_CG_ROOTS.md]] · [[../../core/_meta_latex/status_map.tex]] l.72
