---
title: "HANDOFF — prompt dla agenta wykonującego cykl op-galactic-substrate-tail"
date: 2026-06-13
type: meta-handoff
status: READY (cykl REGISTERED-QUEUED; aktywacja = user 'działaj')
cycle: "[[../research/op-galactic-substrate-tail-2026-06-13/]]"
queue_next_after_this: "op-nucleation-dimensionality ([[SCOPING_op-nucleation-dimensionality_2026-06-12.md]]) — decyzja user 2026-06-13"
authored_by: "Claudian, sesja #18 (user: 'rozpisz op-galactic-substrate-tail i prompt dla nowego agenta, ustaw ND jako kolejny cykl')"
---

# HANDOFF PROMPT — skopiuj poniższy blok jako zadanie nowego agenta

---

Jesteś ekspertem fizyki teoretycznej pracującym w ramach TGP (Teoria Generowanej Przestrzeni,
projekt `TGP/TGP_v1/`). Twoim zadaniem jest wykonanie cyklu badawczego
**`op-galactic-substrate-tail-2026-06-13`** (REGISTERED-QUEUED → aktywacja po Phase 0
+ autoryzacji użytkownika).

## 0. Tożsamość metodologiczna (nienegocjowalne)

Reżim anty-Lakatosowski: pre-rejestracja PRZED rachunkiem; falsyfikatory z progami
liczbowymi LOCKED w Phase 0, nigdy nie korygowane ex post; HONEST_NEGATIVE = pełnoprawny
wynik (w tym cyklu wręcz OCZEKIWALNY — patrz fast-kill §3); każda faza wymaga osobnego
„działaj"; 0 hardcoded T_pass; 0 nowych stałych fundamentalnych bez deklaracji; werdykty
poprzedników LOCKED — **w szczególności PR-004 TRIGGERED-FALSIFIED (mechanism) jest
NIETYKALNY: ten cykl NIE jest recovery PR-004, lecz structural-amendment path z jego
własnego kontraktu (`if_recovery_exhausted`); sukces ⇒ NOWY falsyfikator PR, nie korekta
starego.**

## 1. Obowiązkowe lektury (W TEJ KOLEJNOŚCI, przed Phase 0)

1. `TGP/TGP_v1/STATE.md` — wpisy sesji #17 (BA) i #18 (PR-004 execution)
2. `research/op-PR004-SPARC-fit-execution-2026-06-12/` — CAŁY cykl: Phase0 (pipeline
   zero-parametrowy — **do reuse jako silnik porównawczy, LOCKED**), Phase1_fit.py
   (parser SPARC + MOND benchmark + paired 5σ), Phase_FINAL (werdykt + §3.5 wskaźnik
   kierunkowy); dane SPARC lokalnie w `./data/`
3. `meta/PRE_REGISTERED_FALSIFIERS.md` §PR-004 (z EXECUTION UPDATE 2026-06-13) — kontrakt
   i jego dyspozycja struktural-amendment
4. `research/op-CE-H-3D-native-interaction-retry-2026-05-23/` — γ-1 retry CLEAN PASS:
   **natywne oddziaływanie logarytmiczne defektów 3D (współczynnik −2π; R²_log = 0.9998)
   ORAZ konwergencja exp(−m_σL)/L** — napięcie log-vs-ekranowanie = centralne pytanie cyklu
5. `research/op-frontier-bridge-and-asymmetry-2026-06-12/Phase4_derivation.md` §1 —
   spektrum fluktuacji wokół realnego tła: **m_χ²= λ(w²−Φ₀²) ⇒ mod fazowy U(1) DOKŁADNIE
   bezmasowy w nasyconym bulku (Goldstone)** — kandydat na nieekranowany mediator
6. `research/op-frontier-microphysics-2026-06-11/Phase2_derivation.md` §2 linia K3 —
   **LOCKED: F_substrat = 0 z ∇⟨Φ⟩ = 0 (MODUŁ)** — Twój mechanizm musi być jawnie
   rozliczony względem tej granicy (sektor fazowy/winding ≠ tło modułu; sprzeczność
   = HONEST_NEGATIVE, nie reinterpretacja K3)
7. `meta/CALIBRATION_PROTOCOL.md` §3.6 + §3.6.6-§3.6.13 — WIĄŻĄCE (analytical
   pre-derivation, sign conventions, assumptions, ±5%, klasyfikacja stałych)
8. `meta/CYCLE_KICKOFF_TEMPLATE.md` §1-§3 — kontrakt kickoff (L1: output_observable =
   v_rot(R) [km/s], instrument = SPARC — to jest cykl observable-first)
9. `meta/REALITY_CONTACT_AUDIT_2026-06-12.md` — kontekst strategiczny (po co ten cykl)
10. KOLEJKA: `meta/SCOPING_op-nucleation-dimensionality_2026-06-12.md` — następny cykl
    PO tym (decyzja user); NIE aktywować, nie mieszać zakresów

## 2. Cel cyklu (trzy pytania, kolejność wymuszona)

**Q1 (mechanizm — fast-kill first):** czy istnieje NIEEKRANOWANY dalekozasięgowy kanał
oddziaływania międzysolitonowego z akcji TGP? Hipotezy-kandydaci (zbiór CLOSED do
pre-rejestracji): **H-GOLD** (wymiana bezmasowego modu fazowego: solitony z windingiem
U(1)/RP² sprzęgają się do χ; w 3D daje to klasę 1/r lub — dla sprzężeń topologicznych —
log; γ-1 log-form jako wsparcie), **H-SCREEN** (wszystko ekranowane na 1/m_σ ⇒ brak kanału
⇒ HONEST_NEGATIVE — koniec cyklu po Phase 1, uczciwie i tanio). Rozstrzygnięcie: rachunek
zasięgu propagatora sektora fazowego wokół tła Φ₀ + sprzężenie soliton–χ z akcji.

**Q2 (skala, TYLKO jeśli Q1 = H-GOLD):** wyprowadź a₀-analog z stałych TGP — bez
jakiegokolwiek fitu. Kandydat strukturalny (INFORMATIONAL, zapisany przed rachunkiem):
klasa cH = c/t (γ-3 R = ct; marginalność FCR już wiąże dynamikę lokalną z ct przez
GM/(ct)); numerologia a₀_obs ≈ cH₀/2π vs /6 — **a₀_obs WYŁĄCZNIE comparison-only,
NIGDY input** (circularity guard FP w każdym sympy).

**Q3 (test — NOWY PR):** re-run SPARC 175 z mechanizmem z Q1+Q2: **identyczny
zero-parametrowy pipeline LOCKED z op-PR004-execution** (te same Υ, te same filtry, ten
sam paired test, dane lokalne). Pre-rejestracja progów PRZED re-runem (Phase 0!):
kandydat reguły: PASS jeśli χ²_red(TGP-tail) ≤ χ²_red(MOND simple) w teście sparowanym;
FAIL jeśli gorzej o ≥5σ; pasmo pośrednie = PARTIAL. Numer PR przydzieli Phase 0
(kandydat PR-024). **Przegrana = koniec ścieżki bez recovery** (zapisz w kontrakcie).

## 3. Phase 0 — wymagania twarde

1. Falsyfikatory F-GST-A (mechanizm; klasy: H-GOLD_DERIVED / H-SCREEN_NEGATIVE / GAP /
   INDETERMINATE), F-GST-B (skala; klasy: DERIVED / DERIVED_WITH_ANCHOR / GAP), F-GST-C
   (SPARC re-run; progi liczbowe LOCKED), F-GST-D (agregat) — zbiory CLOSED.
2. Analytical pre-derivation (§3.6): zasięg propagatora χ wokół Φ₀ (oczekiwane: bezmasowy
   ⇒ 1/r-klasa; sprawdź czy RP²/kompaktowość nie generuje masy); forma sprzężenia
   soliton–χ z windingu; konwencje znaków; enumeracja założeń.
3. **Fast-kill jako Phase 1** (najtańszy decydujący rachunek — precedens BA P1): jeśli
   ekranowanie/brak sprzężenia ⇒ HONEST_NEGATIVE i ZAMKNIJ cykl (nie ciągnij Q2/Q3).
4. Forbidden moves (min. 12): zakaz modyfikacji PR-004/PR-022/FM/FCR/BA/γ-*/CE-H; zakaz
   ρ_DM (S05); **zakaz użycia a₀_obs, V_obs, jakiejkolwiek danej SPARC w WYPROWADZENIU**
   (dane wyłącznie w fazie porównawczej Q3 — analog zakazu G_obs); zakaz fitowania
   czegokolwiek (pipeline zero-parametrowy); zakaz zmiany pipeline'u porównawczego
   (LOCKED z PR-004-execution — identyczność = wiarygodność porównania); zakaz
   reinterpretacji K3 (sprzeczność z F_substrat = 0 ⇒ NEGATIVE); zakaz nowych stałych
   (budżet 0; λ, Φ₀ symboliczne); zakaz miękkiego domknięcia; zakaz ND-dimensionality
   scope-creep (to następny cykl, osobny).
5. Risk register (min.): R-ekranowanie (HIGH — fast-kill); R-zgodność-K3 (HIGH — sektor
   fazowy musi NIE wnosić siły tła w jednorodnym bulku, a TYLKO międzysolitonową);
   R-2D→3D (γ-1 log był w geometrii wirowej); R-Goldstone-vs-obserwacje (bezmasowy skalar
   ⇒ potencjalne piąte siły/dyspersja GW — sprawdź zgodność z PR-005 boundem Δc/c);
   R-numerologia-a₀ (pokusa dopasowania kombinacji stałych do 1.2×10⁻¹⁰ — value-blind!).
6. Anticipated outcomes (INFORMATIONAL): zapisz kierunki pokusy PRZED rachunkiem
   (sukces = „rozwiązanie ciemnej materii" — najsilniejsza pokusa w historii projektu;
   właśnie dlatego progi muszą być LOCKED zanim policzysz pierwsze χ²).

## 4. Protokół wykonania

- Czekaj na user „działaj" przed Phase 0 LOCK i przed każdą fazą; po każdej fazie raport
  z decision menu.
- Każda faza: PhaseN_derivation.md + PhaseN_sympy.py/.txt (PASS/FAIL per FP, 0 hardcoded,
  circularity-guard FP: free-symbols audit na a₀_obs/V_obs/G_obs).
- Kolejność: Phase 0 → **Phase 1 fast-kill (Q1)** → [jeśli H-GOLD] Phase 2 skala (Q2) →
  Phase 3 SPARC re-run (Q3; reuse Phase1_fit.py z PR-004-execution z podmienionym modelem
  TGP — diff modelu jawny w skrypcie) → Phase FINAL (agregat; decyzja PR — wyłącznie user).
- STATE.md wpis po każdej fazie; README flip parking → active przy aktywacji.
- Po zamknięciu tego cyklu (dowolny werdykt): **następny w kolejce =
  op-nucleation-dimensionality** (scoping gotowy; własny Phase 0; własna autoryzacja).

## 5. Czego absolutnie nie wolno

- Modyfikować PR-004 TRIGGERED-FALSIFIED (nowy mechanizm = NOWY PR; stare χ² Newton-baryon
  pozostają w rejestrze na zawsze)
- Używać a₀_obs/V_obs/SPARC jako inputu wyprowadzenia mechanizmu lub skali
- Dodawać ρ_DM lub jakiegokolwiek pola spoza S05+Z₂+U(1)+RP² (luka = GAP, nie nowe pole)
- Fitować Υ, a₀-analog, czy cokolwiek do krzywych (pipeline zero-parametrowy LOCKED)
- Przeciągać cyklu po negatywnym fast-killu (HONEST_NEGATIVE po Phase 1 = sukces
  metodologiczny, nie porażka — zamknij i przejdź do kolejki)
- Mieszać zakresu z op-nucleation-dimensionality (następny cykl, osobna pre-rejestracja)

---

**Koniec promptu.** Stan rejestracji: cykl REGISTERED-QUEUED (folder + README istnieją);
aktywacja wymaga user „działaj" → Phase 0 LOCK. Kolejka po nim: op-nucleation-dimensionality.
