---
title: "HANDOFF — propagacja user-gated po sesji #40–#43 (sektor kwarkowy + parameter-counting): 3 work-packages (honest-framing / re-scopes / PR-025)"
date: 2026-06-25
type: handoff
status: READY
author: "Claudian (sesja #40–#43)"
target: "nowy agent / nowa sesja"
scope: "WYŁĄCZNIE propagacja JUŻ-ZALOCKOWANYCH werdyktów #40–#43; ZAKAZ re-litygacji"
---

# HANDOFF — propagacja po #40–#43 (3 work-packages)

> **Dla nowego agenta:** to są **edycje propagacyjne** rozstrzygnięć, które zostały już
> wyliczone i zamknięte w cyklach #40–#43 (2026-06-25). **NIE prowadzisz nowej fizyki ani
> nie re-litygujesz werdyktów** — przenosisz je do warstwy rdzenia/submission/rejestru.
> Każdy WP = osobne „działaj" (można rozbić na 3 sesje).

## §0 — Lektury obowiązkowe PRZED czymkolwiek (w tej kolejności)
1. [[STATE.md]] wpisy **#40, #41, #42, #43** (góra pliku) — pełen kontekst kampanii.
2. [[research/op-L08-quark-g0-tail-vs-core-audit-2026-06-25/Phase_FINAL_close.md]] (#40 NORM-OVERLOAD; §3 dyspozycja)
3. [[research/op-quark-mass-core-g0-rescue-test-2026-06-25/Phase_FINAL_close.md]] (#41 RESCUE-CONFIRMED; §3 dyspozycja + Phase1_audit §5 PR-025)
4. [[research/op-parameter-counting-balance-sheet-2026-06-25/Phase_FINAL_close.md]] (#42 HEADLINE-OPTIMISTIC; §3 rekomendowany framing)
5. [[research/op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]] (#43 POSTULATE-CONDITIONAL; §4 dyspozycja)
6. [[meta/PRE_REGISTERED_FALSIFIERS.md]] §0–§2 (format append-only; potrzebne do WP3)
7. [[meta/CYCLE_KICKOFF_TEMPLATE.md]] §2.6 (konwencja build) — jeśli edytujesz .tex

## §1 — Reguły wiążące (wszystkie WP)
- **ZAKAZ re-litygacji:** werdykty #40–#43 są IMMUTABLE. Edytujesz tekst, by ODZWIERCIEDLIĆ
  werdykt, nie by go zmienić.
- **Build .tex:** po każdej edycji pliku wchodzącego do `main.tex` uruchom
  `pdflatex main.tex` **×2**, potwierdź **exit 0** i **~553 str.**; sprawdź **brak NOWYCH
  dangling refs** (7 undefined refs = pre-existing residual #32, NIE z Twoich edycji).
- **Edycje addytywne/minimalne:** preferuj prozę/adnotacje nad zrywaniem `\ref`/`\label`.
- **Ścieżki względne** od roota vault (np. `core/sek08b_ghost_resolution/...`).
- Po każdym WP: wpis do [[STATE.md]] (#44/#45/#46) + raport decision-menu.
- **Anti-Lakatos:** honest-framing #42 obniża pozorną parsymonię — to ma być widoczne, nie ukryte.

---

## WP1 — Propagacja honest-framingu (#42 HEADLINE-OPTIMISTIC)

**Cel:** zastąpić headline „40 predykcji z 3 inputów" uczciwym framingiem z #42
(uczciwy N_free=10, N_axiom=6, vs SM ~19; korona = leptony 1→3).

**Źródło prawdy:** #42 Phase_FINAL §3 (rekomendowany tekst).

**Wzorzec tekstu (do dostosowania per plik):**
> „~30–40 predykcji obserwacyjnych z **~10 wolnych parametrów numerycznych + 6 aksjomatów
> selekcji strukturalnej** (N=3, Z₂, α=2, β=γ, φ-FP, Koide), w porównaniu do ~19 wolnych
> parametrów Modelu Standardowego. Najmocniejszy wynik: **sektor leptonowy — trzy masy
> naładowanych leptonów z jednego sprzężenia g₀^e (via φ-FP + Koide) do 0,006%.**"

**Pliki docelowe (sprawdź aktualne brzmienie przed edycją):**
- `README.md` — sekcja Abstract + Highlights + istniejąca „Parameter-counting note (2026-06-22)"
  (rozszerz ją o wynik #42: konkretny N_free=10 vs „3", N_axiom=6, vs SM 19).
- `tgp_letter.tex` (PRL submission) — zdanie headline „40 predictions from 3 inputs".
- `tgp_companion.tex` (PRD) — analogicznie.
- (opcjonalnie) `core/sek00_summary/sek00_summary.tex` — jeśli głosi „3 inputy" jako fakt.

**Uwaga anti-Lakatos:** NIE chowaj Φ_0/c₀/g̃/quark-anchors. NIE licz Z₂/α=2 jako parametrów
(to aksjomaty, osobno). Konwencja symetryczna z SM (symetrie nieliczone po obu stronach).

**Build:** README = markdown (poza buildem); tgp_letter/tgp_companion/sek00 = `pdflatex` ×2 exit 0.

**Definicja done:** headline we wszystkich 4 miejscach spójny z #42; build exit 0; wpis STATE.

---

## WP2 — Re-scopes (#40 / #41 / #43 + α_s registry)

**Cel:** zaktualizować 4 miejsca, które noszą nieaktualny/przeszacowany status.

### WP2.1 — sek08b:529 (z #40 NORM-OVERLOAD, potwierdzone dodatekX l.1207)
- **Plik:** `core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex` l.528–529.
- **Problem:** „Universalność kwarkowa: ten sam ODE działa na leptony i kwarki
  (g₀ ∈ [0,817; 0,891])" — myląco sugeruje JEDNĄ domenę.
- **Faktycznie (dodatekX res:X-phiFP-universal, l.1207):** [0,817; 0,891] to **per-sektorowe
  kotwice φ-FP bazowego g₀^(1)** {down 0,817; lepton 0,870; up 0,891}, NIE domena hierarchii
  (rdzeniowe g₀ rozpina [g₀^e, g₀^τ]=[0,869; 1,730]).
- **Edycja:** przeformułuj na „bazowe g₀^(1) per sektor leżą w [0,817; 0,891] (kotwice φ-FP);
  hierarchia generowana przez rdzeniowe g₀ ∈ [g₀^(1), φ·g₀^(1)]". Cytuj dodatekX res:X-phiFP-universal.

### WP2.2 — audyt/L08 problem #3 (z #41 RESCUE-CONFIRMED)
- **Plik:** `audyt/L08_kink_fermion_closure/README.md` — problem #3 quark sub-component.
- **Edycja statusu:** `INDETERMINATE-PENDING-RESCUE` → **`RESCUE-CONFIRMED`** (z caveatami:
  m_t scheme-zależny pole 0,8%/MS-bar 5,5%; 𝒜 warunkowe na CG). Link do #41 FINAL.
- **Dodaj notę:** HALT-B (op-L08-Phase6-quark-sector-mass-formula) był strawmanem (0/3 składników
  maszynerii dodatekX + błędna domena g₀ #40); werdykt HALT-B IMMUTABLE, ale RE-SCOPED.

### WP2.3 — dodatekX l.1353 (z #43 POSTULATE-CONDITIONAL)
- **Plik:** `partial_proofs/quark_sector/dodatekX_quark_sector.tex` ~l.1353 (rem:audit-response,
  „Luka istotnie domknięta").
- **Edycja (addytywna adnotacja, NIE usuwanie):** dopisz notę #43: luka NIE domknięta, lecz
  **przesunięta** do postulatu `K_geo·m_sp²=π·Φ₀²` (eq:X-K-msp-hypothesis), który redukuje się do
  **niedomkniętego mostu Γ→Φ** (CG-1/CG-3 [SZKIC], ex200 4/8). 𝒜=C_F²α_s² = consistency-check
  warunkowy; α_s 0,03σ NIE first-principles. Status: POSTULATE-CONDITIONAL (analog α=2/c₀).
- **Build:** dodatekX wchodzi do main.tex → `pdflatex` ×2 exit 0.

### WP2.4 — α_s w PREDICTIONS_REGISTRY (z #42/#43)
- **Plik:** `PREDICTIONS_REGISTRY.md` — wiersz α_s.
- **Edycja:** sklasyfikuj α_s(M_Z)=0,1179 jako **consistency-check warunkowy** (z 𝒜=C_F²α_s²,
  pending CG-1/CG-3), **NIE first-principles predykcja**. Link #43.
- (opcjonalnie) `core/sek07_predykcje/sek07_predykcje.tex` R12: nota „recovery m_b/m_t DZIAŁA
  (0,6%/0,8%); otwarte pełne [AN] dla 𝒜 (CG)".

**Definicja done:** 4 miejsca spójne z #40/#41/#43; build exit 0 (dla .tex); wpis STATE.

---

## WP3 — PR-025 (formalizacja, z #41 — user autoryzował)

**Cel:** dopisać PR-025 do append-only [[meta/PRE_REGISTERED_FALSIFIERS.md]] (NOWY wpis na końcu,
NIGDY modyfikacja istniejących).

**UWAGA anti-Lakatos (krytyczne):** masy kwarków były **znane** (PDG) przy konstrukcji maszynerii
dodatekX — więc m_b/m_t to **NIE czysta pre-rejestracja**, lecz **retrospektywny log + forward
falsyfikator** (precedens: PR-001 był „RETROACTIVE LOG"). Zapisz to UCZCIWIE.

**Szkic wpisu (dostosuj do formatu §1 rejestru):**
```yaml
PR-025 (RETROSPECTIVE LOG + FORWARD FALSIFIER): TGP quark mass reproduction via 𝒜=a_Γ/φ
  Cycles: op-quark-mass-core-g0-rescue-test-2026-06-25 (#41) + op-A-derivation-from-CG (#43)
  Pre-registration: ⚠ NIE czysta — masy PDG znane; retrospektywny log maszynerii dodatekX
  Native observable (retrospektywnie odtworzone): m_b, m_t z samozgodnego domknięcia
    {m_0=𝒜·m_3/m_1 ; Q_K(m_i+m_0)=2/3}, 𝒜=a_Γ/φ uniwersalne
  Wynik (#41): m_b 0,59% ; m_t 0,77% (pole) / 5,5% (MS-bar) ; 𝒜 univ 0,33% ; α_s=√𝒜/C_F 0,03σ
  Status 𝒜 (#43): POSTULATE-CONDITIONAL — α_s consistency-check warunkowy (pending CG-1/CG-3)
  FORWARD FALSIFIER (genuine, pre-registered TERAZ):
    "Jeśli (a) precyzyjne wspólne-schemat (MS-bar μ) m_t odbiega od predykcji >X% z 5σ,
     LUB (b) domknięcie CG-1/CG-3 wymusi K_geo·m_sp²≠π·Φ₀² (𝒜≠C_F²α_s²),
     maszyneria 𝒜=a_Γ/φ=C_F²α_s² FALSIFIED."
  Recovery scope: brak post-hoc tuningu 𝒜/schematu; CG closure = osobny track
  Caveats: m_t scheme-zależny; 𝒜 warunkowe na niedomknięty CG (#43)
```
- **Próg X% (forward):** zaproponuj wartość (np. 5%) i ZALOCKUJ w tym wpisie (immutable).
- **Decyzja LOCK = user** w oryginalnym kontrakcie; tu user już autoryzował WP3 — odnotuj
  autoryzację z datą.

**Definicja done:** PR-025 dopisany (append-only, uczciwy retrospektywny framing + forward
falsyfikator); cross-link do #41/#43; wpis STATE.

---

## §2 — Kolejność rekomendowana
WP2 (re-scopes, najpilniejsze — korygują żywe statusy) → WP1 (honest-framing) → WP3 (PR-025).
Każdy WP osobno; po każdym build (jeśli .tex) + STATE + decision-menu.

## §3 — Czego NIE robić
- NIE re-derywować/re-litygować #40–#43 ani dodatekX/status_map (read-only fizyka).
- NIE domykać CG-1/CG-3 (wieloletni track; poza zakresem).
- NIE liczyć Z₂/α=2 jako wolnych parametrów (aksjomaty, osobno).
- NIE głosić „α_s wyprowadzone" (jest consistency-check warunkowy, #43).
- NIE chować Φ_0/c₀/g̃ w honest-framingu (bias dwustronny).
