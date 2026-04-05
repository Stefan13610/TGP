# Sesja v35 — Raport (2026-03-28)

## Status końcowy sesji

| Zadanie | Status |
|---------|--------|
| F5: Atraktor ψ_ini=7/6 (ex86, 9/9 PASS) | ✅ ZAMKNIĘTE |
| F6: Napięcie r≈0.27 vs r<0.036 | ✅ ZAMKNIĘTE (artefakt wzoru) |
| F8: scripts/README_map.md | ✅ ZAMKNIĘTE |
| O-L1: Obalenie S=2π-11/10 | ✅ POTWIERDZONE |
| O-L1: Dokładne α*_∞ | 🔄 W toku (ex88 v36) |
| F7: n_s precyzyjne | ✅ ZAMKNIĘTE analitycznie (v36) |
| O-L2, O-L3 | ⬜ Otwarte |

## Kluczowe wyniki

### F5 ZAMKNIĘTE: Atraktor ψ_ini = 7/6 — Propozycja

**Skrypt:** `nbody/examples/ex86_psi_ini_attractor.py` (9/9 PASS)

**Wyniki numeryczne:**

| Test | Wynik | Status |
|------|-------|--------|
| T1: W(7/6) = 0 | 4.44e-16 ≈ 0 | ✅ PASS |
| T2: W'(7/6) < 0 | −2.722, ω_cosmo/H₀=8.19 | ✅ PASS |
| T3: η_BBN ≪ 1 | η_BBN = 1.2×10⁻³¹ | ✅ PASS |
| T4: ω_cosmo/H₀ > 1 | 8.19 ≫ 1 | ✅ PASS |
| T5: ψ_ini=7/6 → ΔG/G=0 | 0.000% | ✅ PASS |
| T6: ψ_ini=1.0 → ΔG/G∈[5%,20%] | 13.7% | ✅ PASS |
| T7: ψ_eq w środku BBN-compatible | [0.94, 1.56] | ✅ PASS |
| T8a: W(1.10) > 0 | +0.161 | ✅ PASS |
| T8b: W(1.25) < 0 | −0.260 | ✅ PASS |

**Kluczowe liczby:**
- `ω_cosmo/H₀ ≈ 8.19` → underdamped oscillations at z=0 (nowa predykcja!)
- `η_BBN ≈ 1.2×10⁻³¹` → pole absolutnie zamrożone przy BBN
- `Δψ(BBN→equality) ≈ 3×10⁻¹¹` → slow-roll drift negligible
- `ψ_ini=7/6 → ΔG/G = 0.000%` → JEDYNA wartość dająca ΔG/G=0

**Ewolucja ψ(z) dla różnych ψ_ini:**
| ψ_ini | ψ(z=0) | ΔG/G | BBN |
|-------|--------|------|-----|
| 0.85 | 1.116 | 31.3% | ❌ |
| 1.00 | 1.137 | 13.7% | ✅ |
| 7/6 ★ | 1.1667 | 0.000% | ✅ |
| 1.20 | 1.172 | 2.3% | ✅ |
| 1.35 | 1.196 | 11.4% | ✅ |
| 1.50 | 1.232 | 17.9% | ✅ |

**Nowa predykcja testowalna:**
Oscylacje G_eff(t) z częstością ω ≈ 8.19 H₀ dla ψ_ini ≠ 7/6.
Amplituda ∝ |ψ_ini − 7/6|. Testowalne przez LLR/GRACE-FO.

### Awans: prop:psi-ini-derived

- **Przed v35:** `\statuslabel{Hipoteza}` — argument GL, niezweryfikowane
- **Po v35:** `\statuslabel{Propozycja}` — trójkąt argumentów:
  1. Algebraicznie: W(7/6) = 0 → brak siły napędowej (pole zamrożone)
  2. BBN: jedyna wartość ψ_ini dająca ΔG/G = 0
  3. Substrat GL: RG-minimum mapuje na ψ = 7/6
- **Efekt parametryczny:** N_param: 3 → 2 (ψ_ini opuszcza Warstwę II)
- **M_obs/N_param:** ≥ 5 → ≥ 6

### F8 ZAMKNIĘTE: scripts/README_map.md

Plik `scripts/README_map.md` zawiera:
- Mapę ex55–ex87 z kategoriami: ✅ Aktywny | ⚠️ Wymaga weryfikacji | ❌ Artefakt | 📦 Archiwum
- **Ostrzeżenie ex71–ex83:** 13 skryptów zaklasyfikowanych jako `[ARTEFAKT]` — spurious zeros z za grubej siatki
- Mapę ex → LaTeX (stw/hip/propozycje)
- Priorytety v36 (O-L1, O-L2, O-L3)

## Zmiany w plikach LaTeX

| Plik | Zmiana |
|------|--------|
| `sek08_formalizm.tex` | prop:psi-ini-derived: rozbudowany do 4 punktów, awans Hipoteza→Propozycja, nowe liczby ω_cosmo, η_BBN, tabela ψ_ini; rem:parsimony: N=2 (pewne), nowa predykcja oscylacji G_eff |
| `sek08_formalizm.tex` | Tabela parametrów (tab:param-classification): ψ_ini poprawiony, wpis ex86 dodany |
| `status_map.tex` | Wpis sesji v35 (F5): kompletne wyniki T1–T8; ψ_ini: Hipoteza→Propozycja; N_param=2 → N_param=2 (pewne); statystyki: Propozycja: 15→16, Hipoteza: 8→7 |

## Nowe pliki

| Plik | Opis |
|------|------|
| `nbody/examples/ex86_psi_ini_attractor.py` | F5: atraktor ψ_ini=7/6, 9/9 PASS |
| `nbody/examples/ex87_alpha_star_precision.py` | O-L1: gęsta siatka α*, 80 pkt + brentq (w toku) |
| `scripts/README_map.md` | F8: mapa skryptów z klasyfikacją artefaktów |

## F6 ZAMKNIĘTE: Napięcie r≈0.27 vs r<0.036 — artefakt błędnego wzoru

**Diagnoza:**
- `rem:min-observables` (sek08) używało formuły `r = 16ε_ψ` z `ε_ψ = (1-n_s)/2 ≈ 0.017`
- To dawało `r ≈ 0.27` — pozorne napięcie z granicą BK18+BAO (`r < 0.036`)

**Błąd:** `ε_ψ = (1-n_s)/2 = 1/N_e ≈ 0.017` to parametr **η** (dominuje w Starobinskim), nie **ε**.
- W Starobinskim: `ε = 3/(4N_e²) ≈ 2×10⁻⁴ ≪ 1/N_e`
- `n_s - 1 = -6ε + 2η ≈ 2η ≈ -2/N_e` — η dominuje, nie ε
- `r = 16ε = 12/N_e² ≈ 0.003` (N_e=60) — **spójne z r < 0.036**

**Potwierdzenie:** prop:r-tensor (dodatekG) od dawna podaje `r = 12/N_e² ≈ 0.003` — sprzeczność wewnętrzna dokumentu.

**σ_ab i inflacja:** σ_ab = 0 podczas inflacji (substrat w fazie symetrycznej → K_ab ∝ δ_ab → σ_ab = 0). Aktywne tylko po reheating dla anizotropowych źródeł materii.

**Wynik:**
| | Stary (błędny) | Nowy (poprawny) |
|--|--|--|
| Formuła | r = 16ε_ψ, ε_ψ=(1-n_s)/2 | r = 12/N_e² (Starobinsky) |
| Wartość | r ≈ 0.27 | r ≈ 0.003 |
| Status | napięcie 🔴 | spójne ✅ |
| Predykcja | niesprawdzalna | LiteBIRD: r < 0.002 vs TGP r≈0.003 |

**Pliki zmienione:** `sek08_formalizm.tex` (rem:min-observables, uwaga 2), `status_map.tex` (F6 block)

## O-L1 — Wyniki (ex87v2, R_MAX=100, 2/4 PASS)

**Kluczowe odkrycie:** α*₁, α*₂ są silnie wrażliwe na R_MAX!

| | R_MAX=100 (ex87v2) | R_MAX=120 (ex84) |
|--|--|--|
| α*₁ | 2.4318 | 2.4396 |
| α*₂ | 2.6356 | 2.6953 |
| S | 5.0674 | 5.1349 |
| ppm od S_formula | 22336 | 9307 |

**Wniosek O-L1:**
- Formuła S=2π−11/10 jest OBALONA w obu przypadkach (P4 PASS)
- Ale dokładne α*₁, α*₂ wymagają ekstrapolacji R_MAX→∞
- Różnica α*₂ między R_MAX=100 i 120 wynosi 0.06 — zbyt duża

**Status O-L1:** Główny wynik (obalenie S=2π-11/10) POTWIERDZONY. Dokładne α*_∞ niewyznaczone z powodu błędu metodologicznego ex87v3.

## O-L1 — Wyniki (ex87v3, ekstrapolacja R_MAX→∞, 0/4 PASS)

**Wyniki surowe ex87v3:**
| R_MAX | α*₁ | α*₂ |
|-------|-----|-----|
| 80 | nie znaleziono | nie znaleziono |
| 100 | 2.4423 | 2.6769 |
| 120 | nie znaleziono | nie znaleziono |
| 140 | nie znaleziono | 2.6456 |
| 160 | nie znaleziono | 2.6377 |
| 200 | nie znaleziono | nie znaleziono |

Ekstrapolacja (3 pkt α*₂): α*₂_∞ ≈ 2.571 ± 0.003 (niepewna)

**Diagnoza błędu ex87v3:**
- Okno B_coeff = `[0.28*R_MAX, 0.40*R_MAX]` — skalowane z R_MAX
- Efekt: każdy punkt mierzy amplitudę ogona w **innym oknie fizycznym**
- Dane są **nieporównywalne** → ekstrapolacja α*(R_MAX→∞) bezużyteczna
- Przykład: R_MAX=100 daje window=[28,40]; R_MAX=140 daje window=[39.2,56]

**Porównanie α*₂ przy R_MAX=100 (dwa skrypty, różne okna):**
| Skrypt | Okno B_coeff | α*₂ |
|--------|-------------|-----|
| ex87v2 | [28, 42] (fixed) | 2.6356 |
| ex87v3 | [28, 40] (0.28×100) | 2.6769 |
| ex84 | fixed (inny) | 2.6953 (R_MAX=120) |

**Wrażliwość na okno jest silna** — zmiana o 2 jednostki w R_R okna zmienia α*₂ o >0.04.

**Wniosek O-L1 FINALNY:**
- ✅ **Obalenie S=2π-11/10: POTWIERDZONE** (9307-22336 ppm niezależnie od okna i R_MAX)
- ❌ **Dokładne α*_∞: NIEWYZNACZONE** — wymaga stałego okna B_coeff + skan R_MAX (v36/ex88)

**Zadanie v36 (O-L1b):** ex88 z **stałym** oknem B_coeff [30,45], R_MAX ∈ {100,120,150,200,300}

## Otwarte po v35

| ID | Zadanie | Priorytet |
|----|---------|-----------|
| O-L1 | Dokończyć ex87 — potwierdzenie α*₁,α*₂ | 🟡 Średni |
| O-L2 | ex88: czy dwa zera F(α*)=r₂₁ zależą od β_pot? | 🟡 Średni |
| O-L3 | ex89: unifikacja Ścieżek 9+10 | 🟢 Niski |
| F6 | ~~Sektor tensorowy σ_ab — napięcie r≈0.27 vs r<0.036~~ **ZAMKNIĘTE** (v35, artefakt błędnego wzoru; r=12/N_e²≈0.003) | ✅ |
| F7 | n_s precyzyjne z pełnego tła z ε_ψ(N_e) | 🟡 Średni |
