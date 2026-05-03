---
title: "FINDINGS — Galaxy Scaling — galaktyka jako płaska studnia potencjału"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/galaxy_scaling
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — Galaxy Scaling — galaktyka jako płaska studnia potencjału

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## TL;DR / Wynik / Verdict sections (cytaty)

### `README.md` — Wyniki gs9a-gs9b: przejście wymiarowe 3D→2D

> ### Wyniki gs9a-gs9b: przejście wymiarowe 3D→2D
> 
> **Struktura matematyczna (gs9a)**:
> - Ogon solitonu zanika jak r^(-(d-1)/2) — potwierdzone numerycznie
> - Informacja Fishera g'²/g jest konforemnie niezmiennicza **tylko** w d=2
> - ALE: konf. inw. termu kinetycznego nie pomaga (term potencjałowy dominuje przy dużych r)
> - Kluczowe: Poisson (∇²) vs Helmholtz (∇²+1) — sprężyna wymusza oscylacje
> 
> **Propaga...(truncated)

### `README.md` — Wyniki gs9c: porównanie z SPARC (2693 punkty danych)

> ### Wyniki gs9c: porównanie z SPARC (2693 punkty danych)
> 
> Dane: Lelli, McGaugh, Schombert, Pawlowski (2017) — RAR z 175 galaktyk.
> 
> **Ranking χ²/N (z optymalizacją a₀)**:
> 
> | # | Model | χ²/N | a₀ (best-fit) | RMS (dex) |
> |---|-------|------|----------------|-----------|
> | 1 | **McGaugh formula** | **2.72** | 1.16×10⁻¹⁰ | 0.1328 |
> | 2 | **MOND simple** | **2.74** | 1.14×10⁻¹⁰ | 0.1328 |
> | 3 | Hybrid...(truncated)

### `README.md` — Wyniki gs9d: mechanizm przejścia wymiarowego

> ### Wyniki gs9d: mechanizm przejścia wymiarowego
> 
> **MECHANIZM ŚREDNIEJ GEOMETRYCZNEJ**:
> 
> Substrat TGP ma dwie skale:
> - **r_S = GM/c²** — skala źródła (promień Schwarzschilda)
> - **r_H = c/H₀** — skala kosmologiczna (horyzont korelacji)
> 
> Efektywna „głębokość" interakcji grawitacyjnej = średnia geometryczna:
> ```
> H = √(r_S × r_H) = √(GM/(c·H₀)) = r_MOND / √(2π)
> ```
> 
> - r ≪ H: pełne 3D → Newton: F ~ GM/...(truncated)

### `README.md` — Wyniki gs9e: konfrontacja z danymi i rewizja

> ### Wyniki gs9e: konfrontacja z danymi i rewizja
> 
> **KRYTYCZNY TEST: a₀(z)**
> 
> | Obserwacja | a₀=const | a₀~H(z) | Status |
> |------------|----------|---------|--------|
> | Lokalny BTFR | OK | OK | oba OK |
> | SPARC RAR | OK | OK | oba OK |
> | Milgrom z~2 (a₀<4×) | OK | ~3× | MARGINALNY |
> | DLA0817g z=4.26 | OK | ~7× | **NAPIĘCIE** |
> 
> **a₀ ~ H(z) jest DISFAVOROWANE** przez BTFR non-evolution (McGaugh 20...(truncated)

### `README.md` — Wyniki gs10: HARD NUMERICS — galaxy-by-galaxy fitting (171 galaktyk)

> ### Wyniki gs10: HARD NUMERICS — galaxy-by-galaxy fitting (171 galaktyk)
> 
> **Dane**: 171 galaktyk SPARC, 3375 punktów danych, indywidualne krzywe rotacji.
> **Metoda**: differential evolution + L-BFGS-B polish na każdej galaktyce.
> 
> **Phase 1 — Individual fits (a₀, Yd, shape free per galaxy)**:
> 
> | Model | χ²/ndof | median χ²_red | median a₀ (×10⁻¹⁰) | median Yd |
> |-------|---------|---------------|---...(truncated)

### `README.md` — Wyniki gs11: analityka zwycięskiego równania

> ### Wyniki gs11: analityka zwycięskiego równania
> 
> **Równanie**: ν(y) = 1 + exp(-y^α) / y^γ, z α = 4/5, γ = 2/5 (rationalnie)
> 
> **Elegancka forma**: podstawiając u = y^(2/5):
> ```
> ν = 1 + exp(-u²) / u  ≈  1 + √π · erfc(u)    (dla u > 1)
> ```
> **Związek z funkcją erfc** (komplementarna funkcja błędu) — sugeruje mechanizm progowy!
> 
> **Asymptotyki**:
> | Reżim | g_obs | Efektywny wymiar |
> |-------|-------|--...(truncated)

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa