# Narzędzia i toolchain dla TGP — 2026-04-14

> **Cel dokumentu**: Inwentaryzacja dostępnych narzędzi dla agenta (Claude Code)
> i użytkownika. Aktualizować przy każdej zmianie toolchainu.
>
> **Ostatnia weryfikacja**: 2026-04-14

---

## 0. Status instalacji (AKTUALNE)

### ✅ ZAINSTALOWANE I DZIAŁAJĄCE

| Narzędzie | Wersja | Ścieżka / import | Zastosowanie |
|-----------|--------|-------------------|-------------|
| **Python** | 3.14 | `C:\Users\Mateusz\AppData\Local\Programs\Python\Python314\python.exe` | Główny język |
| **SymPy** | 1.14.0 | `import sympy` | CAS symboliczny |
| **NumPy** | 2.4.3 | `import numpy` | Obliczenia numeryczne |
| **SciPy** | 1.17.1 | `import scipy` | ODE, optymalizacja, algebra |
| **matplotlib** | 3.10.8 | `import matplotlib` | Wykresy |
| **pandas** | 3.0.1 | `import pandas` | Tabele danych |
| **galois** | 0.4.10 | `import galois` | **Ciała skończone, GL(n,GF(2))** — R1, R3 |
| **networkx** | 3.6.1 | `import networkx` | Teoria grafów, struktury grupowe |
| **emcee** | 3.1.6 | `import emcee` | MCMC sampling |
| **Lean 4** | 4.29.1 | `C:\Users\Mateusz\.elan\bin\lean.exe` | Dowodzenie twierdzeń — R2 |
| **Lake** | 5.0.0 | `C:\Users\Mateusz\.elan\bin\lake.exe` | Build system Lean |
| **elan** | — | `C:\Users\Mateusz\.elan\` | Lean toolchain manager |
| **Cadabra2** | 2.5.15 | CLI: `C:\Program Files (x86)\Cadabra\cadabra2-cli.exe` | Algebra tensorowa — R4 |
| **Git** | — | w PATH | Wersjonowanie |
| **GitHub Actions** | — | `.github/workflows/` | CI: testy + LaTeX |
| **Claude Code** | — | Agent w terminalu | Agentic coding, refactoring |

### ⚠️ ZAINSTALOWANE, OGRANICZENIA

| Narzędzie | Problem | Workaround |
|-----------|---------|------------|
| **Cadabra2 Python** | `import cadabra2` nie działa (initialization failed) | Używać CLI: `cadabra2-cli.exe < skrypt.cdb` |
| **Lean 4** | Nie w PATH basha (Claude Code) | Wywoływać pełną ścieżką lub via Python subprocess |

### ❌ NIEDOSTĘPNE (ograniczenia sprzętowe)

| Narzędzie | Powód | Alternatywa |
|-----------|-------|-------------|
| Lokalne LLM (QwQ-32B, DeepSeek-R1) | Za mało VRAM | Claude Code (API), ChatGPT |
| Goedel-Prover-V2-8B | Za mało VRAM | Ręczne dowodzenie w Lean + Claude |
| PyDSTool | Nie wspiera Python 3.14 (numpy.distutils) | SciPy.integrate + SymPy |

---

## 1. Jak agent (Claude Code) powinien używać narzędzi

### Python — standardowy workflow
```python
# Dostępne bezpośrednio:
import sympy, numpy, scipy, galois, networkx, matplotlib, pandas, emcee
# Uruchamianie skryptów:
python TGP/TGP_v1/scripts/nazwa_skryptu.py
python TGP/TGP_v1/nbody/examples/ex200_xxx.py --quick
```

### galois — GL(n, GF(2)) do R1 i R3
```python
import galois
import numpy as np
GF2 = galois.GF(2)
# Tworzenie macierzy nad GF(2)
M = GF2([[1,0,1],[1,1,0],[0,1,1]])
det = int(np.linalg.det(M.view(np.ndarray).astype(float))) % 2
# Enumeracja GL(3,GF(2)): 168 elementów
```

### Cadabra2 — algebra tensorowa do R4
```bash
# Uruchamianie skryptu Cadabra2 (.cdb):
"C:\Program Files (x86)\Cadabra\cadabra2-cli.exe" skrypt.cdb

# Konwersja Cadabra → Python:
"C:\Program Files (x86)\Cadabra\cadabra2python.exe" skrypt.cdb > skrypt.py

# Konwersja Cadabra → LaTeX:
"C:\Program Files (x86)\Cadabra\cadabra2latex.exe" skrypt.cdb > output.tex
```

Przykład skryptu `.cdb` dla R4 (metryka):
```
{r, t, \theta, \phi}::Coordinate;
{a, b, c, d}::Indices(values={r, t, \theta, \phi});
\Phi::Depends(r);
g_{a b}::Metric;
g^{a b}::InverseMetric;

# Ansatz: g_ij = Phi^p * delta_ij
# Oblicz tensor Ricciego, warunki ghost-free, etc.
```

### Lean 4 — dowodzenie twierdzeń do R2
```bash
# Lean via Python subprocess (bo nie w PATH basha):
python -c "
import subprocess
lean = r'C:\Users\Mateusz\.elan\bin\lean.exe'
lake = r'C:\Users\Mateusz\.elan\bin\lake.exe'
# Sprawdzenie pliku .lean:
r = subprocess.run([lean, 'plik.lean'], capture_output=True, text=True)
print(r.stdout, r.stderr)
"

# Nowy projekt Lean z Mathlib:
# (uruchomić w PowerShell/cmd, nie w bash Claude Code)
# lake init tgp_proofs math
# lake build
```

---

## 2. Systemy algebry symbolicznej (CAS)

### Zainstalowane
- **SymPy 1.14.0** — pełne CAS w Pythonie (ODE, algebra, calculus)
- **Cadabra2 2.5.15** — algebra tensorowa, field theory (CLI)
- **galois 0.4.10** — ciała skończone, GL(n,GF(q))

### Do rozważenia w przyszłości

| Narzędzie | Do czego | Dla którego problemu | Koszt | Priorytet |
|-----------|----------|---------------------|-------|-----------|
| **SageMath** | Pełny CAS + GAP (grupy) w jednym | R2, R3, R6 | Darmowy | ŚREDNI |
| **Mathematica** | Mocne analityczne ODE | R5, R6 | ~$395/rok | NISKI |
| **FORM** | Wielkoskalowe obliczenia QFT | R7 | Darmowy | NISKI |

---

## 3. Dowodzenie twierdzeń

### Zainstalowane
- **Lean 4.29.1** + **Lake 5.0.0** — proof assistant
- Toolchains: `leanprover--lean4---v4.29.0`, `v4.29.1`

### Jak zacząć (dla użytkownika)
1. [Theorem Proving in Lean 4](https://leanprover.github.io/theorem_proving_in_lean4/)
2. [Mathematics in Lean](https://leanprover-community.github.io/mathematics_in_lean/)
3. [PhysLean](https://physlean.com/) — fizyka w Lean 4

### AI-wspomagane dowodzenie (przez API, nie lokalnie)
- **Claude** (obecny agent) — może pomagać formułować twierdzenia w Lean
- **Goedel-Prover-V2** — wymaga VRAM, ale można użyć przez API/cloud
- **LeanDojo** — framework integrujący LLM z Lean

---

## 4. Lokalne LLM — NIEDOSTĘPNE

Za mało VRAM na lokalny model. Alternatywy:
- **Claude Code** (obecny agent) — do rozumowania, pisania kodu, eksploracji
- **Claude API** (cloud) — dłuższe sesje rozumowania
- **ChatGPT / Gemini** (cloud) — drugie opinie, weryfikacja

---

## 5. Pluginy Obsidian (rekomendowane)

| Plugin | Co robi | Status |
|--------|---------|--------|
| **LaTeX Suite** | Snippety: `@a`→`\alpha`, `//`→`\frac{}{}` | REKOMENDOWANY |
| **Theorem & Equation Referencer** | Indeks twierdzeń + cross-ref + import/export LaTeX | REKOMENDOWANY |
| **LaTeX Math (SymPy)** | Inline CAS: `$2+3=$` → automatycznie oblicza | OPCJONALNY |
| **Calctex** | Auto-obliczanie formuł LaTeX | OPCJONALNY |
| **Dataview** | Zapytania na frontmatter (status R1–R7) | REKOMENDOWANY |
| **Kanban** | Tablice kanban dla postępu | OPCJONALNY |
| **Excalidraw** | Diagramy, schematy dowodów | OPCJONALNY |
| **LaTeX OCR** | Obraz → LaTeX | OPCJONALNY |

---

## 6. Narzędzia numeryczne

### Zainstalowane (Python)
- **SciPy** `solve_ivp`, `solve_bvp` — ODE/BVP solver
- **SciPy** `optimize` — fsolve, minimize, root
- **NumPy** — algebra liniowa, FFT
- **galois** — arytmetyka GF(2), macierze nad ciałami skończonymi
- **networkx** — grafy, struktury kombinatoryczne

### Do rozważenia

| Narzędzie | Co daje | Dla problemu | Priorytet |
|-----------|---------|-------------|-----------|
| **DifferentialEquations.jl** (Julia) | Najszybszy solver ODE/BVP | R5, R6 (jeśli Python za wolny) | NISKI |
| **AUTO-07p** | Kontynuacja numeryczna, bifurkacje | R6 (B=√2 w przestrzeni param.) | ŚREDNI |

---

## 7. Workflow agenta — szybka referencja

```
┌─────────────────────────────────────────────────┐
│  OBSIDIAN (notatki + pluginy)                   │
│  ↕ wikilinki do plików .tex i .py               │
├─────────────────────────────────────────────────┤
│  Claude Code (agent: edycja, testy, git)        │
├─────────────────────────────────────────────────┤
│  Python 3.14                                    │
│  ├── SymPy 1.14   (symboliczne)                 │
│  ├── NumPy 2.4    (numeryczne)                  │
│  ├── SciPy 1.17   (ODE, optymalizacja)          │
│  ├── galois 0.4   (GL(n,GF(2)), ciała)          │
│  ├── networkx 3.6 (grafy)                       │
│  ├── matplotlib   (wykresy)                     │
│  └── emcee        (MCMC)                        │
├─────────────────────────────────────────────────┤
│  Cadabra2 2.5.15 (CLI — algebra tensorowa)      │
│  Ścieżka: "C:\Program Files (x86)\Cadabra\"    │
├─────────────────────────────────────────────────┤
│  Lean 4.29.1 + Lake 5.0.0                      │
│  Ścieżka: C:\Users\Mateusz\.elan\bin\          │
├─────────────────────────────────────────────────┤
│  GitHub Actions (CI) + Zenodo (DOI)             │
└─────────────────────────────────────────────────┘
```

## 8. Narzędzie → problem (lookup table)

| Problem | Narzędzie #1 | Narzędzie #2 | Narzędzie #3 |
|---------|-------------|-------------|-------------|
| **R1** Cabibbo | **galois** (GL(3,GF(2))) | SymPy (algebraiczne) | Claude (eksploracja) |
| **R2** CG-1/3/4 | **Lean 4 + Mathlib** | SymPy (analiza funk.) | Claude (formułowanie) |
| **R3** N=3 | **galois** (GL(N,GF(2))) | networkx (struktury) | Lean 4 (jeśli dowód) |
| **R4** Metryka | **Cadabra2** (tensory) | SymPy (weryfikacja) | — |
| **R5** m∝A⁴ | SymPy (analityczne ODE) | SciPy (numeryczne) | Claude (intuicja) |
| **R6** B=√2 | SciPy (ODE sweeps) | SymPy (symetrie Lie) | Lean 4 (formalizacja) |
| **R7** UV | SymPy (β-functions) | SciPy (running) | — |

---

> *Ostatnia aktualizacja: 2026-04-14. Zaktualizować po instalacji nowych narzędzi.*

## Źródła

- [SageMath](https://www.sagemath.org/)
- [Cadabra2](https://cadabra.science/)
- [Lean 4](https://lean-lang.org/)
- [PhysLean](https://physlean.com/)
- [Goedel-Prover-V2](https://github.com/Goedel-LM/Goedel-Prover-V2)
- [LeanDojo](https://leandojo.org/)
- [galois (Python)](https://github.com/mhostetter/galois)
- [FORM](https://github.com/vermaseren/form)
- [AUTO-07p](https://github.com/auto-07p/auto-07p)
- [LaTeX Suite (Obsidian)](https://github.com/artisticat1/obsidian-latex-suite)
- [Theorem & Equation Referencer](https://github.com/RyotaUshio/obsidian-latex-theorem-equation-referencer)
- [Best LLM for Math 2026](https://pricepertoken.com/leaderboards/math)
