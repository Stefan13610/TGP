# Narzędzia i toolchain dla TGP — 2026-04-14

> Przegląd narzędzi, pluginów, modeli i programów które mogą
> usprawnić pracę nad programami badawczymi R1–R7.

---

## 1. Systemy algebry symbolicznej (CAS)

### Obecnie używane
- **Python + SymPy** — w skryptach walidacyjnych (scripts/, nbody/)
- **NumPy / SciPy** — obliczenia numeryczne, ODE, optymalizacja

### Rekomendowane do dodania

| Narzędzie | Do czego | Dla którego problemu | Koszt |
|-----------|----------|---------------------|-------|
| **SageMath** | Zunifikowany CAS (integruje SymPy, GAP, Singular, Maxima) | R2 (teoria grup), R3 (GL(N,𝔽₂)), R6 (analiza ODE) | Darmowy |
| **Cadabra2** | Algebra tensorowa, obliczenia w OW/GR, notacja LaTeX-native | **R4 (metryka z Einsteina)** — idealny do G_μν z g_ij=Φ^p | Darmowy |
| **Mathematica / Wolfram** | CAS ogólnego zastosowania, mocne w analitycznych ODE | R5 (k=4 z wiriału), R6 (B=√2) | Płatny (~$395/rok) |
| **FORM** | Wielkoskalowe obliczenia algebraiczne (QFT) | R7 (UV completion, β-functions) | Darmowy |
| **xAct** (Mathematica) | Pakiety do rachunku tensorowego | R4 (jeśli Mathematica dostępna) | Wymaga Mathematica |

**Konkretna rekomendacja:**
- **R4 (metryka)**: Zainstaluj Cadabra2. Pozwala wpisać g_ij = Φ^p·δ_ij w notacji LaTeX i obliczyć G_μν, warunki ghost-free, energię — dokładnie to czego potrzebujesz.
  ```bash
  # Instalacja (conda)
  conda install -c conda-forge cadabra2
  # Lub z pip
  pip install cadabra2
  ```
- **R3 (dlaczego N=3)**: SageMath ma wbudowane GL(n,GF(2)), pozwala eksplorować algebraiczną strukturę dla różnych N.
  ```python
  # SageMath
  G = GL(3, GF(2))
  print(G.order())  # 168
  print(G.conjugacy_classes_representatives())
  ```

### Linki
- [SageMath](https://www.sagemath.org/) — darmowy, interfejs Pythona
- [Cadabra2](https://cadabra.science/) — field theory CAS, notacja LaTeX
- [SymPy ODE docs](https://docs.sympy.org/latest/modules/solvers/ode.html)
- [FORM](https://github.com/vermaseren/form) — high-energy physics algebra

---

## 2. Dowodzenie twierdzeń (theorem provers)

### Lean 4 + Mathlib + PhysLean

**Dla problemów:** R2 (CG-1/3/4), R5 (m∝A⁴), R6 (B=√2)

Lean 4 to najszybciej rozwijający się system dowodzenia twierdzeń. Kluczowe:
- **Mathlib**: >200 000 sformalizowanych twierdzeń, >750 kontrybutorów
- **PhysLean**: nowa biblioteka formalizująca fizykę w Lean 4 (analog Mathlib dla fizyki)
- W 2025 użyto Lean do znalezienia błędu w opublikowanym artykule o potencjale 2HDM

**Jak to pomaga TGP:**
- CG-1 (kontrakcja Banacha): Lean ma formalizację analizy funkcjonalnej
- B=√2: Sformalizowanie dowodu uodparnia go na błędy
- Publikowalność: artykuł z formalnym dowodem w Lean ma wyższy prestiż

**Jak zacząć:**
1. [Theorem Proving in Lean 4](https://leanprover.github.io/theorem_proving_in_lean4/) — tutorial
2. [Mathematics in Lean](https://leanprover-community.github.io/mathematics_in_lean/) — kurs
3. [PhysLean](https://physlean.com/) — fizyka w Lean

### AI-wspomagane dowodzenie

| Narzędzie | Model | Co robi | Jak pomaga |
|-----------|-------|---------|------------|
| **Goedel-Prover-V2** | 8B / 32B (open-source) | Automatyczne dowodzenie w Lean 4 | Sugeruje taktyki, generuje dowody |
| **Kimina-Prover** | 72B | Reasoning-driven proving w Lean 4 | Skaluje z rozmiarem modelu |
| **LeanDojo** | Framework | AI-driven exploration w ekosystemie Lean | Integruje LLM z Lean |

**Konkretna rekomendacja:**
Zainstaluj Lean 4 + Mathlib. Użyj Goedel-Prover-V2-8B (działa lokalnie na GPU 24GB)
do eksploracji dowodów CG-1. Nawet jeśli nie zamknie problemu, pomoże sformułować
go precyzyjnie.

```bash
# Instalacja Lean 4
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
# Nowy projekt z Mathlib
lake init my_project math
```

### Linki
- [Lean 4](https://lean-lang.org/)
- [Goedel-Prover-V2](https://github.com/Goedel-LM/Goedel-Prover-V2) — HuggingFace models
- [LeanDojo](https://leandojo.org/) — AI + Lean
- [PhysLean](https://physlean.com/) — fizyka w Lean 4

---

## 3. Lokalne modele LLM (do rozumowania matematycznego)

### Najlepsze modele do matematyki/fizyki (2026)

| Model | Rozmiar | MATH benchmark | RAM | Zastosowanie |
|-------|---------|----------------|-----|-------------|
| **DeepSeek-R1-0528** | 671B (MoE, ~22B aktywne) | Najlepszy open-source | 48GB+ VRAM | Rozumowanie krok-po-kroku, fizyka |
| **Qwen3-235B-A22B** | 235B (MoE) | Top 3 | 48GB+ VRAM | Research, scientific reasoning |
| **Phi-4** | 14B | 80.4% | 16GB | Najlepszy stosunek jakość/RAM |
| **QwQ-32B** | 32B | Top 3 math | 24GB | Math reasoning, chain-of-thought |
| **DeepSeek-R1-Distill-32B** | 32B | Bardzo dobry | 24GB | Dystylat z R1, dostępniejszy |

**Konkretna rekomendacja:**
- **Jeśli masz GPU 24GB** (RTX 4090 / 3090): **QwQ-32B** lub **DeepSeek-R1-Distill-32B**
- **Jeśli masz GPU 16GB** (RTX 4080): **Phi-4** (14B)
- **Jeśli masz GPU 8GB** (RTX 4070): **Llama 3.3 8B** (ogólny) lub **Goedel-Prover-V2-8B** (dowody)

**Narzędzia do uruchamiania lokalnie:**
| Program | Opis | Platform |
|---------|------|----------|
| **Ollama** | Najprostszy — `ollama run qwq` | Win/Mac/Linux |
| **LM Studio** | GUI + API, łatwy wybór modelu | Win/Mac/Linux |
| **llama.cpp** | Najszybszy inference (GGUF) | Wszędzie |
| **vLLM** | Serwer API, batching | Linux (GPU) |

```bash
# Szybki start z Ollama
curl -fsSL https://ollama.com/install.sh | sh
ollama run qwq           # 32B reasoning model
ollama run deepseek-r1    # jeśli masz dużo VRAM
```

### Jak używać z TGP
1. **Sprawdzanie dowodów**: Wklej szkic dowodu → model sprawdza logikę
2. **Eksploracja algebraiczna**: "Pokaż jak GL(3,F₂) działa na przestrzeni CKM"
3. **Debugowanie ODE**: "Dlaczego to ODE daje B≈√2? Jakie symetrie ma?"
4. **Pisanie LaTeX**: Generowanie szkiców sekcji z notatek

---

## 4. Pluginy Obsidian

### Już w użyciu (prawdopodobnie)
- Standardowy renderer LaTeX (MathJax)

### Rekomendowane do zainstalowania

| Plugin | Co robi | Dlaczego przydatny |
|--------|---------|-------------------|
| **LaTeX Suite** | 100+ snippetów, szybkie wpisywanie równań | Pisanie notatek z fizyką 3x szybciej |
| **LaTeX Math** | Zintegrowany CAS (SymPy) w notatkach | `$x^2 + 3x + 2 =$` → automatycznie liczy wynik |
| **Theorem & Equation Referencer** | Indeksowanie twierdzeń i równań z cross-reference | Twierdzenia TGP jako baza danych z linkami |
| **Calctex** | Auto-obliczanie formuł LaTeX inline | Szybka weryfikacja numeryczna w notatkach |
| **Dataview** | Zapytania SQL-like na frontmatter | Status tracking problemów R1–R7 |
| **Kanban** | Tablice kanban | Wizualizacja postępu badawczego |
| **Git** (Obsidian Git) | Auto-commit vault | Backup notatek |
| **LaTeX OCR** | Obraz równania → LaTeX | Skanowanie z papierowych notatek |
| **Excalidraw** | Rysunki diagramów | Schematy dowodów, flow teorii |

**Konkretna rekomendacja:**
1. **LaTeX Suite** — od razu. Pisanie `@a` zamienia się w `\alpha`, `//` w `\frac{}{}`.
2. **Theorem & Equation Referencer** — kluczowy. Import z LaTeX (arXiv papers), export do LaTeX.
   Oznaczasz twierdzenia w notatkach → automatyczny indeks → cross-referencje.
3. **LaTeX Math (SymPy)** — inline CAS. Wpisujesz równanie → od razu wynik.

---

## 5. Narzędzia LaTeX i publikacji

### Kompilacja i edycja

| Narzędzie | Opis | Zastosowanie |
|-----------|------|-------------|
| **VS Code + LaTeX Workshop** | Lokalny edytor z preview, synctex | Edycja .tex plików TGP |
| **Typst** | Nowoczesna alternatywa LaTeX (200-400ms kompilacja vs 2-4min) | Szybkie drafty, notatki |
| **WebLaTeX** | VS Code + web + Git + Copilot w jednym | Edycja online z GitHub sync |
| **latexdiff** | Diff między wersjami .tex | Śledzenie zmian w manuskrypcie |
| **Overleaf CE** (self-hosted) | Własna instancja Overleaf | Jeśli chcesz kolaborację |

### CI/CD (już masz GitHub Actions)

| Narzędzie | Co dodaje |
|-----------|----------|
| **tectonic** | Szybsza kompilacja LaTeX (Rust-based) — zastępuje pdflatex |
| **arXiv submission validator** | Sprawdza zgodność z wymaganiami arXiv przed submisją |

---

## 6. Narzędzia numeryczne (rozszerzenie obecnego stacku)

### Dla ODE solitonowych (R5, R6)

| Narzędzie | Co daje | Kiedy użyć |
|-----------|---------|------------|
| **DifferentialEquations.jl** (Julia) | Najszybszy solver ODE/BVP na świecie | Jeśli Python za wolny dla parameter sweeps |
| **AUTO-07p** | Analiza bifurkacji, kontynuacja | Mapowanie przestrzeni rozwiązań ODE solitonowego |
| **PyDSTool** | Dynamical systems toolkit (Python) | Analiza stabilności, bifurkacje |
| **Dedalus** | Spektralne PDE solver (Python) | Jeśli potrzebujesz PDE (nie tylko ODE) |

**Konkretna rekomendacja dla R6 (B=√2):**
AUTO-07p pozwala śledzić rozwiązania ODE solitonowego jako funkcję parametrów (g₀, α, d).
Można zobaczyć **dokładnie** gdzie B=√2 pojawia się w przestrzeni parametrów.

```bash
# AUTO-07p
pip install auto-07p
# Lub Julia (szybsze ODE)
curl -fsSL https://install.julialang.org | sh
julia -e 'using Pkg; Pkg.add("DifferentialEquations")'
```

### Dla teorii grup (R3)

| Narzędzie | Co daje |
|-----------|---------|
| **GAP** | System algebry grup (dedykowany). GL(n,GF(2)), klasy koniugacji, podgrupy |
| **SageMath** (zawiera GAP) | Interfejs Pythona do GAP |

```python
# GAP (standalone)
G := GL(3, GF(2));
Size(G);                          # 168
ConjugacyClasses(G);              # klasy
NormalSubgroups(G);                # podgrupy normalne
```

---

## 7. Workflow: jak to wszystko połączyć

### Proponowany stack

```
┌─────────────────────────────────────────────────┐
│  OBSIDIAN (notatki + LaTeX Suite + Theorem Ref) │
│  ↕ wikilinki do plików .tex i .py               │
├─────────────────────────────────────────────────┤
│  VS Code + LaTeX Workshop (edycja .tex)         │
│  + GitHub Copilot (inline suggestions)          │
│  + Claude Code (agentic tasks, refactoring)     │
├─────────────────────────────────────────────────┤
│  Python (SymPy + SciPy + NumPy) — obecny stack  │
│  + Cadabra2 (tensory, R4)                       │
│  + SageMath/GAP (grupy, R3)                     │
├─────────────────────────────────────────────────┤
│  Lean 4 + Mathlib + Goedel-Prover (dowody, R2)  │
├─────────────────────────────────────────────────┤
│  Ollama + QwQ-32B (lokalny LLM do rozumowania)  │
├─────────────────────────────────────────────────┤
│  GitHub Actions (CI: testy + LaTeX → PDF)       │
│  + Zenodo (archiwizacja DOI)                    │
└─────────────────────────────────────────────────┘
```

### Co zainstalować w jakiej kolejności

| Priorytet | Narzędzie | Czas instalacji | Dla problemu |
|-----------|-----------|-----------------|-------------|
| 1️⃣ | LaTeX Suite (Obsidian plugin) | 2 min | Codzienne notatki |
| 2️⃣ | Cadabra2 (`pip install cadabra2`) | 5 min | R4 (metryka) |
| 3️⃣ | Ollama + QwQ-32B | 10 min + download | Rozumowanie lokalne |
| 4️⃣ | SageMath | 15 min | R3 (GL(n,𝔽₂)) |
| 5️⃣ | Lean 4 + Mathlib | 30 min | R2 (dowody formalne) |
| 6️⃣ | Theorem & Equation Referencer (Obsidian) | 2 min | Indeks twierdzeń |
| 7️⃣ | AUTO-07p | 5 min | R6 (bifurkacje ODE) |
| 8️⃣ | Goedel-Prover-V2-8B (HuggingFace) | download | AI-assisted proofs |

---

## 8. Podsumowanie: narzędzie → problem

| Problem | Narzędzie #1 | Narzędzie #2 | Narzędzie #3 |
|---------|-------------|-------------|-------------|
| **R1** Cabibbo | SageMath (GAP: GL(3,𝔽₂)) | SymPy (algebraiczne) | QwQ-32B (eksploracja) |
| **R2** CG-1/3/4 | **Lean 4 + Mathlib** | Goedel-Prover | SageMath (analiza funk.) |
| **R3** N=3 | **GAP / SageMath** | QwQ-32B (heurystyki) | Lean 4 (jeśli dowód) |
| **R4** Metryka | **Cadabra2** | SymPy (weryfikacja) | xAct (alternatywa) |
| **R5** m∝A⁴ | SymPy (analityczne ODE) | AUTO-07p (kontynuacja) | QwQ-32B (intuicja) |
| **R6** B=√2 | **AUTO-07p** (bifurkacje) | SymPy (symetrie Lie) | Lean 4 (formalizacja) |
| **R7** UV | **FORM** (β-functions) | SymPy (running) | — |

---

> *Dokument utworzony 2026-04-14. Aktualizować w miarę zmian toolchainu.*

## Źródła

- [SageMath](https://www.sagemath.org/)
- [Cadabra2](https://cadabra.science/)
- [Lean 4](https://lean-lang.org/)
- [PhysLean](https://physlean.com/)
- [Goedel-Prover-V2 (GitHub)](https://github.com/Goedel-LM/Goedel-Prover-V2)
- [LeanDojo](https://leandojo.org/)
- [FORM](https://github.com/vermaseren/form)
- [AUTO-07p](https://github.com/auto-07p/auto-07p)
- [Ollama](https://ollama.com/)
- [LaTeX Suite (Obsidian)](https://github.com/artisticat1/obsidian-latex-suite)
- [Theorem & Equation Referencer](https://github.com/RyotaUshio/obsidian-latex-theorem-equation-referencer)
- [LaTeX Math (Obsidian/SymPy)](https://github.com/zarstensen/obsidian-latex-math)
- [WebLaTeX](https://github.com/sanjib-sen/WebLaTex)
- [DeepSeek-R1](https://huggingface.co/deepseek-ai/DeepSeek-R1)
- [QwQ-32B](https://huggingface.co/Qwen/QwQ-32B)
- [GAP System](https://www.gap-system.org/)
- [Best LLM for Math 2026](https://pricepertoken.com/leaderboards/math)
- [Best Open Source LLM for Scientific Research 2026](https://www.siliconflow.com/articles/en/best-open-source-llm-for-scientific-research-academia)
