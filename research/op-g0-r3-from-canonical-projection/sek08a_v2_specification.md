---
title: "Sek08a v2.0 specification (G.0 closure draft, NIE actual core mod)"
date: 2026-05-02
phase: 3
sub-task: P31
parent: "[[Phase3_setup.md]]"
predecessors:
  - "[[Phase1_results.md]]"
  - "[[Phase2_results.md]]"
  - "[[phase3_P32_newton_limit_rederivation.txt]]"
status: SPECIFICATION COMPLETE — gotowe do user review
scope_constraint: "DRAFT only. Actual sek08a/sek08a_akcja_zunifikowana.tex modification = Phase 4."
---

# Sek08a v2.0 — full specification draft (G.0 closure)

> **Cel:** Specyfikuj DOKLADNIE jakie zmiany w `core/sek08a_akcja_zunifikowana/`
> sa wymagane po G.0 closure (Phase 1+2+P32).
>
> **Scope:** Ten dokument jest SPECIFICATION + DRAFT TEXT. NIE jest aktualna
> modyfikacja core. Actual mod = Phase 4 po user approval.

---

## 0. Summary changes (TL;DR)

| Element | sek08a v1.x | sek08a v2.0 (G.0) |
|---|---|---|
| Kinetic K(ψ) | ψ⁴ | **ψ⁴** (unchanged) |
| Volume √(-g) | c₀·ψ (M9.1) | **c₀·ψ/(4-3ψ)** (M9.1'') |
| Potential V(ψ) | (β/3)ψ³ - (γ/4)ψ⁴ | **-γ·ψ²(4-3ψ)²/12** (V_M911) |
| Effective U_eff | ψ V/(4-3ψ) - bug | **γ·(ψ⁴/4 - ψ³/3)** |
| Vacuum ψ_vac | claim ψ=1 (lecz derivation gives 16/15) | **ψ=1** (rigorous) |
| m_sp² | claim γ (lecz derivation gives -γ tachion) | **+γ** (rigorous, stable) |
| Φ-EOM (static) | β/γ-dependent EOM (mismatch z R3) | **R3 ODE EXACT** |
| Source coupling | 2q·ρ/Φ_0 | **5q·ρ/Φ_0** |
| Newton-limit | q·c²/Φ_0 = 2πG_0 | **q·c²/Φ_0 = (4/5)πG_0** |
| κ (operational, po re-fit) | 4πG_0/(3H_0²) = 3/(4Φ_0) | **4πG_0/(3H_0²) = 3/(4Φ_0)** (INVARIANT po re-fit) |

**Backwards compatibility:** Wszystkie observable predictions IDENTYCZNE (po re-fit
Φ_0). Tylko parameterization się zmienia.

---

## 1. Sekcje sek08a wymagajace update

Z grep audit (P33) sek08a v1.x ma 4 occurrences V_orig, 38 occurrences kappa.
Update plan poniżej.

### 1.1 prop:K_psi-uniqueness — UNCHANGED ✓

K(ψ) = K_geo · ψ⁴ (T-D-uniqueness, α=2) jest zachowane. Phase 1 G0a confirmed.

### 1.2 NEW: prop:V-M911-canonical (REPLACE prop:V-tgp-canonical)

```latex
\begin{proposition}[Kanoniczny potencjał V_M911 dla M9.1'']
\label{prop:V-M911-canonical}
W TGP-canonical akcji $S_{\rm TGP}$ z metryką M9.1''
($\sqrt{-g}=c_0\psi/(4-3\psi)$) i kinetic $K(\psi)=\psi^4$,
unikalnym (mod stała) potencjałem zapewniającym że static spherically-symmetric
EOM jest tożsame z R3 ODE jest:
\begin{equation}\label{eq:V-M911}
  V_{M911}(\psi) = -\frac{\gamma}{12}\,\psi^2 (4-3\psi)^2
\end{equation}
\end{proposition}

\begin{proof}[Szkic dowodu]
Statyczne EOM po wariacji $\delta S/\delta\psi=0$ z volume $4\pi r^2 c \psi/(4-3\psi)$:
\begin{equation}
\psi'' + \frac{2}{r}\psi' + \frac{2}{\psi}(\psi')^2 = -\frac{1}{K(\psi)}\frac{dU_{\rm eff}}{d\psi}
\end{equation}
gdzie $U_{\rm eff}(\psi) = \psi V(\psi)/(4-3\psi)$.

Wymaganie tożsamości z R3 ODE:
\begin{equation}
\psi'' + \frac{2}{r}\psi' + \frac{2}{\psi}(\psi')^2 = \frac{1-\psi}{\psi^2}
\end{equation}
daje $-U_{\rm eff}'(\psi)/\psi^4 = (1-\psi)/\psi^2$, czyli
$U_{\rm eff}'(\psi) = -\psi^2(1-\psi) = \psi^3 - \psi^2$.
Całkując: $U_{\rm eff}(\psi) = \gamma(\psi^4/4 - \psi^3/3) + C$.
Stąd $V(\psi) = (4-3\psi) U_{\rm eff}/\psi = -\gamma\psi^2(4-3\psi)^2/12$. $\square$
\end{proof}

\begin{remark}[Backwards compatibility z sek08a v1.x]
Stare wyrażenie $V(\psi) = (\beta/3)\psi^3 - (\gamma/4)\psi^4$ było derived
przy zal. $\sqrt{-g}=c_0\psi$ (forma I metryki, M9.1, FALSIFIED 2026-04-25,
Cassini 3$\sigma$). Po canonical M9.1'' $\sqrt{-g}=c_0\psi/(4-3\psi)$
i wymogu R3-equivalence (G.0 closure), poprawne $V$ to $V_{M911}$.

Zauważ: $V_{M911}(1) = -\gamma/12 \neq -\gamma/12$ jak stary $V(1) = \beta/3-\gamma/4$
przy $\beta=\gamma$. Ale efektywny $U_{\rm eff}(1) = -\gamma/12$ jest invariant!
Stąd vacuum energy density (Λ) zachowana — sek05 ciemna energia OK.
\end{remark}
```

### 1.3 NEW: prop:psi-EOM (REPLACE old prop:psi-EOM)

```latex
\begin{proposition}[Statyczne EOM solitonu = R3 ODE]
\label{prop:psi-EOM-R3}
W zunifikowanej akcji $S_{\rm TGP}$ z $K(\psi)=\psi^4$, $V=V_{M911}$,
$\sqrt{-g}=c_0\psi/(4-3\psi)$, statyczne sferycznie-symetryczne równanie pola
po wariacji jest tożsame z R3 ODE:
\begin{equation}\label{eq:psi-EOM-R3}
\psi'' + \frac{2}{r}\psi' + \frac{2}{\psi}(\psi')^2 = \frac{1-\psi}{\psi^2}
\end{equation}
zwane R3 ODE (research/why_n3/, $\alpha=2$).
\end{proposition}

\begin{remark}[Konsekwencje strukturalne]
Z prop.~\ref{prop:psi-EOM-R3} wynika, ze \textbf{wszystkie sukcesy R3} (N=3
generacji, mass spectrum lepton z PDG <0.01\%, Koide K=2/3, spin-1/2,
4-th gen forbidden) są \textbf{automatycznie} sukcesami sek08a v2.0.

Konkretne anchor predictions (Phase 2 P22):
\begin{itemize}
\item $m_\mu/m_e = 206.766$ (PDG 206.768, diff $-0.0013\%$)
\item $m_\tau/m_e = 3477.40$ (PDG 3477.23, diff $+0.0049\%$)
\item Bariera $g_{\rm crit} = 1.874 \equiv \psi_{\rm horizon} = 4/3$ (M9.1'' Lorentzian)
\item 4-th generacja: $g_0^{(4)} = \phi g_0^\tau \approx 2.87 > g_{\rm crit}$ FORBIDDEN
\end{itemize}
\end{remark}
```

### 1.4 UPDATE: prop:vacuum-stability (FIXED — bug w sek08a v1.x)

```latex
\begin{proposition}[Stabilnosc tla prozniowego — G.0 corrected]
\label{prop:vacuum-stability-G0}
W akcji $S_{\rm TGP}$ G.0-canonical (z $V=V_{M911}$, $K=\psi^4$,
$\sqrt{-g}=c_0\psi/(4-3\psi)$):

\begin{enumerate}
\item Vacuum: $\psi_{\rm vac} = 1$ (jedyny realny pierwiastek $U_{\rm eff}'(\psi)=0$
  o additional condition $\psi>0$).
\item Stability: $U_{\rm eff}''(1) = +\gamma > 0$.
\item Yukawa mass: $m_{\rm sp}^2 = U_{\rm eff}''(1)/K(1) = \gamma$ (stable, scalar field).
\end{enumerate}
\end{proposition}

\begin{remark}[Korekta sek08a v1.x bug]
Sek08a v1.x \emph{deklarowało} $\psi_{\rm vac}=1$ i $m_{\rm sp}^2 = \gamma$,
ale derivacja z $V_{\rm orig}$ + $\sqrt{-g}=c_0\psi$ (M9.1 falsified) faktycznie
dawała:
\begin{itemize}
\item $U_{\rm eff,old}(\psi) = \gamma\psi^4(4-3\psi)/12$ z vacuum przy $\psi=16/15$ (NIE 1!)
\item $U''_{\rm eff,old}(1) = -\gamma$ tachionowy bug (Phase 2 P21 uncovered)
\end{itemize}
G.0 closure z $V_{M911}$ + M9.1'' $\sqrt{-g}$ \textbf{naprawia} oba bugi:
vacuum poprawnie przy $\psi=1$ z stable $m_{\rm sp}^2 = +\gamma$.
\end{remark}
```

### 1.5 UPDATE: prop:kappa-corrected (Newton-limit re-derivation)

```latex
\begin{proposition}[$\kappa$ po G.0 re-fit Phi_0]
\label{prop:kappa-corrected-G0}
W zunifikowanej akcji G.0-canonical, operational kappa zachowuje swą wartość:
\begin{equation}\label{eq:kappa-corrected-G0}
\kappa = \frac{q\,c^2}{\Phi_0}\cdot\frac{c_{\rm src}}{3H_0^2} = \frac{4\pi G_0}{3 H_0^2}
       = \frac{3}{4\Phi_0}
\end{equation}
gdzie:
\begin{itemize}
\item $c_{\rm src} = 5$ (source coefficient w G.0, vs $c_{\rm src}=2$ w sek08a v1.x)
\item Po re-fit Newton-limit: $qc^2/\Phi_0 = (4/5)\pi G_0$ (G.0)
       vs $qc^2/\Phi_0 = 2\pi G_0$ (v1.x)
\item Iloczyn $(qc^2/\Phi_0) \cdot c_{\rm src} = 4\pi G_0$ jest INVARIANT
\end{itemize}
\end{proposition}

\begin{remark}[Re-fit Phi_0]
G.0 closure wymaga jednej z dwóch re-calibrations:
\begin{itemize}
\item Jesli $q$ fundamentalne: $\Phi_0^{\rm new} = (5/2) \Phi_0^{\rm old}$
\item Jesli $\Phi_0$ fundamentalne: $q^{\rm new} = (2/5) q^{\rm old}$
\end{itemize}
W obu wypadkach observable Newton's $G_0$, $\kappa$, $|dG/G|/H_0$, BBN, LLR, CMB
SA IDENTYCZNE z sek08a v1.x — re-calibration jest gauge-equivalence.
\end{remark}
```

### 1.6 UPDATE: eq:cosmo-linearized-unified

```latex
% W rozdziale FRW perturbations:
\begin{equation}\label{eq:cosmo-linearized-unified-G0}
\ddot{\delta\psi} + 3H\dot{\delta\psi} + m_{\rm sp}^2 c^2 \delta\psi = -\frac{5q c^2}{\Phi_0}\delta\rho
\end{equation}

% Note: source coefficient = 5 (NIE 2 jak v1.x). To wynika z volume update
% √(-g)=c₀ψ/(4-3ψ).
```

### 1.7 UPDATE: eq:newton-limit

```latex
% Newton-limit relation:
\begin{equation}\label{eq:newton-limit-G0}
\frac{q\,c^2}{\Phi_0} = \frac{4\pi G_0}{5}
\end{equation}

% Comment: faktor 5 w mianowniku odpowiada source-coefficient z G.0
% (vs faktor 2 w v1.x).
```

### 1.8 UPDATE: All Φ-EOM derivations using √(-g)

WSZYSTKIE derivacje Φ-EOM z założeniem √(-g)=c·ψ wymagają re-derivation
z poprawnym √(-g)=c·ψ/(4-3ψ). Output będzie R3 ODE (Phase 1 G0a verified).

Konkretne miejsca (z P33 grep):
- sek08_formalizm L418, L4107: Element objętości
- sek08c_metryka L239, L413: Derivation chain
- sek08a_akcja L616: $\sqrt{-g}\sim e^{4U}$ approximation (też wymaga aktualizacji)

---

## 2. Section header changes

```latex
% Update title:
\title{Sekcja 8a — Zunifikowana akcja TGP-canonical (v2.0, G.0 closure)}

% Add introductory disclaimer:
\begin{remark}[Wersja 2.0 — G.0 closure]
Ta wersja sekcji 8a jest wynikiem programu \textbf{G.0} 
(\href{research/op-g0-r3-from-canonical-projection/README.md}{research/op-g0-r3-from-canonical-projection/}),
w którym ustalono, że pełna formal consistency między TGP-canonical
$\Phi$-EOM a R3 ODE wymaga update potencjału $V$ do $V_{M911}$ oraz
$\sqrt{-g}$ do M9.1'' canonical form. G.0 Phase 1+2+3 (4/4 PASS each)
zweryfikowano że obserwowalne predictions są INVARIANT po re-fit Φ₀.

Sek08a v1.x: deprecated po 2026-05-02, zachowane w git history.
\end{remark}
```

---

## 3. Cross-section impact

### 3.1 sek08c (metric chain) — P34 separate

A1, A2, A3 audit annotations w sek08c stają się closeable:
- A1 (Φ-EOM mismatch): RESOLVED przez G.0 prop:psi-EOM-R3
- A2 (√(-g)=c·ψ obsoleted): RESOLVED przez G.0 √(-g)=c·ψ/(4-3ψ) adopted
- A3 (4 metric forms): RESOLVED — M9.1'' UNIQUE canonical

Patrz `sek08c_A1_A2_A3_closure_draft.md` (P34).

### 3.2 sek08 (formalizm) — main re-write

Sek08 ma 9 V_orig + 107 kappa + 12 vacuum-stability references. Re-write:
- Sekcje "Akcja zunifikowana" — pointer do sek08a v2.0
- prop:vacuum-stability — replace z prop:vacuum-stability-G0
- Wszystkie kappa = 3/(4Φ_0) — pozostaje wartość, update derivation comment

### 3.3 sek09 (cechowanie) — adjust kappa context

Sek09 ma 8 kappa references. Wartość zachowana, ale comment "z prop:kappa-corrected"
może wymagać updated reference do nowego prop:kappa-corrected-G0.

### 3.4 dodatekH (lancuch wyprowadzeń) — chain update

Łańcuch derywacji A11b: $S_{\rm TGP} \to \kappa = 3/(4\Phi_0)$ wymaga
add G.0 closure step (V_M911, sqrt(-g) update, source coef 5).

### 3.5 status_map.tex — version bump

Sek08a status: "v1.x" → "v2.0 (G.0 closure)"
Add G.0 cycle entry.

---

## 4. Backwards compatibility statement

**Sek08a v2.0 jest backwards-compatible w predykcjach z sek08a v1.x** pod
warunkiem re-calibration jednego free parameter (Φ_0 lub q).

Wszystkie observational tests, ktore pass'owały z v1.x, pass'ują z v2.0:
- Newton's G_0 ✓
- PPN γ_PPN=1 (Cassini) ✓
- PPN β_PPN=1 (Mercury) ✓
- BBN |dG/G| ≤ 0.15 ✓
- LLR |dG/G|/H_0 ≤ 0.02 ✓
- CMB n_s, r ✓
- Mass spectrum lepton (PDG <0.01%) ✓ (przez R3 ODE jako effective EOM)
- Koide K=2/3 ✓
- 4-th generation forbidden ✓

**Strukturalne ulepszenia v2.0 vs v1.x:**
- Φ-EOM = R3 ODE (jednoznaczny chain Φ → fermiony)
- Vacuum mass m_sp² = +γ (FIXED tachion bug)
- Vacuum w ψ=1 (FIXED bug w sek08a v1.x derivation)
- N=3 generacji + mass spectrum derivable z fundamental aksjomatu
  (NIE niezalezny formalizm — promoted z structural theorem do TGP-canonical theorem)

---

## 5. Versioning notes

```
sek08a v1.0 (2024-..) — initial action S_TGP
sek08a v1.1 (2025-..) — kappa = 3/(4*Phi_0) corrected
sek08a v1.x (2026-04-..) — kinetic K(ψ)=ψ⁴ T-D-uniqueness lock
sek08a v2.0 (2026-05-..) — G.0 closure: V_M911 + M9.1'' + R3 EOM unification
```

Version bump: minor (1.x) → major (2.0) bo:
- V(ψ) form structurally changed
- √(-g) form structurally changed
- Φ_0 (lub q) re-calibrated (free parameter)
- Φ-EOM identified jako R3 ODE (new theorem)
- Vacuum derivation cleaned (bug fix)

---

## 6. Phase 4 implementation plan

(Phase 4 = actual core mod, separate task po user approval)

### 6.1 Order

1. **Pre-step:** create git branch `sek08a-v2.0-g0-closure`
2. Update `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`:
   - Replace prop:V-tgp-canonical → prop:V-M911-canonical
   - Replace prop:psi-EOM → prop:psi-EOM-R3
   - Replace prop:vacuum-stability → prop:vacuum-stability-G0
   - Update prop:kappa-corrected → prop:kappa-corrected-G0
   - Update eq:cosmo-linearized-unified, eq:newton-limit
   - Add intro remark (v2.0 versioning)
3. Update `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex`:
   - Close A1, A2, A3 (per P34 draft)
   - Update derivation steps using new √(-g)
4. Update `core/sek08_formalizm/sek08_formalizm.tex`:
   - 12 vacuum-stability references → vacuum-stability-G0
   - 9 V_orig references → V_M911 (or pointer to sek08a v2.0)
   - kappa references comment update (value invariant)
5. Update `core/_meta_latex/status_map.tex`:
   - Sek08a v1.x → v2.0
   - Add G.0 cycle entry
6. Update other affected files (sek09, dodatekH, dodatekO, sek07_predykcje, ...).
7. Run latex compile + verify no broken references
8. Commit: "sek08a v2.0 — G.0 closure (V_M911, M9.1'', R3 EOM unification)"

### 6.2 Estimate

- Sek08a main update: 4-6 godzin (z verification of all derivations)
- Sek08c A1/A2/A3 closure: 1-2 godziny
- Sek08 + sek09 + dodatki updates: 4-6 godzin
- Status map + meta: 0.5 godziny
- Verification + recompile: 1-2 godziny

**Total Phase 4 estimate: 11-16 godzin (1-2 dni intense pracy).**

---

**Status:** Specification complete. Ready dla user review + approval do Phase 4.
