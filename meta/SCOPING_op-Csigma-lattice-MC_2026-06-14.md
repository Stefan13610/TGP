---
title: "SCOPING — op-Csigma-lattice-MC: liczbowe wyznaczenie T=C_σσ_0² (tensor stiffness) z kierunkowego bąbla ⟨O_ab O_cd⟩ na sieci 3D Ising (Monte Carlo / FRG), by przypiąć κ_E=8πG_0T/c³ liczbowo i rozstrzygnąć sektor radiacyjny (FALSIFIED / SURVIVE / GAP). Klasa dodatekQ CG-2 (numeryczny), wielosesyjny."
date: 2026-06-14
type: meta-scoping
status: "PRE-PHASE-0 NOTE (wymaga własnego Phase 0 + 'działaj'; zero werdyktów)"
origin: "User 2026-06-14 (sesja #31): 'zarejestruj op-Csigma-lattice-MC' — spawn z op-Csigma-coarse-graining FINAL (residual GAP = liczbowe T)"
parent_cycle: "[[../research/op-Csigma-coarse-graining-2026-06-14/]] (PARTIAL; UNDERDETERMINED-fine-tuned węższy; lean FALSIFIED; residual GAP = wartość T)"
grandparent_cycle: "[[../research/op-sigma-kinetic-Csigma-2026-06-14/]] (UNDERDETERMINED-fine-tuned; C_σ=GAP)"
target_open_problem: "rem:sigma-params (wartość C_σ/T 'niezobliczona') + rem:param-counting (redukcja 3→2 — domknięcie LICZBOWE); dodatekQ CG-2 (LPA'/MC numeryczny)"
reuse: "dodatekQ CG-2 (tgp_erg_lpa_prime.py, 8/8 PASS; ρ₀*=0.03045, ν=0.749, K_IR/K_UV=1.000); op-Csigma-coarse-graining Phase 2 (bąbel) + Phase 3 (T=C_σσ_0²)"
anti_lakatos_note: "Pre-derywacja = oczekiwanie (κ_E=O(1)≠5/6 generycznie → FALSIFIED-hard), rachunek nadrzędny. Zakaz strojenia T do 5/6. GAP (continuum nie zbiega / unit-bridge nieokreślony) = pełnoprawny wynik. Systematyki raportowane jawnie z błędami."
tags: [scoping, C-sigma, lattice-MC, tensor-stiffness, kappa-E-pinning, numerical, multi-session, anti-Lakatos]
---

# SCOPING — liczbowe T=C_σσ_0² z lattice-MC (domknięcie sektora radiacyjnego, hard)

## §1 — Pytanie i obiekt (precyzyjnie)
op-Csigma-coarse-graining (PARTIAL) zredukował sektor radiacyjny do **jednego** fizycznego parametru
$T=C_\sigma\sigma_0^2$ (tensor stiffness), z $\kappa_E=8\pi G_0 T/c^3$; survival ⟺ κ_E=5/6 (miara zero,
niechroniona); naturalna κ_E=1 → 7/6 FALSIFIED. **Jedyny residual GAP: liczbowa wartość $T$.**

> **Czy $T=C_\sigma\sigma_0^2$ jest wyznaczalne liczbowo z kierunkowego bąbla
> $\Pi_{ab,cd}(p)=\langle O_{ab}(p)O_{cd}(-p)\rangle_{\rm c}$, $O_{ab}=(\partial_a\hat s\,\partial_b\hat s)_{\rm TF}$,
> na sieci 3D Ising ($d=3,n=1,\mathbb Z_2$) — a złożone w $\kappa_E=8\pi G_0T/c^3$ daje $5/6$ (SURVIVE)
> czy $\neq5/6$ (FALSIFIED)?**

**Obiekt LOCKED (z parent):** $T=C_\sigma\sigma_0^2$ — JEDEN parametr (redundancja przeskalowania
$\sigma\to\lambda\sigma$, Phase 3, sympy 3/3). NIE liczyć $C_\sigma$ i $\sigma_0$ osobno (gauge). $C_\sigma$
= współczynnik $p^2$ w $\Pi$ (Phase 2: bąbel analityczny $\Pi(p)=\tfrac1{8\pi m}-\tfrac{p^2}{96\pi m^3}$ —
wzorzec do reprodukcji liczbowej z realnym, nie-Gaussowskim substratem).

## §2 — Trasa wyprowadzenia (fazy; klasa dodatekQ CG-2 numeryczny)

**Phase 1 — Setup sieci + operator + obserwable.** Sieć 3D Ising/GL ($H_\Gamma$ eq:B-H), parametry z
dodatekQ CG-2 ($\rho_0^*$, $\nu=0.749$). Definicja sieciowa $O_{ab}=(\partial_a\hat s\,\partial_b\hat s)_{\rm TF}$
(różnice skończone; rozstrzygnąć notację skalar-vs-wektor ŝ, flaga Phase 2 §1). Obserwable: $\Pi_{ab,cd}(p)$
w przestrzeni pędów (FFT konfiguracji). Algorytm klastrowy (Wolff/Swendsen-Wang) — R-critical-slowing.

**Phase 2 — Pomiar bąbla + ekstrakcja C_σ (RDZEŃ numeryczny).** MC na rosnących sieciach $L$; pomiar
$\Pi(p)$ przy małych $p$; dopasowanie $\Pi(p)=\Pi(0)-A p^2+\dots$ → $A$ → $C_\sigma$ (po renormalizacji
operatora złożonego). Ekstrapolacja continuum ($a_{\rm sub}\to0$, $b\to\infty$) i skończony rozmiar (FSS).
Walidacja względem analitycznego bąbla Gaussowskiego (Phase 2 parent) w reżimie słabego sprzężenia.

**Phase 3 — Unit-bridge → T w jednostkach c³/G_0 + złożenie κ_E.** Most jednostkowy: $G_0\sim J\mu$
(rem:param-counting) + samospójność $a_\Gamma\cdot\Phi_0=1$ (dodatekQ Q.4). $\sigma_0$ ze skali korelatora
(~Φ_0). Złożyć $T=C_\sigma\sigma_0^2$ w jednostkach $c^3/G_0$ → $\kappa_E=8\pi G_0T/c^3$ z błędem
statystycznym + systematycznym (continuum, unit-bridge, renormalizacja).

**Phase FINAL — κ_E vs 5/6 → werdykt liczbowy.** Porównanie value-blind (próg 5/6 pre-LOCKED Phase 0):
$\kappa_E=5/6\pm\sigma_{\rm tot}$ → SURVIVE; $\neq$ (poza pasmem błędu) → **FALSIFIED hard**; szerokie pasmo
obejmujące 5/6 i 1 → PARTIAL (sektor pozostaje UNDERDETERMINED-fine-tuned, ale z liczbowym przedziałem).

## §3 — Falsyfikatory (szkic do LOCK w Phase 0)

| ID | Test | Reguła |
|---|---|---|
| F-LMC-A | Czy $\Pi(p)$ ma mierzalny, dodatni współczynnik $p^2$ (C_σ>0 skończone), zbieżny w continuum? | TAK → przejdź; NIE (brak skończonej sztywności / niezbieżność) → GAP strukturalny |
| F-LMC-B | Ekstrakcja $T=C_\sigma\sigma_0^2$ w jednostkach $c^3/G_0$ (unit-bridge) | wartość ± błąd / PARTIAL (systematyka unit-bridge dominuje) / GAP |
| F-LMC-C | $\kappa_E=8\pi G_0T/c^3$ vs 5/6 (próg LOCKED) | $=5/6\pm\sigma_{\rm tot}$ SURVIVE; ≠ FALSIFIED; pasmo szerokie → PARTIAL |
| F-LMC-D | Agregat | DERIVED+κ_E≠5/6 → **FALSIFIED hard**; +5/6 → SURVIVE; GAP/PARTIAL → UNDERDETERMINED-fine-tuned (liczbowy przedział) |

## §4 — Reuse machinery (LOCKED)
- dodatekQ CG-2: `tgp_erg_lpa_prime.py` (8/8 PASS; ρ₀*=0.03045, ν=0.749, K_IR/K_UV=1.000) — wzorzec numeryczny + parametry punktu stałego.
- op-Csigma-coarse-graining Phase 2 (bąbel analityczny — walidacja) + Phase 3 (T=C_σσ_0², unit-bridge struktura).
- dodatekQ Q.4 (a_Γ·Φ_0=1 samospójność, `tgp_cg5_phi0_self_consistency.py` 8/8 PASS) — most jednostkowy.

## §5 — Forbidden moves (szkic)
1. **Zakaz strojenia $T$ do κ_E=5/6** (value-blind; survival musi WYNIKNĄĆ z MC, nie być dobrane).
2. Zakaz liczenia $C_\sigma$, $\sigma_0$ osobno jako fizycznych — tylko niezmienniczy $T=C_\sigma\sigma_0^2$ (redundancja, parent Phase 3).
3. Zakaz rewizji parent/PR-025/rdzeń LOCKED.
4. **GAP/PARTIAL jawnie** jeśli continuum nie zbiega lub unit-bridge nieokreślony — zakaz fabrykowania wartości (anti-Lakatos).
5. **Systematyki raportowane z błędami** (continuum, FSS, renormalizacja operatora, unit-bridge) — zakaz ukrywania niepewności.
6. Budżet nowych stałych 0 (T istnieje; {J,a_sub,μ,...} = parametry substratu).

## §6 — Risk register
- **R-continuum (HIGH):** ekstrapolacja continuum bąbla operatora ZŁOŻONEGO (renormalizacja $O_{ab}$) — realny GAP, jeśli mixing/divergencje. dodatekQ CG-3/CG-4 (continuum) są OTWARTE analitycznie.
- **R-unit-bridge (HIGH):** połączenie jednostek substratu z $c^3/G_0$ wymaga $G_0\sim J\mu$ + $a_\Gamma\Phi_0=1$; systematyka może zdominować i dać tylko PARTIAL (przedział).
- **R-critical-slowing (MED):** blisko krytyczności (ξ_corr duże, gdzie bąbel ~1/m³ wzmocniony) — wymóg algorytmów klastrowych.
- **R-tensor-projection (MED):** definicja sieciowa $O_{ab}$ bezśladowego + notacja ŝ skalar-vs-wektor (flaga Phase 2 §1) — wpływa na prefaktor tensorowy.

## §7 — Oczekiwane wyniki (INFORMATIONAL; rachunek nadrzędny)
Najbardziej prawdopodobne: **DERIVED → κ_E=O(1)≠5/6 → FALSIFIED hard** (zgodnie z lean parent: brak symetrii
chroniącej 5/6). Drugi: **PARTIAL** (unit-bridge/continuum systematyka daje przedział obejmujący 5/6 i 1 —
sektor pozostaje UNDERDETERMINED-fine-tuned, ale liczbowy). Trzeci (mało prawdopodobny): κ_E=5/6±σ SURVIVE
(spisek). Czwarty: **GAP** (continuum operatora złożonego niezbieżne).

## §8 — Rejestracja
**`op-Csigma-lattice-MC`** — REGISTERED-QUEUED (NOT activated; własny Phase 0 + „działaj"). Wielosesyjny
(klasa dodatekQ CG-2 numeryczny). Rozstrzyga OSTATNI residual GAP sektora radiacyjnego (liczbowe $T$ →
twardy werdykt FALSIFIED / SURVIVE). [[../research/op-Csigma-lattice-MC-2026-06-14/]] (README).
