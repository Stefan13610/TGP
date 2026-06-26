---
title: "Phase 2 — op-Csigma-coarse-graining: RDZEŃ. Metoda = propagator KOMPOZYTU (bąbel 2-punktowy) bo σ_ab=⟨ŝŝ⟩_TF jest operatorem złożonym. WYNIK: C_σ>0 z dodatniego, skończonego członu p² bąbla (sztywność przestrzenna emerguje z H_Γ); skaling C_σ~J·a_sub^p/(m·a_sub)^k DERIVED, prefaktor O(1)=GAP (anti-Lakatos). κ_E=8πG_0C_σσ_0²/c³ złożone; BRAK symetrii lockującej (det J≠0)→κ_E generyczne O(1), 5/6=miara zero. WERDYKT: PARTIAL, lean FALSIFIED; sektor UNDERDETERMINED-fine-tuned (węższy: κ_E O(1)-bounded)."
type: phase_result
status: PHASE2_COMPLETE (RDZEŃ — PARTIAL)
phase: 2
cycle: op-Csigma-coarse-graining-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'działaj z fazą 2'"
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
sympy: "PASS: bąbel 3D Π(0)=1/(8πm), coeff p²=−1/(96πm³) (calki beta-funkcyjne EXACT); sztywność/masa=1/(12m²); κ_E=8πG_0C_σσ_0²/c³; det J=−ξ_eff/C_σ≠0"
flag_F_CG_A: "PASS-CONDITIONAL (sztywność przestrzenna z H_Γ via bąbel; człon czasowy/c_0 dziedziczony z metryki emergentnej rem:cGW-tensor — nie de novo)"
flag_F_CG_B: "PARTIAL (metoda+skaling+znak DERIVED; prefaktor O(1)+wykładniki = GAP, schemat-zależne)"
flag_F_CG_C: "PARTIAL (κ_E=8πG_0C_σσ_0²/c³ złożone; σ_0 vs C_σ = Phase 3)"
flag_F_CG_D: "PARTIAL-lean-FALSIFIED (5/6 niechronione, miara zero; naturalna O(1)≠5/6)"
flag_F_CG_E: "UNDERDETERMINED-fine-tuned (WĘŻSZY: κ_E O(1)-bounded) ; lean FALSIFIED"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — RDZEŃ: ekstrakcja $C_\sigma$ z coarse-grainingu (PARTIAL)

## §0 — Werdykt fazy w skrócie
| Falsyfikator | Wynik |
|---|---|
| **F-CG-A** (Lorentz/kinetyka) | **PASS-CONDITIONAL** — bąbel daje dodatni skończony człon $p^2$ ⇒ sztywność **przestrzenna** $C_\sigma>0$ emerguje z $H_\Gamma$; człon **czasowy** $(\partial_t\sigma)^2$ i $c_0$ **dziedziczone** z metryki emergentnej (rem:cGW-tensor), nie wyprowadzone de novo |
| **F-CG-B** (forma $C_\sigma$) | **PARTIAL** — metoda (bąbel kompozytu) + skaling $C_\sigma\sim c_{\rm pref}J\,a_{\rm sub}^p/(m a_{\rm sub})^k$ + znak $>0$: **DERIVED**; prefaktor O(1) + wykładniki: **GAP** (schemat-/normalizacja-zależne) |
| **F-CG-C** (σ_0, złożenie) | **PARTIAL** — $\kappa_E=8\pi G_0C_\sigma\sigma_0^2/c^3$ złożone; rozdz. $C_\sigma$ vs $\sigma_0$ = Phase 3 |
| **F-CG-D** (κ_E vs 5/6) | **PARTIAL-lean-FALSIFIED** — brak symetrii lockującej (det J≠0) ⇒ κ_E generyczne O(1); **5/6 = miara zero (niechronione)** |
| **F-CG-E** (agregat) | **UNDERDETERMINED-fine-tuned (węższy)**; lean **FALSIFIED** |

## §1 — Metoda: kinetyka operatora ZŁOŻONEGO = propagator kompozytu
Phase 1 ustaliło: $\sigma_{ab}=\langle\hat s_i\hat s_j\rangle_{\rm TF}$ jest **operatorem złożonym** (bilinowym), nie polem fundamentalnym. Stąd jego człon kinetyczny **nie jest postulowany** — emerguje jako **odwrotność dwupunktowego korelatora kompozytu** (bąbla):
$$\Pi_{ab,cd}(p)=\big\langle\,O_{ab}(p)\,O_{cd}(-p)\,\big\rangle_{\rm c},\qquad O_{ab}=(\partial_a\hat s\,\partial_b\hat s)_{\rm TF},$$
a efektywne działanie σ to $S_\sigma=\tfrac12\int\sigma\,[g^{-1}-\Pi(p)]\,\sigma$ (Hubbard–Stratonovich na członie wiązania). Współczynnik przy $p^2$ w $\Pi$ daje $C_\sigma$; człon $p^0$ daje $m_\sigma^2$. **To jest standardowa, niefabrykowana procedura** (kinetyka kondensatu z bąbla), tu zastosowana do anizotropowego sektora substratu.

> **Nota o notacji źródła:** rdzeń miesza $\hat s_i\in\mathbb R$ (skalar, eq:B-H) z $\hat s_i^a$ (komponenty, fbox ssec:tensor-substrate). Przyjmuję interpretację spójną ze „spin-2, 5 komponentów": $K_{ab}$ z kierunkowych przesunięć skalara $\hat s$, którego część anizotropowa to $\langle\partial_a\hat s\,\partial_b\hat s\rangle_{\rm TF}$ (najniższy rank-2 bilinowy w pochodnych skalara). Alternatywna lektura (ŝ wektorowe) zmienia prefaktor tensorowy, **nie** zmienia werdyktu (skaling+brak-locka identyczne). Flaga do ujednolicenia w rdzeniu.

## §2 — Jądro rachunku (bąbel 3D, sympy EXACT)
Substrat: $d=3$, $n=1$, $\mathbb Z_2$ (dodatekQ). Propagator $G(q)=1/(q^2+m^2)$, $m=1/\xi_{\rm corr}$. Bąbel:
$$\Pi(p)=\int\frac{d^3q}{(2\pi)^3}\frac{1}{(q^2+m^2)\,((p+q)^2+m^2)}.$$
Rozwinięcie małego $p$ (zweryfikowane **całkami beta-funkcyjnymi**, `Phase2_sympy.txt`):
$$\boxed{\ \Pi(0)=\frac{1}{8\pi m},\qquad \Pi(p)=\Pi(0)-\frac{p^2}{96\pi m^3}+O(p^4)\ }$$
$$\frac{-[\text{coeff }p^2]}{\Pi(0)}=\frac{1}{12\,m^2}\quad(\text{sztywność/masa}).$$
**Konsekwencja:** bąbel jest **skończony i analityczny** w $p=0$ ⇒ emergentny człon $p^2$ **istnieje i jest dodatni**, z $C_\sigma\sim(\text{sprzężenie})^2\cdot\frac{1}{96\pi m^3}$. To **realizuje przestrzenną część F-CG-A**: H_Γ generuje dodatnią sztywność σ_ab.

## §3 — F-CG-A: Lorentz (PASS-CONDITIONAL, uczciwa granica)
- **Część przestrzenna $(\nabla\sigma)^2$:** DERIVED z bąbla (§2), $C_\sigma>0\Leftarrow J>0$. ✓
- **Część czasowa $(\partial_t\sigma)^2/c_0^2$:** $H_\Gamma$ jest **statyczny** — sam nie zawiera dynamiki czasowej. Człon czasowy i prędkość $c_0$ są **dziedziczone z metryki emergentnej** (rem:cGW-tensor: σ_ab propaguje z $c_0$ jako perturbacja NA geometrii ustanowionej przez Φ), **dokładnie jak w sektorze skalarnym** ($\Box\Phi$ z tego samego mechanizmu).
- **Dyspozycja:** **PASS-CONDITIONAL** — kinetyka Lorentza emerguje, ale przez **złożenie** (sztywność przestrzenna z H_Γ ⊕ struktura czasowa z emergentnej metryki). Lorentz **nie** jest wyprowadzony de novo z samego H_Γ. Nie deklaruję mocniejszego PASS niż to, na co pozwala rachunek. Brak sprzeczności z GW170817 (rem:cGW-tensor LOCKED).

## §4 — F-CG-B: forma $C_\sigma$ (PARTIAL — skaling DERIVED, prefaktor GAP)
Analiza wymiarowa + struktura bąbla (3D, $1/m^3$):
$$C_\sigma\;\sim\;c_{\rm pref}\;J\;a_{\rm sub}^{\,p}\;\Big(\tfrac{\xi_{\rm corr}}{a_{\rm sub}}\Big)^{k},\qquad C_\sigma>0.$$
- **DERIVED:** (i) metoda (bąbel kompozytu), (ii) znak $C_\sigma>0$, (iii) skaling — w tym **wzmocnienie ~ $\xi_{\rm corr}^3$** przy zbliżaniu do krytyczności ($m\to0$, z $1/m^3$ bąbla). To strukturalne, nie liczbowe.
- **GAP (residual):** prefaktor $c_{\rm pref}=O(1)$ oraz wykładniki $(p,k)$ zależą od: precyzyjnej dymensji operatora $O_{ab}$, normalizacji pola HS, schematu regularyzacji sieci, projekcji tensorowej. **Nie podaję jednej liczby** — to byłoby fabrykowanie (lekcja #30; §5 forbidden moves Phase 0). Zgodne z PARTIAL wg LOCKED reguły F-CG-B.

## §5 — F-CG-C/D: złożenie $\kappa_E$ i test 5/6 (PARTIAL-lean-FALSIFIED)
Z $h^{TT}_{ij}=2\sigma_{ij}/\sigma_0$ (eq:h-TT-from-sigma) podstawienie $\sigma=\sigma_0 h/2$ w $\tfrac{C_\sigma}{2}(\partial\sigma)^2$ daje sztywność tensorową TGP $\tfrac{C_\sigma\sigma_0^2}{8}(\partial h)^2$; grawiton GR ma $\tfrac{c^3}{64\pi G_0}(\partial h)^2$. Stosunek (= stosunek strumieni Isaacsona):
$$\boxed{\ \kappa_E=\frac{C_\sigma\sigma_0^2/8}{c^3/(64\pi G_0)}=\frac{8\pi G_0\,C_\sigma\,\sigma_0^2}{c^3}\ }$$
**Test value-blind (próg LOCKED Phase 0):**
- survival ⟺ $\kappa_E=5/6$ EXACT (miara zero);
- „naturalna" $\kappa_E=1$ ⟹ total $7/6$ = gałąź B PR-025 (2646σ FALSIFIED).

**Dlaczego PARTIAL, a nie liczbowy werdykt:** rozstrzygnięcie $=5/6$ wymaga prefaktora $C_\sigma$ (§4 GAP) **oraz** $\sigma_0$ (Phase 3). Skaling daje $\kappa_E=O(1)$, ale **nie odróżnia 5/6 od 1** bez prefaktora. Survival ani potwierdzone ani wykluczone liczbowo.

## §6 — Wynik strukturalny (nowy vs parent): BRAK symetrii lockującej $\kappa_E$
Najważniejszy nieliczbowy rezultat Phase 2:
- W **GR** jedno $G$ lockuje amplitudę↔strumień (1 parametr) ⇒ tensorowa relacja jest wymuszona.
- W **TGP**: amplituda zlockowana ($\xi_{\rm eff}=4\pi G_0\sigma_0\Phi_0$, thm:amplitude-matching), ale strumień $\propto C_\sigma\sigma_0^2$ **niezależny**: $\det J(\text{amp,flux}\,|\,C_\sigma,\sigma_0)=-\xi_{\rm eff}/C_\sigma\neq0$ (sympy). **2 parametry tam, gdzie GR ma 1.**
- $G_0\sim J\mu$ (sektor skalarny, $\langle\hat s^2\rangle$) i $C_\sigma\sim J a_{\rm sub}^p$ (sektor tensorowy, $\langle\hat s\hat s\rangle$) to **różne projekcje tego samego ŝ** — **brak tożsamości Warda** lockującej $C_\sigma\sigma_0^2$ do sztywności grawitonu.

> **Wniosek:** $\kappa_E$ jest generycznie $O(1)$, ale **ani 1, ani 5/6 nie są wymuszone** żadną symetrią. $5/6$ pozostaje **miarą zero (niechronioną)**. To **zaostrza** wynik rodzica: tam $\kappa_E$ „swobodne nad continuum"; tu $\kappa_E$ **ograniczone do O(1)** (przez strukturę bąbla) — węższy przedział, dalej obejmujący 5/6 i 1, ale „naturalna" $O(1)\neq5/6$ ⇒ **lean FALSIFIED** wzmocniony.

## §7 — Adwersarialny self-check (R-conspiracy, uczciwie obie strony)
Czy istnieje argument za $\kappa_E=5/6$ (SURVIVE)? **Tak, niebłahy:** jeśli TGP poprawnie redukuje się do GR, **całkowity** strumień powinien = GR, więc tensor musiałby nieść dokładnie $1-1/6=5/6$ (skalar konforemny zabiera 1/6). To „emergent-GR-consistency". **Kontra:** $\det J\neq0$ pokazuje, że amplitude-matching **nie** wymusza flux-matching — emergent-GR w amplitudzie nie przenosi się na strumień. Zatem $5/6$ nie jest wymuszone, lecz byłoby **szczęśliwym zbiegiem** dwóch niezależnych sektorów (R-conspiracy). **Obie hipotezy (κ_E=1 FALSIFIED, κ_E=5/6 SURVIVE) klamrują wynik; skaling sam nie rozstrzyga** → uczciwie PARTIAL, lean (nie hard) FALSIFIED.

## §8 — Anti-Lakatos (checklist Phase 2)
✓ Bąbel zweryfikowany EXACT (sympy, całki beta-funkcyjne) — metoda niefabrykowana.
✓ **Prefaktor $C_\sigma$ NIE podany** — tylko skaling+znak; jawny GAP (§4). Zero strojenia do 5/6.
✓ F-CG-A nie zadeklarowane mocniej niż PASS-CONDITIONAL (część czasowa dziedziczona — flaga §3).
✓ Obie strony R-conspiracy przedstawione (§7) — brak framework-protection bias.
✓ Lean FALSIFIED oznaczony jako **strukturalny (nie-twardy)**, bo $\det J\neq0$ ⇒ κ_E formalnie wolne.
✓ Budżet nowych stałych 0; liczby poprzedników LOCKED.

## §9 — Handoff do Phase 3
**Phase 3:** normalizacja $\sigma_0$ z $h^{TT}=2\sigma/\sigma_0$ + rozstrzygnięcie, czy $C_\sigma$ i $\sigma_0$ to **jeden** parametr (rem:sigma-params: „$C_\sigma$ lub równoważnie $\sigma_0$") czy **dwa** — to może zredukować $\kappa_E$ do jednej kombinacji i zaostrzyć test. **Residual GAP** (prefaktor O(1) $C_\sigma$): docelowo **lattice-MC kierunkowego bąbla** $\langle O_{ab}O_{cd}\rangle$ na sieci 3D Ising — jedyna droga do liczbowego $\kappa_E$ i twardego werdyktu. Phase 3 wymaga „działaj".
