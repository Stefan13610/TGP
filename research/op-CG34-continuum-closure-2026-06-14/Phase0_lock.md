---
title: "Phase 0 — op-CG34-continuum-closure: LOCK kryteriów numerycznego domknięcia CG-3 (zbieżność Φ_B→Φ w H¹, homogenizacja) i CG-4 (identyfikacja K_hom=K_TGP). Cel: podnieść CG-3/CG-4 z [OTWARTY] do [ZAMKNIĘTY NUM] (standard CG-2/CG-5), naprawiając błędne wcześniejsze numeryki (H¹ norm=0, c*→0 artefakt ⟨|∇Φ|²⟩, CG-4 fit R²=0.33). Kryteria z PLAN_NUMERYCZNY_CG3_CG4 + lematy A1–A5 (dodatekQ2). Anti-Lakatos."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-CG34-continuum-closure-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'działaj z domknięciem CG-3/CG-4'"
target: "dodatekQ CG-3 [OTWARTY, AN] + CG-4 [OTWARTY, AN+NUM] → [ZAMKNIĘTY NUM]"
reuse: "dodatekQ2 lematy A1–A5; PLAN_NUMERYCZNY_CG3_CG4 (N1–N5); CG-2 (K_IR=ρ, ν=0.749); CG-5 (a_Γ·Φ_0=1); substrate_mc_cg3.py (SUPERSEDED — bug K₁)"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-CG34-continuum-closure)

> **Cel:** domknąć **numerycznie** CG-3 (homogenizacja: $\Phi_B\to\Phi$ w $H^1$) i CG-4 (identyfikacja
> $K_{\rm hom}=K_{\rm TGP}$), podnosząc je z **[OTWARTY]** do **[ZAMKNIĘTY NUM]** (standard CG-2/CG-5).
> Wszystkie kryteria PASS/FAIL zalockowane TUTAJ, przed pierwszą liczbą (anti-Lakatos).

## §1 — Co jest otwarte (z rdzenia dodatekQ, tabela statusu)
| Krok | Treść | Typ | Status (przed cyklem) |
|---|---|---|---|
| CG-3 | Zbieżność $\Phi_B\to\Phi$ w $H^1$ (homogenizacja) | [AN] | **[OTWARTY]** |
| CG-4 | Identyfikacja $K_{\rm hom}=K_{\rm TGP}$ | [AN+NUM] | **[OTWARTY]** |

Lematy A1–A5 (dodatekQ2) dają **słabą** wersję (kompaktowość/podsekwencja zamiast pełnej zbieżności;
częściowa identyfikacja). Pełne CG-3/CG-4 wymagają zamknięcia numerycznego (jak CG-2/CG-5).

## §2 — DIAGNOZA błędów wcześniejszych numeryk (do naprawienia)
1. **`cg_strong_numerical.py` CG-3:** $\|\Delta\Phi\|_{L^2}=\|\Delta\Phi\|_{H^1}=0.000000$ — **zepsuty pomiar** (zera). Do naprawy: poprawny H¹ na wspólnej siatce.
2. **`substrate_mc_cg3.py` c*:** $c_*(L{=}16)=6.3\times10^{-4}\to c_*(L{=}32)=2.1\times10^{-5}$ — **spadek 30× z L** ⟹ pozorne $c_*\to0$ (groziłoby upadkiem A1). **Przyczyna:** $K_1$ utożsamione z $\langle|\nabla\Phi_B|^2\rangle$, które maleje gdy fluktuacje blokowe się uśredniają — **artefakt**, nie fizyka.
3. **`cg_strong_numerical.py` CG-4:** fit $K\sim a\Phi+b$ z $R^2=0.33$, intercept $-3.6\times10^7$ — **garbage** (1D, nieznormalizowany, $\xi$ jako proxy K).

**Poprawka rdzeniowa (LOCKED):** sztywność $K_1$ mierzymy z **statycznego structure factor**
$S_\Phi(k)=T/(2K_1\hat k^2+m^2)$ ⟹ $1/S_\Phi(k)$ liniowe w $\hat k^2$, nachylenie $\to K_1$. To jest
**scale-correct** (niezależne od rozmiaru bloku w continuum) — w przeciwieństwie do $\langle|\nabla\Phi|^2\rangle$.

## §3 — Kryteria domknięcia LOCKED (PASS/FAIL)

### CG-3 (homogenizacja, $H^1$)
| ID | Test | Reguła PASS |
|---|---|---|
| **C3-A** | Uniform bound $\sup_b\|\Phi_B\|_{L^2}<\infty$ (założenie A1(i)) | $\|\Phi_B\|_{L^2}$ ograniczone, stabilne w $b$ (rozrzut <20%) |
| **C3-B** | Zanik korelacji $\langle\Phi\Phi\rangle_c\sim e^{-r/\xi}$ (W4) | dopasowanie eksponencjalne, $\xi$ skończone |
| **C3-C** | **$H^1$ Cauchy**: $\|\Phi_B-\Phi_{2B}\|_{H^1}$ maleje z $b$ | monotoniczny spadek (naprawione zera); norma NIEZEROWA |
| **C3-D** | Rate vs A5: $\|\Phi_B-\Phi\|_{L^2}\le C_1(L_B/\xi)+C_2(a/L_B)^{1/2}$ | dopasowanie z $C_1,C_2>0$, zgodność jakościowa |

### CG-4 (identyfikacja $K_{\rm hom}=K_{\rm TGP}$)
| ID | Test | Reguła PASS |
|---|---|---|
| **C4-A** | $c_*=\min_\Phi K_1(\Phi)>0$ **continuum-stabilne** (structure factor) | $c_*>0$ i **NIE maleje** systematycznie z $L$ (rozwiązanie red-flag) |
| **C4-B** | $\alpha=2$: $K_1(\Phi)\cdot\Phi=\text{const}$ (lemat A3, $K_1\propto1/\Phi$) | $K_1\cdot\Phi$ stałe ±20% w skanie $T$ (różne $\Phi_0$) |
| **C4-C** | $\beta_{\rm eff}=\gamma_{\rm eff}$ (lemat A4, $U'(\Phi_0)=0$) | $\|\beta/\gamma-1\|<0.3$ lub argument strukturalny |
| **C4-D** | **Residuum operatora N5**: $\|R\|/\|\nabla^2\Phi\|\to0$ z $b$ | $\|R\|_{\rm rel}<0.2$ dla najlepszego $b$ ORAZ maleje z $b$ |

## §4 — Reguła agregatu LOCKED
- **CG-3 = ZAMKNIĘTY NUM** ⟺ C3-A..C3-C PASS (C3-D wspierające).
- **CG-4 = ZAMKNIĘTY NUM** ⟺ C4-A (c*>0 stabilne) ∧ C4-B (α=2) ∧ C4-D (residuum) PASS (C4-C strukturalny).
- Jeśli C4-A pokazuje $c_*\to0$ w continuum → **CG-4 = GAP** (A1 zagrożone — jawnie).
- Werdykt **wyliczony** z tej reguły w Phase FINAL, nie wybrany.

## §5 — Forbidden moves (anti-Lakatos)
1. **Zakaz utożsamiania $K_1$ z $\langle|\nabla\Phi|^2\rangle$** (błąd prior; tylko structure-factor stiffness).
2. **Zakaz strojenia** parametrów, by $\alpha=2$ czy $\beta=\gamma$ „wyszło" — mierzone, raportowane z błędem.
3. **GAP/PARTIAL jawnie** jeśli c* maleje w continuum lub residuum nie zbiega — zakaz fabrykowania ZAMKNIĘTY.
4. **Walidacja silnika MC OBOWIĄZKOWA** przed pomiarami (faza, ⟨φ²⟩, ξ) — inaczej pomiary nieważne.
5. **Zakaz rewizji rdzenia** bez autoryzacji (status update = rekomendacja, nie edycja core .tex w tym cyklu).
6. Budżet nowych stałych 0; reuse CG-2/CG-5/A1–A5.

## §6 — Reuse (LOCKED)
- Lematy A1–A5 (dodatekQ2) — struktura analityczna (kompaktowość, lokalność, α=2, β=γ).
- PLAN_NUMERYCZNY_CG3_CG4 (N1–N5) — specyfikacja etapów.
- CG-2 (`tgp_erg_lpa_prime.py`): $K_{\rm IR}(\rho)=\rho$, $\nu=0.749$, $\Phi_0=2\rho_0^*=0.0609$.
- CG-5: $a_\Gamma\Phi_0=1$.
- `substrate_mc_cg3.py` — **SUPERSEDED** (bug K₁); reuse tylko struktury (block-averaging, V_eff).

## §7 — Oczekiwanie (INFORMATIONAL, nie wiąże)
Po naprawie pomiaru K₁ (structure factor): $c_*>0$ powinno być **continuum-stabilne** (red-flag był
artefaktem), $\alpha=2$ i $\beta=\gamma$ potwierdzone, residuum malejące ⟹ **CG-3/CG-4 ZAMKNIĘTY NUM**.
Ryzyko: jeśli $c_*$ realnie maleje → GAP. Werdykt rozstrzygnie pomiar.

## §8 — Status
**🔒 PHASE 0 LOCKED.** Anti-Lakatos aktywny. Kryteria C3-A..D, C4-A..D, reguła agregatu, forbidden moves —
zalockowane. Przejście do Phase 1 (silnik MC + poprawny K₁).
