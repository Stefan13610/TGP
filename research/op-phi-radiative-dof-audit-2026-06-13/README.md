---
title: "op-phi-radiative-dof-audit-2026-06-13"
folder_status: closed-resolved
claim_status: "HONEST_NEGATIVE ⟹ PR-025 EXHAUSTIVE-OVER-LIVE (pending user ratification)"
date_opened: 2026-06-13
date_closed: 2026-06-13
parent_pr: PR-025
---

# op-phi-radiative-dof-audit — czy Φ może być zmienną pomocniczą (non-radiative)?

**Pytanie (sformułowane przez autora TGP):** brakujący warunek sektora
radiacyjnego = dowód, że Φ jest tylko zmienną pomocniczą dla σ_ab, a nie
dodatkowym radiacyjnym stopniem swobody. Rozstrzyga P_total = κ_E·P_GR
vs κ_E·P_GR + (1/6)P_GR.

**Wynik (13/13 sympy PASS, fast-kill 1 sesja, 0 danych obserwacyjnych):**
- **F-AUX-A NEGATIVE:** Hessian = K₁ > 0, zero więzów Diraca, DOF_Φ = 1; metoda zwalidowana na EM (A₀ wykryte jako auxiliary).
- **F-AUX-B NEGATIVE:** Lorentz + lock Newtona ⟹ F(s)=1/s ⟹ biegun radiacyjny nieusuwalny; zabicie bieguna = zniszczenie statyki.
- **F-AUX-C NEGATIVE:** Z₂ dyskretna (brak generatora), U(1) działa tylko na fazę — brak więzu na mod oddechowy.
- **F-AUX-D = HONEST_NEGATIVE** (werdykt wyliczony z flag).

**Konsekwencja:** kanał skalarny (1/6)P_GR nieusuwalny w LIVE ⟹
**PR-025 upgrade: EXHAUSTIVE-OVER-LIVE.**

Pliki: [[Phase0_balance.md]] → [[Phase1_sympy.py]]/[[Phase1_sympy.txt]] →
[[Phase1_results.md]] → [[Phase_FINAL_close.md]].
Scoping: [[../../meta/SCOPING_op-phi-radiative-dof-audit_2026-06-13.md]].
