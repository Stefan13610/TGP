# R2: Ciągłe przejście substrat → pole (CG-1/3/4)

## Problem

Trzy otwarte twierdzenia blokują claim "TGP wyprowadzone z pierwszych zasad":

| Twierdzenie | Opis | Status |
|-------------|------|--------|
| **CG-1** | Istnienie i jednoznaczność punktu stałego S* operatora blokowania (kontrakcja Banacha) | OTWARTE |
| **CG-3** | Zbieżność Φ_B → Φ w H¹ (homogenizacja, twierdzenie de Giorgi–Nash–Moser) | OTWARTE |
| **CG-4** | Identyfikacja K_hom = K_TGP | OTWARTE |

Bez nich: TGP jest dobrze motywowaną **teorią efektywną**.
Z nimi: TGP jest **wyprowadzone z pierwszych zasad**.

## Obecny status

- Słabe twierdzenie α=2: ZAMKNIĘTE (Lemma A1–A5, a1_alpha2_frg_synthesis.py: 7/7 PASS)
- CG-2 (numeryczny): K_IR/K_UV = 1.000 (FRG LPA', 8/8 PASS)
- Silne twierdzenie: OTWARTE — czysta matematyka

## Plan ataku

### CG-1: Kontrakcja Banacha operatora blokowania

**Cel:** Dowieść że operator Kadanoffa T: F → F na przestrzeni hamiltonianów
ma jednoznaczny punkt stały S* i T^n(H₀) → S* dla dowolnego H₀ w klasie.

**Narzędzia:**
- Teoria ERG (Polchinski, Wetterinck)
- Bauerschmidt, Brydges, Slade — "Renormalisation Group Analysis" (2019)
- Kontraktywność w normie operatorowej na algebrze operatorów

**Podproblemy:**
1. Zdefiniować przestrzeń Banacha F hamiltonianów na Γ
2. Wykazać kontraktywność: ||T(H₁) - T(H₂)|| ≤ q·||H₁ - H₂||, q < 1
3. Zidentyfikować punkt stały S* z hamiltonianem TGP

### CG-3: Homogenizacja Φ_B → Φ

**Cel:** Dowieść zbieżność pola zgrubnego Φ_B do ciągłego Φ w H¹(R³).

**Narzędzia:**
- Twierdzenie de Giorgi–Nash–Moser (regularność)
- Teoria homogenizacji (Jikov, Kozlov, Oleinik)
- Γ-zbieżność funkcjonałów energii

**Podproblemy:**
1. Sprawdzić warunki: eliptyczność, ograniczoność współczynników
2. Dowieść zbieżność E_B[Φ_B] → E[Φ] w sensie Γ
3. Wykazać H¹ regularność granicy

### CG-4: Identyfikacja K_hom = K_TGP

**Cel:** Dowieść że K(Φ) = Φ^α z α=2 jest jedynym K spójnym z CG-1 i CG-3.

**Zależy od:** CG-1 + CG-3

## Kryterium zamknięcia

Formalne dowody trzech twierdzeń (CG-1, CG-3, CG-4) na poziomie
publikowalnym w J. Math. Phys. / Comm. Math. Phys.

## Pliki do scalenia z rdzeniem

- Rozszerzenie `dodatekQ_coarse_graining_formal.tex`
- Nowy `dodatek_CG_proof.tex` z pełnymi dowodami

## Referencje rdzenia

- `dodatekQ_coarse_graining_formal.tex`, linie 170–193
- `scripts/a1_alpha2_frg_synthesis.py` (słabe twierdzenie)
- `PLAN_ROZWOJU_v3.md`, linie 27–50

## Uwaga

To jest **czysta matematyka**. Może być publikowalne niezależnie od fizyki TGP
jako wynik w teorii renormalizacji / homogenizacji.

## Status

- [ ] CG-1: Definicja przestrzeni Banacha F
- [ ] CG-1: Dowód kontraktywności
- [ ] CG-1: Identyfikacja S*
- [ ] CG-3: Sprawdzenie warunków homogenizacji
- [ ] CG-3: Dowód Γ-zbieżności
- [ ] CG-4: Identyfikacja K z CG-1+CG-3
- [ ] Redakcja artykułu matematycznego
