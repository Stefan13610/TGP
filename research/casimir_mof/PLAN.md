# Efekt Casimira w nietypowych geometriach (MOF, mikropory)

**Data startu:** 2026-04-20
**Status:** plan / open
**Kategoria:** TGP applications → próżnia substratowa → nanogeometria

## 1. Problem fizyczny

Klasyczny efekt Casimira (Casimir 1948) daje siłę przyciągającą między
dwiema doskonale przewodzącymi płytami:
$$F/A = -\frac{\pi^2}{240}\,\frac{\hbar c}{d^4}.$$

W XXI wieku rozszerzenie tej teorii do realnych materiałów (Lifshitz 1956,
Dzyaloshinskii–Pitaevskii 1961) i złożonych geometrii (sfery, klin, Casimir–Polder)
jest żmudne. W szczególności dla **mikroporowatych struktur** (MOF-5, UiO-66, ZIF-8)
z porami rzędu 1–5 nm brak jest zwięzłego analitycznego ujęcia:

- ciśnienie osmotyczne wewnątrz porów,
- kontrybucja próżniowa do izoterm adsorpcji,
- stabilność klatki pod próżnią.

Typowe eksperymenty: Mohideen/Roy 1998 (płyta-sfera), Decca 2003 (torsion),
dla MOF: pomiary mechaniczne (bulk modulus) + DFT. Rozbieżności z teorią na
poziomie 10-30% sugerują brak uniwersalnego modelu.

## 2. Dlaczego TGP

Aksjomat A-I TGP (absolutny stan odniesienia $N_0$) daje **naturalny baseline**
energii próżni: nie "zero arbitralne", lecz "energia substratu bez materii".
Fluktuacje substratu $\delta\Phi$ mają skończoną gęstość energii
$$\epsilon_{\text{vac}} = \frac{1}{2}\big\langle(\partial_t\delta\Phi)^2
  + c_0^2(\nabla\delta\Phi)^2\big\rangle$$
a ograniczenia przez geometrię porów kwantyzują mody fluktuacji, co daje
siłę Casimira **bezpośrednio z substratu**, bez "renormalizacji przez
odejmowanie nieskończoności".

Kluczowa obserwacja: w TGP **przestrzeń jest generowana przez materię** (A-IV),
więc w pustym porze $\Phi \to \PhiZero$ (substrat niezaburzony), a ściany
MOF (atomy metalu + linkera) tworzą gradient $\nabla\Phi$, który deformuje
widmo modów.

## 3. Cele badawcze

### C1 — Mody substratu w pore-geometrii

Znaleźć numerycznie własne mody operatora Laplace'a (liniowy reżim TGP)
w geometrii:
- 2D kanał (klasyczny Casimir),
- kulista wnęka (MOF-5 cage, d ≈ 0.8 nm),
- kanały cylindryczne (SBA-15, d = 2-10 nm),
- 3D MOF siatka pełna (periodic boundary).

### C2 — Ciśnienie Casimira w MOF

Obliczyć
$$P_{\text{Cas}}(d) = -\frac{\partial E_{\text{zpf}}}{\partial V}$$
dla realistycznych wartości $d$ i porównać z:
- klasycznym Lifshitz: $P \sim -\hbar c/d^4$,
- osmotycznym równaniem stanu (gaz w porach — odjąć, żeby zobaczyć czysty Casimir),
- bulk modulus MOF (Yot 2012 dla MOF-5: K ≈ 15 GPa).

### C3 — Nowa skala TGP dla Casimira w gradient-substratie

Pytanie: czy gradient $\Phi$ (gęstość "ścian" porów) wprowadza dodatkową skalę
$\ell_{\Phi} = (\hbar c/\beta q\PhiZero^2)^{1/2}$, która modyfikuje $1/d^4$?
Wyprowadzić poprawkę:
$$P_{\text{Cas}}^{\text{TGP}} = -\frac{\pi^2\hbar c}{240 d^4}\big(1 + c_1 (\ell_\Phi/d)^p\big)$$
z $c_1, p$ przewidzianymi z rdzenia TGP (a nie fitowanych).

### C4 — Test eksperymentalny

Znajdź MOF z kontrolowalnym $d$ (gościnne molekuły zmieniają efektywny
rozmiar pora) i przewiduj trend $P(d)$ niezgodny z czystym Lifshitz, ale
zgodny z TGP-correction.

## 4. Plan numeryczny

- **ps01_casimir_parallel_plates.py** — reprodukcja klasyki; weryfikacja metody.
- **ps02_casimir_cylindrical_pore.py** — mode-sum dla cylindra (SBA-15 like).
- **ps03_casimir_MOF_cage.py** — kulista wnęka + Bloch-periodic siatka.
- **ps04_osmotic_vs_casimir.py** — rozłożenie równania stanu dla argonu w MOF-5.
- **ps05_gradient_phi_correction.py** — poprawka TGP od gradientu substratu.

## 5. Literatura startowa

- Casimir 1948 (oryginał), Lifshitz 1956, Milton *The Casimir Effect* (2001)
- Mohideen & Roy, PRL 81, 4549 (1998) — precision measurement
- Bressi et al., PRL 88, 041804 (2002) — dwa walce
- MOF-5 mechanical: Yot et al., Dalton Trans. 41, 3813 (2012)
- Review: Klimchitskaya–Mostepanenko, Rev. Mod. Phys. 81, 1827 (2009)

## 6. Relacje z innymi sektorami TGP

- **Core paper Eq. (field)** dla linearized regimes — bezpośrednia baza.
- **continuum_limit**: kwantyzacja modów substratu na siatce → continuum
  → Casimir spektrum.
- **mass_scaling_k4**: zależność "siły podstawowej" od skali TGP.

## 7. Falsyfikowalność

- TGP przewiduje konkretny znak i wielkość $c_1$ w poprawce gradient-term.
- Obserwacja: $P(d)$ w MOF-5 (1.1 nm) i MOF-177 (1.5 nm) powinna odchylać się
  od $1/d^4$ o mierzalny procent.
- Jeśli odchylenie się nie pojawi → korekta TGP jest podrzędna, ale to też
  informacja (ograniczenie $\ell_\Phi < 0.1$ nm na przykład).

## 8. Otwarte pytania

1. Jak prawidłowo wprowadzić warunki brzegowe dla $\Phi$ na ścianie MOF
   (Dirichlet? Robin? "soft-wall" z profilu $\Phi$ atomu metalu?)
2. Czy Casimir-Polder dla pojedynczej molekuły gościa w porze dostaje
   poprawkę TGP inną niż ścianowa?
3. Czy istnieje analog "Boyer repulsive Casimir" w specjalnej geometrii
   MOF, wynikający z chiral Z_2 symmetry substratu?

## 9. Link do rdzenia TGP

Core paper [[papers/core/tgp_core.pdf]] § "Effective field theory" + § "Emergent metric".
Casimir to pierwszy test TGP **na próżni** (nie na materii), uzupełniający
grawitacyjny reżim (PPN) i kondensowaną materię (SC).
