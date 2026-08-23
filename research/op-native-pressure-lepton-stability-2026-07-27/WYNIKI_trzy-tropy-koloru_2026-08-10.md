# Trzy tropy koloru — sprawdzone u źródła. Dwa MARTWE, jeden ODSŁANIA ROZWIDLENIE rdzenia

**Data:** 2026-08-10
**Typ:** KONTROLA SKŁADNIKÓW (u źródła + weryfikacja rachunkowa) — **wynik negatywny dla tropów 1 i 2**
**Cel programu:** znaleźć struktury niosące KOLOR; ładunek ma być pochodną koloru.
**Skrypty:** `COLOR_z3_sources_check.py`, `COLOR_rank3_nogo.py`

---

## 0. Wynik w trzech zdaniach

> **Trop 1 (korelator rangi 3) — MARTWY, i to podwójnie**: substrat TGP to skalar Isinga
> z **twardą ℤ₂**, więc korelatory nieparzyste w ŝ znikają; a moment kierunkowy rangi 3
> z korelatora wiązaniowego znika **tożsamościowo** z samej niezmienniczości translacyjnej.
> **Trop 3 — potwierdzony: ℤ₂ jest twarda**, i to ona zabija trop 1.
> **Trop 2 — ℤ₃ w rdzeniu to TRZY RÓŻNE grupy w trzech różnych rolach; BRAK POWIĄZANIA
> z ℤ₃ z abelianizacji 2T.** Przy okazji wyszło, że rdzeń ma **dwa niezgodne substraty**.

---

## 1. Trop 3 (sprawdzam pierwszy, bo rozstrzyga trop 1): ℤ₂ jest TWARDA

`axioms/substrat/dodatekB_substrat.tex` eq:B-H, lin. 51–58 — **dosłownie**:

```
H_Γ = Σ_i [ m₀²/2 · ŝ_i² + λ₀/4 · ŝ_i⁴ ]  −  J Σ_⟨ij⟩ ŝ_i ŝ_j
```
> „z symetrią **ℤ₂: ŝ_i → −ŝ_i**"

**Węzeł substratu niesie RZECZYWISTY SKALAR ŝ_i** (model typu Ising/φ⁴), nie wektor,
nie tensor, nie spinor. ℤ₂ nie jest ozdobnikiem — jest **nośna**:

| Gdzie | Cytat / rola |
|---|---|
| lin. 126 | „Ponieważ **ŝ_i² ≥ 0**, mamy Φ ≥ 0 — aksjomat ax:przestrzen zapewniony **automatycznie**" |
| lin. 228–232 | „H_Γ jest niezmienniczy względem ŝ→−ŝ. Pole Φ=⟨ŝ²⟩ jest **parzyste**, więc potencjał efektywny [nie ma] nieparzystych potęg ŝ" |
| lin. 38 | klasa uniwersalności **3D Ising** (η_3D Ising) |

⟹ **Usunięcie ℤ₂ rozwala fundament sektora metrycznego** (Φ ≥ 0 przestaje być automatyczne,
V_eff dostaje człony nieparzyste). **ℤ₂ jest twarda. Trop 3 potwierdzony.**

## 2. Trop 1: korelator rangi 3 — DWIE niezależne przeszkody

### (a) Przeszkoda ℤ₂ (wewnętrzna)

Rangę tensora buduje się z przesunięć; parzystość — z liczby czynników ŝ:

| Obserwabla | # ŝ | ℤ₂-parzysta? | # przesunięć |
|---|---|---|---|
| Φ = ⟨ŝ²⟩ | 2 | ✅ TAK | 0 |
| σ_ab = ⟨ŝ ŝ₊â⟩^TF | 2 | ✅ TAK | 1 |
| **⟨ŝ ŝ₊â ŝ₊b̂⟩** | **3** | ❌ **NIE → ZNIKA** | 2 |
| ⟨ŝ ŝ₊â ŝ₊b̂ ŝ₊ĉ⟩ | 4 | ✅ TAK | 3 |

W fazie złamanej (v=⟨ŝ⟩≠0) rozkład ŝ = v+δs daje: człon v³ (bezkierunkowy),
człon v·⟨δs δs⟩ (**redukowalny do σ_ab — nic nowego**), oraz ⟨δs δs δs⟩_c,
który wymaga **sprzężenia kubicznego λ₃ŝ³ — zabronionego przez ℤ₂**.
⟹ **żaden nowy niezależny parametr porządku nie powstaje.**

### (b) Przeszkoda inwersyjna — mocniejsza, nie wymaga nawet ℤ₂

Jedyny więz na korelator wiązaniowy to **niezmienniczość translacyjna**:
C(−â) = ⟨ŝ_i ŝ_{i−â}⟩ = ⟨ŝ_{i+â} ŝ_i⟩ = **C(â)**, exact.

`COLOR_rank3_nogo.py` — **dowolne anizotropowe** C spełniające C(−â)=C(â):

| sąsiedztwo | #wiązań | max\|M_ij^TF\| (rank 2) | max\|T_ijk\| (rank 3) |
|---|---|---|---|
| 6 | 6 | 0.580266 | **0.0** |
| 18 | 18 | 4.461614 | 2.8e−16 |
| 26 | 26 | 8.356380 | 3.1e−16 |
| 124 | 124 | 59.593014 | 9.5e−15 |

**rank-2 generycznie NIEZEROWY** (nematyk σ_ab istnieje — filar spinu stoi),
**rank-3 znika tożsamościowo** dla każdego sąsiedztwa i każdego C.
Powód elementarny: T_ijk jest **nieparzysty** w â, C **parzysty** ⟹ pary (â,−â) kasują się
składnik po składniku. **To nie jest artefakt siatki ani wyboru sąsiedztwa.**

> ⛔ **Trop 1 MARTWY.** ℤ₃ nie wypadnie z korelatorów substratu Isinga.
> **Jedyne, co przeżywa: korelator 4-punktowy** ⟨ŝ ŝ₊â ŝ₊b̂ ŝ₊ĉ⟩ (parzysty, 3 przesunięcia).
> Ale: nie ma go w rdzeniu jako obserwabli poziomu 0, a w MFT/Gaussie rozkłada się
> na pary (Wick) ⟹ redukowalny. **To byłby NOWY POSTULAT, nie wniosek.**

## 3. Trop 2: ℤ₃ w rdzeniu to TRZY różne grupy — BRAK POWIĄZANIA

### (a) GL(3,𝔽₂) — sprawdzone rachunkiem exact (`COLOR_z3_sources_check.py`)

```
|GL(3,F_2)| = 168        |Z(G)| = 1        |[G,G]| = 168        |G^ab| = 1
rzędy elementów: {1:1, 2:21, 3:56, 4:42, 7:48}
==> G PERFEKCYJNA;  Hom(G,U(1)) = TYLKO TRYWIALNY
```

**To jest rozstrzygające.** ℤ₃ w GL(3,𝔽₂) istnieje **tylko jako PODGRUPA** (56 elementów rzędu 3),
**nigdy jako iloraz**. A etykieta typu koloru/ładunku to **charakter**, czyli
Hom(G,U(1)) = Hom(G^ab,U(1)) — a **G^ab jest trywialna**.
⟹ **GL(3,𝔽₂) nie może nadać ŻADNEJ etykiety ℤ₃.** Strukturalnie niezdolna.

Dodatkowo `tgp_companion.tex` lin. 279–281 pokazuje, że **N=3 jest WEJŚCIEM, nie wynikiem**:
> „168 = (2N+1)·2^N·N,  **N = 3**"

To jest **faktoryzacja po podstawieniu N=3**, nie wyprowadzenie N=3. (Zgodnie z zakazem
z RAMY §1 — `why_n3 „N=3"` był już wcześniej zdebunkowany.)

### (b) Trzy ℤ₃ obok siebie

| Obiekt | Typ | Daje etykietę? | Rola |
|---|---|---|---|
| **2T^ab** (moja warstwa wstążek) | **ILORAZ** 2T/Q₈ | ✅ TAK → 3 charaktery | π₁ defektu |
| centrum SU(3) (`status_map`:425) | CENTRUM Z(SU(3)) | nie jest ilorazem π₁ | grupa cechowania |
| ℤ₃ < GL(3,𝔽₂) | **PODGRUPA** | ❌ NIE (G perfekcyjna) | symetria dyskretna V(g) |

**Trzy różne grupy, trzy różne role, na trzech różnych poziomach.**

Przeszukanie **całego** `TGP/TGP_v1/**/*.tex` na „tetraedr|trialn|2T|binarn" → **jedyne trafienia
w niepowiązanym artykule n-body**. ⟹ **W rdzeniu NIE MA śladu 2T, trialności ani abelianizacji.**

> ⛔ **Trop 2 MARTWY: BRAK POWIĄZANIA.** Most nie istnieje. Warstwa 2T jest w całości **moja**.

### (c) Ponadto: rdzeń sam nie wyprowadza koloru

`status_map`:473–475 — o łańcuchu Koide:
> `Propozycja [AN+NUM+POST(ℤ₃)]` … „Koide Q_K=3/2 z ℤ₃ (**założona, nie wyprowadzona**)"

Rdzeń **sam się przyznaje**, że ℤ₃ jest tam postulatem.

## 4. Znalezisko uboczne (istotne): rdzeń ma DWA NIEZGODNE SUBSTRATY

`core/sek09_cechowanie.tex` lin. 554–576, `ax:color-substrate` — **dosłownie**:
> „Substrat dwójkowy … **rozszerzamy do trójskładnikowego**. Każdy węzeł i ∈ V nosi
> **trójskładnikowy** stan Ξ_i = (ψ^r, ψ^g, ψ^b) ∈ **ℂ³** … Hamiltonian jest niezmienniczy
> względem globalnej **SU(3)_c**"
> H_Γ^(SU3) = Σ_i [m₀²/2·|Ξ_i|² + λ₀/4·|Ξ_i|⁴] − J_c Σ_⟨ij⟩ Re(Ξ_i†Ξ_j)

| | `dodatekB` (sektor metryczny) | `sek09` (sektor koloru) |
|---|---|---|
| stan węzła | **ŝ_i ∈ ℝ** (skalar) | **Ξ_i ∈ ℂ³** |
| symetria | **ℤ₂** (ŝ→−ŝ) | **SU(3)_c** |
| Φ | ⟨ŝ²⟩ | **niezdefiniowane** |

To **nie jest rozszerzenie — to podmiana**. ℤ₂ nie jest podgrupą SU(3)_c działającą w ten
sam sposób; ℝ nie jest podprzestrzenią ℂ³ z zachowaniem Φ=⟨ŝ²⟩. **Kolor w TGP jest
POSTULOWANY jako osobny substrat**, a nie wyprowadzony z tego, który daje metrykę i spin.

Możliwe pogodzenie (Φ = ⟨|Ξ|²⟩) **nie jest w rdzeniu wykonane** — to otwarte, sprawdzalne
zadanie, nie fakt. **Nie buduję na tym.**

## 5. Co to robi z osią §0 (spin i kolor z jednej topologii)

**Oś §0 jest teraz pod realnym naciskiem — i to jest uczciwa treść tego wyniku.**

| | Nośnik | Status |
|---|---|---|
| **Spin-½** | σ_ab = ⟨ŝŝ₊â⟩^TF → RP² (nematyk uniaxial), poziom 0 | ✅ **stoi** (rank-2 generycznie ≠ 0, potwierdzone) |
| **Kolor ℤ₃** | — | ❌ **brak nośnika w substracie Isinga** |

Z jednego substratu skalarnego z ℤ₂ **da się** dostać spin (rank 2 wystarcza),
**nie da się** dostać koloru (wymaga rank 3, który znika z dwóch niezależnych powodów).

**To nie znaczy, że oś §0 jest fałszywa** — znaczy, że **nie zamyka się na obecnym substracie**.
Rozwidlenie jest ostre i wyliczalne:

1. **Substrat pozostaje skalarny z ℤ₂** ⟹ kolor musi przyjść skądinąd (postulat, jak w sek09)
   ⟹ **oś §0 upada** (spin i kolor z różnych źródeł).
2. **Substrat ma bogatszy stan węzła** (np. Ξ ∈ ℂ³ jako **jedyny**, z Φ=⟨|Ξ|²⟩)
   ⟹ oś §0 **może** stać, ale wymaga przebudowy `dodatekB` i wykazania, że sektor metryczny
   i klasa 3D Ising przeżywają. **Duża ingerencja w rdzeń, dziś niewykonana.**
3. **Kolor nie jest własnością pojedynczego obiektu, tylko relacji** — spójne z Twoją ontologią
   (energia relacyjna, ładunek z ruchu względnego). Wtedy szukanie ℤ₃ w **jednopunktowym**
   parametrze porządku jest **szukaniem w złym miejscu**. **Nietknięte, nieprzetestowane.**

## 6. Czego NIE twierdzę

- Nie twierdzę, że TGP nie może mieć koloru. Twierdzę, że **nie wypada on z korelatorów
  substratu ŝ** i że **żaden ℤ₃ obecny w rdzeniu nie jest ℤ₃ z 2T**.
- Nie twierdzę, że korelator 4-punktowy nie zadziała — **nie testowałem go**. Wiem tylko,
  że w MFT jest redukowalny i że **nie ma go w rdzeniu**.
- Nie twierdzę, że sek09 jest błędny — jest **wewnętrznie spójny w swoim sektorze**,
  ale **nie łączy się** z substratem metrycznym.

## 7. Jednozdaniowo

> Substrat TGP to skalar Isinga z twardą ℤ₂; jego korelatory dają **rangę 2 i tylko rangę 2**
> (nematyk → RP² → spin-½ ✅), a ranga 3 potrzebna do porządku tetraedrycznego znika
> tożsamościowo z dwóch niezależnych powodów; żaden z trzech ℤ₃ obecnych w rdzeniu
> (GL(3,𝔽₂) — perfekcyjna, centrum SU(3) — postulat, Koide — jawnie „założona")
> nie ma powiązania z ℤ₃ z abelianizacji 2T; **kolor w TGP jest dziś postulowany, nie wyprowadzony.**

---

**Źródła (u źródła):** `axioms/substrat/dodatekB_substrat.tex` lin. 48–58 (H_Γ, ℤ₂), 126, 228–232;
`core/sek09_cechowanie/sek09_cechowanie.tex` lin. 554–576 (ax:color-substrate, Ξ∈ℂ³);
`core/_meta_latex/status_map.tex` lin. 422–428 (SU(3) z π₂ + centrum ℤ₃), 471–477 (POST(ℤ₃));
`tgp_companion.tex` lin. 219–222, 268–290 (GL(3,𝔽₂), 168=(2N+1)2^N·N przy N=3)
**Rachunek:** `COLOR_z3_sources_check.py` (GL(3,𝔽₂) exact nad 𝔽₂), `COLOR_rank3_nogo.py` (no-go rank-3)
