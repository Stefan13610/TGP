---
title: "L07 — Warunek zerowej sumy: aksjomat vs derywacja"
date: 2026-05-06
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L07_zero_sum_axiom
tags:
  - audit
  - ontology
  - dark-energy
  - zero-sum
  - axiom
  - L07
  - EXT-2
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-zero-sum-derivation", "op-nonlocal-foundations"]
  depends_on: []
  impacts: ["sek05_ciemna_energia (Λ_eff > 0)", "ax:zero w sek01_ontologia"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# L07 — Warunek zerowej sumy: aksjomat vs derywacja

## Klasa: LUKA ONTOLOGICZNA (zewnętrzna recenzja EXT-2) • Priorytet: **P2**

## Diagnoza (z EXT-2)

`sek05_ciemna_energia.tex` lin. 240–293 (`prop:Lambda-positive`)
opiera **całość mechanizmu Λ_eff** na warunku zerowej sumy:

```
∫_Σ φ √h d³x = 0     (eq:zero-sum-phi, lin. 250)
```

zwany **"zasadą zerowej sumy"** (`ax:zero` w sek01).

Bez tego warunku ⟨φ²_min⟩ niekoniecznie > 0, a Λ_eff = (8πG_0/c_0⁴)
⟨U(φ_min)⟩ niekoniecznie > 0. **Cała sek05 (ciemna energia z fluktuacji)
zawisa na tym jednym wymaganiu.**

Warunek jest **globalny** — więz na całej hypersurface przestrzennej Σ.
Nie wynika z lokalnej dynamiki (Φ-EOM), nie wynika z H_Γ jawnie.

## Pliki dotknięte

- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]
  lin. 146–161 (subsection "Związek z zasadą zerowej sumy"),
  240–293 (prop:Lambda-positive)
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] aksjomat ax:zero
  (linia do zlokalizowania w pełnym pliku 1435 lin.)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] lin. 736+
  (H_Γ definition v2 GL-bond)
- `axioms/substrat/` (substrat directory — referencje do roli Z₂)

## Problemy jakie rodzi

1. **Lokalność vs globalność.** Globalny więz bez Lagrange'a multiplier
   z dynamiki to **klasyczna czerwona flaga w fizyce teoretycznej**.
   Albo jest konsekwencją gauge/redundancji (jak constraint Hamiltoniana
   w GR z dyfeomorfizmu), albo wprowadza nielokalność w postaci
   "system jako całość zna globalny rozkład materii".

2. **Operacyjna nieobserwowalność.** Warunek dotyczy całej hypersurface
   Σ jednocześnie. Lokalny eksperymentator nie może go zweryfikować ani
   złamać. To czyni go **nie-falsyfikowalnym** — niezgodne z duchem
   TGP (`TGP_FOUNDATIONS.md` § 5.3: "GR jako numeryczny analog, nie
   izomorfizm; otwiera drzwi falsyfikacji").

3. **Mechaniczna ad-hoc-owość.** Wprowadzenie ax:zero jako *aksjomat*
   przyznaje, że jest on potrzebny, ale nie wynika z głębszej
   struktury. TGP ma 3 fundamenty:
   - jedno pole Φ z Z₂ (§1 FOUNDATIONS)
   - hamiltonian H_Γ z GL-bond
   - **warunek zerowej sumy ∫ φ √h = 0**

   Trzeci wygląda na *montowanie* zamiast *wyprowadzenia*.

4. **Kompatybilność z lokalnymi rozwiązaniami.** Jeśli rozważam
   izolowaną cząstkę w TGP (radialny kink Φ(r)), to lokalnie φ > 0
   (rośnie wokół źródła). ∫ φ √h > 0 dla każdej skończonej kuli. Więz
   ∫_Σ φ = 0 wymaga, by **gdzieś indziej** φ był ujemny — ale
   `sek02_pole` lin. 116–117 explicit zakazuje φ < 0 ("Φ(r) > 0 wszędzie,
   przestrzeń jest generowana"). **Sprzeczność?** Albo "deficyt φ < 0"
   oznacza coś innego niż dosłowna ujemność — ale wtedy potrzebne jest
   doprecyzowanie operacyjne.

5. **Cosmological constant problem zwija się tylko częściowo.**
   Λ_eff ~ γ/12 ~ Φ_0 H_0²/(12 c_0²). Wymaga γ w skali H_0² — ale
   *dlaczego* γ jest w skali H_0², a nie M_Pl²? Standardowy CC problem
   (10¹²⁰ rozbieżność QFT vacuum vs obs) został przesunięty z "dlaczego
   Λ jest takie małe" na "dlaczego γ jest takie małe". **Bez
   wyprowadzenia γ z UV, problem wraca tylnymi drzwiami.**

## Potencjalne ścieżki domknięcia

### Ścieżka A — derywacja zerowej sumy z Z₂-symetrii substratu

**Hipoteza:** ∫_Σ φ √h = 0 może być konsekwencją tożsamości typu
Bianchi dla pola φ na grafie Γ. Konkretnie: jeśli H_Γ = −J Σ (φ_i φ_j)²
ma chiralną Z₂ (`sek01_ontologia.tex:83`), to operator parzystości
P: φ_i ↔ −φ_i pozostawia H_Γ invariantnym. Suma ∫ φ √h nad
hypersurface zachowującą P powinna być zerem z konstrukcji.

**Cykl proponowany:** `op-zero-sum-derivation/` z explicit derywacją
∫ φ = 0 jako tożsamość Z₂ (nie aksjomat).

### Ścieżka B — Lagrange'a multiplier w akcji zunifikowanej

Dodać do S_TGP człon λ ∫ φ √h. Wariacja po λ daje warunek; po φ daje
EOM z λ jako efektywnym Λ. To by *zlikwidowało* ax:zero jako aksjomat,
redukując go do dynamicznego więzu z explicit mnożnikiem. Cena: λ
jest dodatkowym parametrem, którego skala wymaga wyjaśnienia.

### Ścieżka C — Wprowadzenie φ_eff = φ − ⟨φ⟩_Σ

Reinterpretacja: "zerowa suma" jest *definicją* φ_eff jako odchylenia
od średniej. Wówczas ∫ φ_eff = 0 jest tożsamością. Λ_eff przestaje
wynikać z warunku, a zaczyna wynikać z ⟨φ⟩_Σ — co jest mierzalne
(tło kosmologiczne). Wymaga przekształcenia całej sek05.

### Ścieżka D — domknięcie przez nielokalność

Przyznanie, że TGP jest **nielokalna na skali kosmologicznej** (nie
sprzeczne z relativity jeśli nielokalność jest *spacelike* na
hypersurface), z explicit dyskusją horyzontów i kompatybilności
z FRW. Cykl `op-nonlocal-foundations/`.

## Rekomendowany priorytet

**P2 — wysoki.** Bez ścieżki A lub B, sek05 (główna predykcja
kosmologiczna TGP) opiera się na aksjomacie, którego status nie jest
ani derived, ani falsyfikowalny.

## Powiązanie z istniejącym audytem

**Nowa klasa L07** (proponowana w EXT-2). Brak odpowiednika w obecnym
S/L/D/M. Skala konsekwencji: mechanizm ciemnej energii w TGP.

## Cross-references

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2 — recenzja źródłowa
- [[../README.md]] — indeks audytu
- [[../PRIORITY_MATRIX.md]] — do update z L07 P2
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]
  prop:Lambda-positive
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] ax:zero
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] H_Γ definition
- [[../../TGP_FOUNDATIONS.md]] § 5.3 falsyfikacja
