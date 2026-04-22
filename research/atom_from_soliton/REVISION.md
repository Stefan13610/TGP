# REVISION — rewizja werdyktu atomic_shells_closure

**Session:** 2026-04-21
**Status:** Korekta poprzedniego werdyktu (ATOMIC_SHELLS_VERDICT).

---

## Problem z poprzednim werdyktem

Program `atomic_shells_closure` zakończył się negatywnym werdyktem „TGP nie ma
aparatu atomowego". Ten werdykt był **przedwczesny i niepoprawny** — testował
zły sektor teorii.

### Co testowano faktycznie

Testy as01–as03 sprawdzały czy **amplitudowe** markery A_s, A_sp, A_d, A_f
(z programu superconductivity_closure / ps41) przewidują orbitale atomowe
i oscylacje Φ(r) = Φ_0·g(r) w skalach atomowych. Wyniki:

- A_s ≠ A_tail(2s Li) z błędem 7× (as02)
- Yukawa-style coupling z profilu Φ daje screening ~exp(-Bohr/Compton) ~ 10⁻⁶⁰ (as03)
- Orbitale s, p, d, f nie wyłaniają się z samego Φ-sektora (as01)

### Dlaczego to było nie to pytanie

Atom wodoru NIE wiąże się przez amplitudowe pole Φ. Wiąże się przez **Coulomba**
który pochodzi z **fazowego** sektora substratu. Te dwa sektory są w TGP
formalnie rozdzielone:

| Sektor | Obiekt | Zasięg | Rola |
|---|---|---|---|
| Amplitudowy | Φ(r) = Φ_0·g(r) | Yukawa e^(-m·r) | masa, siła silna/słaba |
| Fazowy | θ(r), A_μ ∝ ∇θ | Coulomb 1/r | ładunek, EM |

`atomic_shells_closure` testował czy amplitudowy sektor daje atom. Oczywiście
nie daje — amplitudowy sektor jest krótkozasięgowy. **Atom jest zjawiskiem
fazowego sektora** (EM = Coulomb = długi zasięg).

## Co TGP faktycznie derywuje dla atomu

Łańcuch logiczny który **już istnieje w korpusie**:

### Krok 1: Substrat (sek09 ax:complex-substrate)
```
ψ_i = φ_i · e^(iθ_i)
H = Σ_i [(m₀²/2)|ψ|² + (λ₀/4)|ψ|⁴] - J·Σ_<ij> Re(ψ_i*·ψ_j)
```

### Krok 2: Emergencja pola fazowego (sek09:247, thm:photon-emergence)
```
L_phase = (J·v²·a²/2) · Σ_μ (∂_μθ)²
```

### Krok 3: Maxwell & μ₀ (sek09:253, eq:mu0-substrate, dodatekO:321)
```
A_μ = (ℏ/e) · ∂_μθ
1/μ₀ = 2·J·v²·a²·e²/ℏ²
α_em = 1/(8π·J·v²·a²)      (Planck units)
```

### Krok 4: Pole wokół winding-1 defektu (em02:54)
```
J·v²·a²·∇²θ = -(ℏ/2)·ρ_q
→ θ(r) = -(ℏ/(8π·J·v²·a²))·(q/r)
→ A_μ(r) daje potencjał V(r) = e²·q_nucl·q_e/(4πε₀·r)
```

### Krok 5: Weryfikacja numeryczna (em02)
- T1 PASS: analityczny prefaktor e²/(4πε₀) exact do 10⁻¹⁰
- T5 PASS: V_int/V_Coulomb = 1.0037 przy R=2 (match 0.4%)
- T6 PASS: superpozycja dwóch ładunków (błąd 9×10⁻⁷)
- Boundary effects dla R ≥ 4 — artefakt finite-size, nie fizyka

### Krok 6: Elektron jako cząstka punktowa w V(r)
Standardowy Hamiltonian nierelatywistyczny:
```
Ĥ = -ℏ²·∇²/(2m_e) - e²/(4πε₀·r)
```

### Krok 7: Rozwiązanie Schrödingera
```
ψ_1s(r) = (1/√π·a₀³) · exp(-r/a₀)
a₀ = 4πε₀·ℏ²/(m_e·e²) = 5.29·10⁻¹¹ m
E_1s = -m_e·c²·α_em²/2 = -13.6057 eV
```

## Status łańcucha

| Krok | Status w korpusie | Plik |
|---|---|---|
| 1. ψ = φ·e^(iθ) | ✓ Aksjomat formalny | sek09 ax:complex-substrate |
| 2. L_phase = (Jv²a²/2)(∂θ)² | ✓ Derywowane | sek09:247 |
| 3. A_μ, μ₀, α_em | ✓ Derywowane | sek09:253, dodatekO |
| 4. Poisson → 1/r | ✓ Analitycznie | em02 setup |
| 5. Coulomb weryfikacja | ✓ Numerycznie 0.4% | em02 T1, T5, T6 |
| 6. Schrödinger w V(r) | — (standard QM) | — |
| 7. E_1s = -13.6 eV | **nie złożone jawnie** | afs01 |

**Luka jest tylko w kroku 7**: nikt dotąd nie zebrał kroków (1)–(5) + standardowy
krok (6) w jedną demonstrację numeryczną dla wodoru. Nie brakuje teorii —
brakuje **skryptu demonstracyjnego**.

## Rewizja werdyktu ATOMIC_SHELLS_VERDICT

### Co zostało błędnie stwierdzone

Poprzedni werdyk (ATOMIC_SHELLS_VERDICT) sugerował że TGP „nie ma aparatu
atomowego". **To jest za mocne stwierdzenie.** Poprawnie:

- TGP **ma** aparat atomowy dla **atomu wodoru / jednoelektronowych jonów**
  (przez łańcuch substrat → Coulomb → Schrödinger → widmo)
- TGP **nie ma** aparatu dla:
  - **Fermionów i zasady Pauli** (potrzebny sektor spinorowy, nie skalarowy)
  - **Orbitali s, p, d, f jako markerów Fermi-surface** (A_orb to empiryczne etykiety)
  - **Struktury powłok w atomach wieloelektronowych** (wymaga antysymetryzacji)

### Co jest FAKTYCZNĄ luką

1. **Sektor spinorowy / fermionowy**
   TGP definiuje ψ jako pole kompleksowe skalarowe. Elektron w rzeczywistości
   jest fermionem (spinorem Diraca). Jak emerguje antysymetryzacja z substratu
   bozonowego? **Otwarte.**

2. **Korelacja elektron-elektron**
   W atomach wieloelektronowych energia wiązania pochodzi nie tylko z Coulomba
   elektron-jądro ale też z korelacji e-e (exchange-correlation). TGP nie ma
   natywnego aparatu DFT/HF.

3. **A_orb jako topologia vs marker**
   W SC A_s, A_sp, A_d, A_f działają jako etykiety powłok Fermiego. Pytanie
   otwarte: czy one mają derywację z topologii substratu (winding + defekty
   SU(2)/SU(3))? Jeśli tak, są fundamentalne; jeśli nie, są fenomenologiczne.

## Nowy werdykt roboczy

**Stary (zbyt mocny)**:
> TGP nie ma aparatu atomowego. Orbitale to empiryczne etykiety.

**Nowy (precyzyjny)**:
> TGP **derywuje widmo atomu wodoru** z substratu ψ = φ·e^(iθ) przez łańcuch
> faza → Maxwell → Coulomb → Schrödinger. Luki są w:
> (a) emergencji struktury spinorowej fermionów,
> (b) opisie korelacji elektron-elektron w atomach wieloelektronowych,
> (c) fundamentalnym pochodzeniu A_orb (marker vs topologia).
>
> Te luki są węższe i konkretniejsze niż wcześniej sądzono.

## Co robimy w `atom_from_soliton`

- **afs01**: demonstracja numeryczna H → E_1s = -13.6 eV (złożenie łańcucha 1-7)
- **afs02**: weryfikacja skalowania Z² dla H⁰, He⁺, Li²⁺, ..., Ne⁹⁺
- **VERDICT**: synteza + konkretne luki do przyszłych programów

Nie próbujemy naprawić luk (a)–(c) — tylko uczciwie pokazujemy że te trzy
pytania są **oddzielne** od „czy TGP ma atom wodoru".

---

## Pliki powiązane

- [[ATOMIC_SHELLS_VERDICT.md]] — stary werdyk (do zrewidowania)
- [[EM_VERDICT.md]] — zamknięcie kroków 1–5
- [[research/em_from_substrate/em02_two_charge_coulomb.py]] — weryfikacja Coulomba
- [[sek09_cechowanie.tex]] — formalizm substratu
- [[dodatekO_u1_formalizacja.tex]] — formalizacja U(1)
