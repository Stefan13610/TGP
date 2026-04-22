# atom_from_soliton — program badawczy

**Session:** 2026-04-21
**Cel:** Pokazać że TGP ma pełny łańcuch derywacji widma atomu wodoru z samego substratu ψ_i = φ_i·e^(iθ_i), i zrewidować przedwczesny werdykt `atomic_shells_closure`.

---

## Motywacja

Użytkownik trafnie wskazał logiczną lukę w moim poprzednim rozumowaniu:

> „Jeśli mamy opis solitonów, to dlaczego jest problem z opisem atomowym?
> Jądro jako studnia potencjału zamkniętych tam mas, równowaga generowanej
> przestrzeni — nawet niestabilność jąder atomowych powinna wynikać z samych
> równań solitonowych."

Przegląd korpusu pokazał że wcześniejszy werdykt `atomic_shells_closure`
testował **zły sektor**:

- **Sektor amplitudowy** (Yukawa, real-field Φ(r)·g(r)): krótki zasięg
  ~exp(-m_sp·r). Dla atomu daje 10⁻⁶⁰ screeningu — tak, to nie daje wiązania
  atomowego.
- **Sektor fazowy** (winding, ψ = φ·e^(iθ), A_μ ∝ ∇θ): długi zasięg 1/r.
  To JEST mechanizm wiązania atomowego. **Zweryfikowany numerycznie**
  w [[research/em_from_substrate/em02_two_charge_coulomb.py]] (T1 exact,
  T5 match 0.4%).

`atomic_shells_closure` zapytało „czy A_orb z SC (empiryczny marker Fermiego)
przewiduje orbitale?". Odpowiedź: nie. **Ale to było złe pytanie.** Prawidłowe
pytanie: „czy widmo atomu emerguje z substratu ψ = φ·e^(iθ)?" — i odpowiedź
jest **tak**, przez już zamknięte łańcuchy em_from_substrate + alpha_em_rg_flow.

## Podstawy teoretyczne (z korpusu TGP)

Łańcuch który chcemy zamknąć numerycznie:

```
(1) Substrat ψ = φ·e^(iθ)                              [sek09 ax:complex-substrate]
        ↓
(2) L_phase = (J·v²·a²/2)·(∂_μθ)²                      [sek09 linia 247]
        ↓
(3) A_μ = (ℏ/e)·∂_μθ   &   1/μ₀ = 2Jv²a²e²/ℏ²          [sek09:253, eq:mu0-substrate]
        ↓
(4) Winding-1 defekt → θ(r) = -(ℏ/8πJv²a²)·(q/r)        [em02:54, poisson]
        ↓
(5) A_μ wokół defektu → Coulomb V = e²/(4πε₀R)         [em02 T1 PASS exact]
        ↓
(6) Elektron-soliton (masa m_e) w tym V → Schrödinger  [standard QM]
        ↓
(7) E_1s = -m_e·c²·α²/2 = -13.6 eV                     [Bohr]
```

Kroki (1)–(5) są **już zamknięte** w korpusie. Brakuje tylko numerycznego
złożenia (6)+(7) w jednym pliku który zademonstruje pełen łańcuch.

## Hipotezy robocze

**H1. Widmo wodoru emerguje z substratu**
Numerycznie: zbuduj pole Coulomba z winding-1 defektu na siatce 3D (jak em02),
włóż elektron (masa m_e, ładunek -e), rozwiąż równanie Schrödingera radialne,
sprawdź E_1s = -13.6 eV. **Oczekiwanie**: PASS (nic więcej nie może wyjść,
skoro Coulomb i m_e są ustalone).

**H2. Skalowanie Z² dla He²⁺, Li³⁺ itp.**
Defekt o windingu n=Z daje potencjał -Ze²/(4πε₀r). E_1s(Z) = -Z²·13.6 eV.
**Oczekiwanie**: PASS trywialnie.

**H3. Struktura powłok z fazy substratu**
Dla atomów wieloelektronowych zasada Pauli musi pochodzić z antysymetryzacji.
Skąd bierze się antysymetria w TGP? To jest **nietrywialna luka** — bo ψ_i to
pole bozonowe (kompleksowe skalarne), nie spinorowe fermionowe.
**Oczekiwanie**: FAIL w obecnym formalizmie — wymaga rozszerzenia o spinor.

**H4. Niestabilność jądrowa z wieloczastkowych równań solitonowych**
Jądro = układ N solitonów (p+n). Z V₂ + V₃ (nbody/) + ewentualnie sektor
fazowy. Czy ³H → ³He + e + ν̄ emerguje jako tunelowanie? **Oczekiwanie**:
otwarte — wymaga sektora słabego, nie tylko EM.

## Skrypty

| # | Skrypt | Cel |
|---|---|---|
| afs00 | REVISION.md | Markdown: rewizja atomic_shells verdict, logiczny łańcuch |
| afs01 | hydrogen_from_substrate.py | Pełen łańcuch (1)→(7) numerycznie dla H |
| afs02 | Z_scaling.py | He²⁺, Li²⁺, ..., Z=1...20, test Z² |
| afs03 | many_electron_gap.py | H3: pokazać że dla He, Li (neutralne) potrzeba więcej |

## Kryteria falsyfikacji

- **H1 PASS** jeśli E_1s numerycznie z substratu daje -13.6 eV z błędem < 1%
- **H2 PASS** jeśli E_1s(Z) skaluje jak Z² dla Z=1..20 z r² > 0.999
- **H3 dokument** — czysto teoretyczna luka, nie test numeryczny
- **H4** — odłożone (to jest oddzielny program `nuclear_from_soliton`)

## Scenariusze

**„Triumfalny"**: H1+H2 PASS → TGP DERYWUJE wodór + jednoelektronowe jony.
`atomic_shells_closure` werdyk zrewidowany z uczciwym komentarzem.
To jest most między atomic_shells (FAIL dla orbitali) a SC (PASS dla Tc).

**„Pokorny"**: H1 PASS ale tylko trywialnie (bo używamy α_em kalibrowanego).
To prawda — ale wciąż wartościowe bo zamyka łańcuch logiczny.

**„Otwarty"**: H3 potwierdza że fermiony/spinor to REALNA luka TGP.
Trzeba będzie osobnego programu `spinor_emergence_closure`.

---

## Pliki sesji (do utworzenia)

- `PLAN.md` — ten plik
- `REVISION.md` — rewizja atomic_shells verdict
- `afs01_hydrogen_from_substrate.py` — wodór numerycznie
- `afs02_Z_scaling.py` — skalowanie Z²
- `ATOM_FROM_SOLITON_VERDICT.md` — werdyk końcowy
