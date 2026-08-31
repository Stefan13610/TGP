# Phase FINAL — zamknięcie: op-substrate-fluctuation-channel-2026-08-23

**Status cyklu: CLOSED-EXECUTED (2026-08-23, jedna sesja: LOCK → Phase 1–3).**
Kryteria: [[Phase0_balance.md]] (LOCK + Amendment A1 zamknięte przed kodem).

---

## Werdykty

### QF — kanał fluktuacyjny: **PASS** (QF-1 ∧ QF-2 ∧ QF-3 ∧ QF-4)

Na poziomie Gaussowskim substratu (poziom 0) istnieją DOKŁADNIE trzy kanały
oddziaływania defekt–defekt i mają rozdzielne sygnatury (Phase 1 sympy 9/9;
Phase 2 siatka L=64/128 FFT):

| Kanał | Forma (exact) | Znak | Zasięg | Krytyczność |
|---|---|---|---|---|
| M-S źródłowy | −q₁q₂·G_m(d) | ~sgn(q₁q₂) — **ładunkowy** | μ(m) | ~−1/d |
| M-D klasyczny | ≈ −(v₁v₂/G₀²)·G_m(d) | ~sgn(v₁v₂) — **ładunkowy** | μ(m) | ~−1/d |
| **M-D fluktuacyjny** | **F_fl = ½ln(1−(G(d)/G₀)²)** | **< 0 ZAWSZE (uniwersalne przyciąganie, v-niezależne)** | **2μ(m)** | **~−1/d²** |

Numerycznie: κ_S/μ i κ_D/2μ zgodne do 0.8–4.0% (R²≥0.9998); krytyczne
p = −1.077/−1.029 (L=64/128, dryf 0.048), slope F_fl = −2.15/−2.06;
tabela znaków 6/6 zgodna (kanały klasyczne flipują, fluktuacyjny nie);
kontrole negatywne przeszły (zanik m=2: ratio 8.6e−8; obrazy periodyczne
6.1e−4; pozytywność G_m).

**Incydent QF-4c (udokumentowany, wzorzec A3):** pierwszy przebieg dał
FAIL-as-implemented — G_fft(d=22; m=2, L=64) = −2.33e−18, tj. wartość POD
progiem zaokrągleń FFT (eps·G(0)=2.38e−17). Diagnoza Phase 2b: niezależna
suma spektralna mpmath dps=40 daje G(22) = +5.397e−20 > 0 (sanity metody:
zgodność z FFT dla d=4 do 3.9e−14). Błąd implementacji TESTU (porównanie
szumu maszynowego do zera) ZNALEZIONY ⟹ korekta dopuszczalna; kryterium
LOCK (G_m>0) nietknięte. Pełny ślad: Phase2_output.txt (FAIL pierwotny,
niezamazany) + Phase2b_output.txt (rozstrzygnięcie).

**Interpretacja (inwentarz, zero claimów obserwacyjnych):** kanał
fluktuacyjny jest **jedynym kanałem poziomu 0 o uniwersalnym znaku
przyciągania** — sygnatura, której brakowało kanałom przetestowanym
w programie „most do grawitacji" (amplitudowy: masywny; Goldstone:
ładunkowy). Zastrzeżenie jawne (LOCK §6): dla defektów PUNKTOWYCH na
krytyczności daje potencjał ~−1/d², NIE Newtonowskie −1/d — skalowanie dla
obiektów rozciągłych i N-ciałowość = NEEDS N1/N2. Kontekst skali z korpusu:
m_eff² = γ(1+T_Γ) (cor:entropy-potential), γ ≈ 12Λ_eff/Φ₀ ⟹ ξ = 1/√γ
rzędu horyzontu ⟹ wewnątrz horyzontu substrat jest w reżimie efektywnie
krytycznym (potęgowym), spójnie z warunkiem continuum
prop:continuum-conditions — połączenie tych faktów z inwentarzem kanałów
było celem cyklu (NEEDS N3).

### QB-1 — znak wkładu wiązania: **ΔC_bond = +8zJ·s_b⁶ ≥ 0 ZAWSZE** (sympy exact)

Wiązanie gradientowe v2 STABILIZUJE punkt jednorodny dla każdej gęstości
kąpieli. **Znak tachionowy NIE MOŻE emergować z wiązania gradientowego
substratu na poziomie MFT** — wynik negatywny zgłoszony wprost; zawęża
pochodzenie znaku tachionowego W (pcha na poziom 1 — testowany przez
op-bath-two-sectors Q2 — lub poza-MFT, NEEDS N4).

### QB-2 — próg rozrzedzenia: **ISTNIEJE (spinodala substratu)**

C(σ)/J = r + 3uσ + 48σ³ (z=6, σ=Φ_b): jedyny korzeń (dC/dσ = 3u+144σ² > 0).
Przy punkcie WF (r*=−2.251, u*=3.917, eq:B-WF): **Φ_c/Φ_vac = 0.298** —
kąpiel o gęstości próżniowej stabilna (C=13.61>0), kąpiel rozrzedzona
poniżej ~30% gęstości próżniowej ⟹ krzywizna tachionowa (źródło: część
ON-SITE, nie wiązanie; bond zawęża obszar niestabilny: σ_c=0.171 < 0.192
bez bondu). Skan r∈[−3,−1]×u∈[2,6] (41×41): Φ_c/Φ_vac ∈ [0.197, 0.331]
wszędzie w (0,1). Phase 3: 9/9 PASS (w tym kontrola negatywna flip-J;
incydent QB-2a pierwszego przebiegu = literówka stałej oczekiwanej
w teście, 288→144, udokumentowana w skrypcie; pochodna i wnioski bez zmian).

### QB-3 — odczyt (deskryptywnie, zero claimów poziomu 1)

Wzór {pełna kąpiel → stabilnie; rozrzedzenie → tachionowo} jest
strukturalnie równoległy do hipotezy dwóch sektorów N2, ale z **odwrotną
rolą gęstości**: na poziomie 0 destabilizuje NISKA gęstość tła (spinodala),
podczas gdy hipoteza N2 wiąże tachion z OBECNOŚCIĄ kąpieli. To napięcie
jest informacją, nie sprzecznością (inne maszynerie: statystyka on-site
vs dynamika solitonowa ODE); rozstrzygnięcie należy do poziomu 1
(op-bath-two-sectors), nie do tego cyklu.

## Anti-Lakatos

✓ LOCK + Amendment A1 zamknięte PRZED kodem; kryteria/okna nietknięte po
starcie obliczeń. ✓ Dwa incydenty testowe (QF-4c precyzja, QB-2a literówka)
zdiagnozowane z błędem implementacji ZNALEZIONYM przed werdyktem; pierwotne
outputy FAIL zachowane niezamazane. ✓ Wynik negatywny (QB-1: tachion nie
z wiązania) zgłoszony wprost. ✓ Kontrole negatywne wykonane (T5, QF-4a/b/c,
QB-1d). ✓ Każdy skrypt z testami zdolnymi dać FAIL (2 faktyczne FAIL-e po
drodze to potwierdzają). ✓ Zakaz dublowania op-bath-two-sectors dotrzymany
(poziom 0 vs 1). ✓ Rdzeń .tex NIETKNIĘTY (NEEDS user-gated). ✓ NIE commitowano.

## Deliverables

Phase0_balance.md (LOCK+A1) · Phase1_analytic.py/.txt (9/9) ·
Phase2_lattice.py/.txt (8/9 + werdykt maszynowy) · Phase2b_positivity_check.py/.txt
(rozstrzygnięcie QF-4c) · Phase3_bath_curvature.py/.txt (9/9) · NEEDS.md · README.
