---
title: "Phase_method_decisions — decyzje metodyczne FROZEN cyklu op-collapse-matter-source (formy CYTAT, człon materii 𝒰_mat wg P1-H1, lista λ̃, ρ̂, klasyfikatory Q-H1/Q-H2, okna pomiarowe, gate'y P2)"
date: 2026-09-15
type: phase-method-decisions
tgp_owner: research/op-collapse-matter-source-2026-09-14
status: FROZEN
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_method_decisions.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_correction_note_2_energy_eval.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase_method_decisions (FROZEN przed jakimkolwiek kodem cyklu)

**Status: FROZEN 2026-09-15, ZERO obliczeń cyklu przed zapisem tego
dokumentu.** Jedyny dopuszczony uzupełniany element: §2 — forma 𝒰_mat
zostanie wpisana WYNIKIEM wyprowadzenia P1-H1 (sympy, analityka — nie
bieg numeryczny) PRZED pierwszym biegiem numerycznym; forma OCZEKIWANA
jest już zamrożona w LOCK §1. Wszystkie pozostałe decyzje zamrożone;
zmiana po pierwszym biegu produkcyjnym = forbidden move (LOCK §6).
Pozycje wykraczające poza literę LOCKa oznaczone **[INPUT-MD]**.

---

## 1. Formy sektora pola (CYTATY — dziedziczone bez modyfikacji; LOCK §1)

Źródło cytatów: `../op-r3-stationary-states-2026-09-14/Phase1_output.txt`
(formy zweryfikowane sympy, Q-D1-PASS łańcucha; konwencja |g^tt|):

- „M(psi)   = psi^6/(4-3psi)^2 = psi**6/(3*psi - 4)**2" ⟹
  **M(ψ) = ψ⁶/(4−3ψ)²**; „weryfikacja M' vs cytat
  12psi^5(2-psi)/(4-3psi)^3: simplify = 0" ⟹ **M′ = 12ψ⁵(2−ψ)/(4−3ψ)³**
- „K(psi)   = psi**4" ⟹ **𝒦(ψ) = ψ⁴**, 𝒦′ = 4ψ³
- „U(psi)   = psi**4/4 - psi**3/3" ⟹ **𝒰(ψ) = ψ⁴/4 − ψ³/3**;
  „U'  = psi**3 - psi**2  (= psi^2(psi-1): 0)"; „U'' = 3*psi**2 - 2*psi;
  U''(1) = 1"
- **K_geo = γ = c₀ = 1 [LOCK §1]**; dziedzina ψ∈(0,4/3); pas
  klasyfikacyjny ψ ≥ 4/3−1e−6 (górny) / ψ ≤ 1e−6 (dolny); ZERO
  podłóg/barier poza pasem.
- Dziedziczona tożsamość bez kancelacji (correction note 2 poprzednika,
  `Phase_correction_note_2_energy_eval.md`):
  **𝒰(ψ)−𝒰(1) = (ψ−1)²(3ψ²+2ψ+1)/12** — używana w ewaluatorze energii.

EOM z materią (LOCK §1): Mψ̈ + ½M′ψ̇² = (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′²
− 𝒰′ − ∂𝒰_mat/∂ψ.

## 2. Człon materii (z RDZENIA; forma = WYNIK P1-H1)

Zapis literalny rdzenia (sek08a, TYLKO ODCZYT): `eq:L-mat-unified`
L_mat = −(q/Φ₀)ψρ (czynnik ψ = konsekwencja elementu objętości,
decyzja post-audit A4 2026-05-01, Option 2 — preserve axiom);
ρ ≡ −T^μ_μ/c₀² ≥ 0 (L01 formal definition 2026-05-04);
√−g = c₀ψ/(4−3ψ) (`eq:vol-element-M911`, canonical M9.1'').

**Forma 𝒰_mat — wpisana z WYNIKU P1-H1 (sympy, gate 1e−12) przed
pierwszym biegiem numerycznym:**

> **[WYNIK P1-H1, 2026-09-15 — `Phase1_output.txt`, wszystkie gate'y
> PASS (simplify=0, kontrola numeryczna ≤3.55e−14):**
> **𝒰_mat(ψ,r) = λ̃·ρ̂(r)·ψ²/(4−3ψ)** (tożsame z −√−g·L_mat przy
> λ̃=(qc₀/Φ₀)ρ₀); **∂𝒰_mat/∂ψ = λ̃ρ̂·ψ(8−3ψ)/(4−3ψ)²**;
> tożsamość próżniowo odjęta (bez kancelacji, do ewaluatora energii):
> **𝒰_mat(ψ,r)−𝒰_mat(1,r) = λ̃ρ̂·(ψ−1)(ψ+4)/(4−3ψ)**.
> Wzorzec liniowy (P1-H2): δψ_lin(r) = −(5λ̃/r)∫₀^∞ r′ρ̂(r′)
> [e^{−|r−r′|}−e^{−(r+r′)}]/2 dr′; δψ_lin(0) = −5λ̃·I₀,
> I₀ = 0.776061934827. Uzupełnienie wykonane PRZED pierwszym biegiem
> numerycznym; formy od tego punktu FROZEN (forbidden move).

λ̃ ≡ (q c₀/Φ₀)·ρ₀ > 0 (bezwymiarowa siła sprzężenia; λ̃<0 poza
zakresem, bo ρ≥0). **Lista λ̃ [LOCK, FROZEN]:** Q-H1:
{0.01, 0.05, 0.2, 0.5}; Q-H2: {0.05, 0.2, 0.5}.
**ρ̂(r) = exp(−r²/18) [LOCK, FROZEN]** (gauss σ_ρ=3), statyczne —
ZERO pól dynamicznych (S05: ρ = źródło zewnętrzne).

## 3. Protokół numeryczny (dziedziczony; LOCK §3)

Silnik: kopia `../op-r3-stationary-states-2026-09-14/engine_core.py`
(z cytatem w nagłówku) — **jedyna zmiana merytoryczna: człon
−∂𝒰_mat/∂ψ w RHS (w F, punktowo: dH_i += h·r_i²·∂𝒰_mat/∂ψ|_{ψ_i,r_i})
oraz człon [𝒰_mat(ψ,r)−𝒰_mat(1,r)] w ewaluatorze energii (tożsamość
sfaktoryzowana §2, bez kancelacji — dziedziczona lekcja correction
note 2).** Dla λ̃=0 silnik jest bitowo identyczny z poprzednikiem
(dodawanie dokładnych zer).

- Siatka: r_i=(i+½)h, h=0.05 (główna) / h=0.025 (potwierdzenia),
  R=200, N=R/h; zerowy strumień w r=0 i r=R.
- Integrator: uogólniony Störmer–Verlet na (ψ,π), π=Mψ̇, dt=0.005
  (dt/2=0.0025 w potwierdzeniach — na h=0.05 **[INPUT-MD,
  dziedziczone]**); iteracje punktu stałego ze stopem stagnacji
  maszynowej (correction note 1(b) poprzednika); NonConvergence ⟹
  BREAKDOWN.
- Sponge: smootherstep γ₀=1.0 na [160,200], działa na π; ON we
  wszystkich biegach produkcyjnych Phase 3 i w P2a/P2b; OFF w P2c
  (pomiar zachowania energii = układ zamknięty) **[INPUT-MD,
  dziedziczone]**.
- Energia: próżniowo odjęta (E[ψ≡1,π≡0]=0 również przy λ̃>0 dzięki
  odjęciu 𝒰_mat(1,r)), czynnik 4π; **E_core** = suma członów dla
  r_i≤80 (gradientowe: r_{i+½}≤80) — **gęstość PEŁNA z członem
  materii** (to ona jest zachowana) **[INPUT-MD — operacjonalizacja]**;
  deskryptywnie raportowane też E_core^field (bez członu materii).
  E_ref ≔ E_core(50).
- Zapis: ψ(0,t)≔ψ(h/2,t), E_core(t), max_{r≤80}|ψ̇| co dt_out=0.1;
  pełny profil co 5 j.cz. (Q-H1) / co 50 j.cz. (Q-H2); checkpoint
  npz co 100 j.cz. **[INPUT-MD]**.
- Starty Q-H2 (LOCK, FROZEN — IDENTYCZNE definicje jak u poprzednika,
  MD §6 tamtego cyklu): gauss ψ=1+a·exp(−r²/2σ²); quasi-R3
  ψ=1+a·sinc(r/π)·exp(−r²/(2·15²)) (sin(r)/r przez np.sinc);
  π₀≡0; źródło włączone od t=0. Reprezentanci: qR3 a=−0.20,
  qR3 a=+0.20, gauss a=−0.30 σ=3, gauss a=+0.15 σ=6.
- Baseline λ̃=0 (poprzednik, `Phase3_output.txt` + `Phase3_results/`,
  wszystkie COLLAPSE=BREAKDOWN-BOUNDARY górny, czasy t_end h05/h025):
  qR3 a=−0.20: 4.74/4.74; qR3 a=+0.20: 2.275/2.27;
  g a=−0.30 σ3: 6.64/6.635; g a=+0.15 σ6: 18.185/10.34.

## 4. Phase 1 — analityka pre-rejestrowana (sympy; zero ewolucji)

- **P1-H1:** 𝒰_mat z literalnego −√−g·L_mat = √−g·(q/Φ₀)ψρ,
  √−g=c₀ψ/(4−3ψ), ρ=ρ₀ρ̂(r), λ̃=(qc₀/Φ₀)ρ₀; gate: simplify(𝒰_mat −
  λ̃ρ̂ψ²/(4−3ψ)) = 0 oraz simplify(∂𝒰_mat/∂ψ − λ̃ρ̂ψ(8−3ψ)/(4−3ψ)²)=0;
  kontrola numeryczna |Δ|≤1e−12·max(1,|v|) w ψ∈{0.5,1,1.2}
  **[INPUT-MD — punkty kontrolne]**; dodatkowo gate tożsamości
  próżniowo odjętej (ψ−1)(ψ+4)/(4−3ψ) (§2).
- **P1-H2:** linearyzacja wokół ψ=1 ze źródłem: (−∇²+1)δψ = −5λ̃ρ̂
  (wyprowadzenie sympy z EOM; 5 = ∂²𝒰_mat/∂ψ∂λ̃|₁·ρ̂⁻¹); rozwiązanie
  δψ = −5λ̃(G_Yuk∗ρ̂), G_Yuk=e^{−r}/(4πr); redukcja sferyczna:
  **δψ(r) = −(5λ̃/r)∫₀^∞ r′ρ̂(r′)·[e^{−|r−r′|}−e^{−(r+r′)}]/2 dr′**
  (wyprowadzenie kątowe w skrypcie); ewaluacja kwadraturą
  (scipy.integrate.quad, granica górna 60 **[INPUT-MD]**) = wzorzec
  gate'u P3-H1a; kontrola rezydualna: dyskretny operator (−∇²+1) na
  δψ_lin vs −5λ̃ρ̂ (deskryptywnie). PREDYKCJA pre-rejestrowana
  [LOCK §2]: **δψ<0** wszędzie, ogon e^{−r}/r.
- **P1-H3:** granice sympy: lim 𝒰_mat przy ψ→4/3⁻ = +∞ (odpychanie
  od sufitu przy ρ̂>0); lim przy ψ→0⁺ = 0 (brak bariery od materii
  na podłodze). Pre-rejestrowana ASYMETRIA: stabilizacja łatwiejsza
  dla kolapsów górnych niż dolnych — konfrontacja w Q-H2.

## 5. Phase 2 — gate'y (FROZEN; FAIL ⟹ STOP)

- **P2a (próżnia):** λ̃=0, konfiguracja produkcyjna (sponge ON), ψ≡1,
  π≡0, 100·T₀ (T₀=2π), obie siatki; gate ‖ψ−1‖∞ ≤ 1e−10 cały bieg.
- **P2b (regresja):** λ̃=0, start qR3 a=−0.20, h=0.05, dt=0.005,
  sponge ON (konfiguracja produkcyjna poprzednika); gate: COLLAPSE
  (BREAKDOWN-BOUNDARY) z czasem zdarzenia 4.74 ±2%.
- **P2c (energia ze źródłem):** λ̃=0.05, start gauss a=+0.05 σ=3,
  sponge OFF, t=700, obie siatki **[INPUT-MD — obie siatki]**;
  dryf ≔ |⟨E⟩_{t∈[90T₀,100T₀]} − ⟨E⟩_{t∈[0,10T₀]}| / |⟨E⟩_{[0,10T₀]}|
  ≤ 1e−6 (okna dziedziczone z MD §7 poprzednika **[INPUT-MD]**);
  deskryptywnie max|E−E(0)|.

## 6. Phase 3 — Q-H1: klasyfikatory i okna (FROZEN)

Biegi: start ψ≡1, π₀=0, źródło od t=0; λ̃∈{0.01,0.05,0.2,0.5},
t_max=300, h=0.05; λ̃∈{0.01,0.5} dodatkowo h=0.025 [LOCK].

Definicje pomiarowe **[INPUT-MD — operacjonalizacje]**:
- **ψ̄(r)** ≔ średnia próbek profilowych z okna t∈[250,300] (profil
  co 5 j.cz. ⟹ 11 próbek); **δψ(0)** ≔ ψ̄(r₀)−1, r₀=h/2.
- **Osiadłość („‖ψ̇‖∞ → poziom szumu")**: V ≔ max_{t∈[250,300]}
  max_{r≤80}|ψ̇(r,t)| (ψ̇=π/M, próbki co dt_out=0.1);
  D ≔ max_{r≤80}|ψ̄(r)−1|; **osiadły ⟺ V ≤ 0.01·max(D, 1e−12)**.
- **COLLAPSE** (nadkategoria; priorytet 1): pas graniczny
  (ψ≥4/3−1e−6 / ψ≤1e−6) lub niefinityczność/NonConvergence
  w dowolnym t; podtyp raportowany (BOUNDARY górny/dolny/BREAKDOWN);
  czas zdarzenia t_end.
- **THRESHOLD-PULL** (priorytet 2): osiadły ORAZ trwale za progiem
  M911 w oknie [250,300]: min_r ψ(r,t) < 5/6 dla WSZYSTKICH próbek
  profilowych okna (lub analogicznie max_r ψ(r,t) > 7/6).
- **DEFORMATION** (priorytet 3): osiadły ORAZ wszystkie próbki
  profilowe okna w [5/6, 7/6] wszędzie.
- **INCONCLUSIVE-RUN:** reszta (nieosiadły do t_max; stany mieszane).
- **Zbieżność klasy:** dla λ̃∈{0.01,0.5} klasa h05 = klasa h025
  (dla COLLAPSE dodatkowo t_end ±10%); dla λ̃∈{0.05,0.2} klasa z h05
  (LOCK nie przewiduje tam drugiej siatki) **[INPUT-MD]**.
- **Gate liniowy (λ̃=0.01, h=0.05):** err ≔ max_{r≤40}|(ψ̄−1) −
  δψ_lin| / max_{r≤40}|δψ_lin| ≤ 0.05; δψ_lin z kwadratury P1-H2
  na węzłach siatki.
- Deskryptywnie obowiązkowo: δψ(0) zmierzone vs −5λ̃(G_Yuk∗ρ̂)(0)
  dla wszystkich λ̃ (odchylenie nieliniowe = wynik).

## 7. Phase 3 — Q-H2: klasyfikatory (FROZEN)

Biegi: 4 starty §3 × λ̃∈{0.05,0.2,0.5}, t_max=1000, h=0.05, dt=0.005.

- **COLLAPSE** (nadkategoria, podtyp raportowany): jak §6.
- **RADIATED:** bez zdarzenia brzegowego ORAZ E_core(t_max) <
  0.05·E_ref (gęstość pełna §3; wymaga E_ref>0, patrz niżej).
- **STABILIZED:** przeżywa t_max bez zdarzenia brzegowego i bez
  zaniku: E_core(t_max) ≥ 0.05·E_ref; podtypy deskryptywne:
  osiadły statycznie (kryterium osiadłości §6 na oknie [950,1000])
  vs oscylujący.
- **INCONCLUSIVE-RUN:** reszta; w tym przypadek E_ref ≤ 0 bez
  zdarzenia brzegowego (progi energetyczne tracą sens)
  **[INPUT-MD — zabezpieczenie]**.
- **τ (do zgodności potwierdzeń):** t_hold = największe t:
  E_core(t′) ≥ 0.5·E_ref ∀ próbek t′∈[50,t]; τ = t_hold−50;
  t_hold=t_max ⟹ τ cenzurowane („≥950") **[INPUT-MD, dziedziczone]**.
- **Potwierdzenia [LOCK, FROZEN]:** każdy bieg z kategorią ≠ COLLAPSE
  (baseline λ̃=0) → DWA biegi potwierdzające: h=0.025 (dt=0.005)
  ORAZ dt/2 (dt=0.0025, h=0.05); zgodność: ta sama kategoria; dla
  COLLAPSE t_end ±10%; dla STABILIZED/RADIATED τ ±10% (cenzurowane:
  zgodność ⟺ druga ≥855; obie cenzurowane ⟹ zgodne **[INPUT-MD,
  dziedziczone]**). Kontrola negatywu obowiązkowa: jeden bieg
  COLLAPSE — para (qR3 a=+0.20, najniższe λ̃ dające COLLAPSE dla
  tego startu) → h=0.025 (zgodność kategorii + t_end ±10%).

## 8. Werdykty — litera LOCK §4 (bez zmian, przywołanie)

Q-H1-DEFORMATION / Q-H1-PULL / Q-H1-INCONCLUSIVE;
Q-H2-PASS (≥1 para STABILIZED potwierdzona h/2 i dt/2, baseline
COLLAPSE) / Q-H2-FAIL (wszystkie 12 par COLLAPSE zbieżnie) /
Q-H2-INCONCLUSIVE. INCONCLUSIVE ≠ pozytyw.

## 9. Rejestr WEJŚĆ (flagowane)

[LOCK]: formy M,𝒦,𝒰 (cytaty §1); 𝒰_mat z eq:L-mat-unified (forma =
wynik P1-H1); λ̃: {0.01,0.05,0.2,0.5} (Q-H1) / {0.05,0.2,0.5} (Q-H2);
ρ̂=exp(−r²/18); starty-reprezentanci (4); h∈{0.05,0.025}; R=200;
dt=0.005; sponge [160,200] smootherstep; t_max=300 (Q-H1) / 1000
(Q-H2); progi 5/6, 7/6; pas 4/3−1e−6 / 1e−6; E_core r≤80;
E_ref=E_core(50); progi 0.05·E_ref i 0.5·E_ref; gate liniowy ≤5%
na r≤40 (λ̃=0.01); P2a 1e−10; P2b 4.74±2%; P2c ≤1e−6/100T₀;
potwierdzenia h/2+dt/2, ±10%; kontrola negatywu qR3+0.20.
[INPUT-MD]: γ₀=1.0 (dziedziczone); sponge OFF w P2c; P2c obie
siatki; okna dryfu [0,10T₀]/[90T₀,100T₀]; dt_out=0.1; profil co 5
(Q-H1) / 50 (Q-H2) j.cz.; checkpoint co 100 j.cz.; ψ(0,t)≔ψ(h/2,t);
gęstość E_core pełna (z 𝒰_mat, próżniowo odjęta) + E_core^field
deskryptywnie; osiadłość V≤0.01·max(D,1e−12) i okna [250,300] /
[950,1000]; trwałość progu = wszystkie próbki okna; priorytety
klasyfikacji; zbieżność Q-H1 tylko na siatkach przewidzianych
LOCKiem; kwadratura quad do 60; punkty kontrolne P1 {0.5,1,1.2};
τ/t_hold i reguła cenzurowania (dziedziczone); dt/2 na h=0.05;
E_ref≤0 ⟹ INCONCLUSIVE-RUN; tolerancja punktu stałego = stagnacja
maszynowa (dziedziczona z correction note 1(b) poprzednika).

**FROZEN. Zmiany po starcie obliczeń = forbidden move (poza
uzupełnieniem §2 wynikiem P1-H1 przed pierwszym biegiem numerycznym).**
