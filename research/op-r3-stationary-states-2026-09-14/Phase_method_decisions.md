---
title: "Phase_method_decisions — decyzje metodyczne FROZEN cyklu op-r3-stationary-states (integrator 2. rzędu na (ψ,π), sponge, detektor oscylonu, definicje gate'ów)"
date: 2026-09-14
type: phase-method-decisions
tgp_owner: research/op-r3-stationary-states-2026-09-14
status: FROZEN
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-action-audit-spectrum-insert-2026-09-13/Phase1_output.txt]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase3_relax_M911.py]]"
---

# Phase_method_decisions (FROZEN przed jakimkolwiek kodem cyklu)

**Status: FROZEN 2026-09-14, ZERO obliczeń cyklu przed zapisem tego
dokumentu.** Wszystkie decyzje poniżej są zamrożone; zmiana po
pierwszym biegu = forbidden move (LOCK §3 b,c). Pozycje wykraczające
poza literę LOCKa oznaczone **[INPUT-MD]**.

---

## 1. Formy modelu (CYTATY — dziedziczone bez modyfikacji; LOCK §1, forbidden move a)

Źródło: `../op-action-audit-spectrum-insert-2026-09-13/Phase1_output.txt`
(Q-D1-PASS; formy są tam WYNIKIEM wyprowadzenia z jednej akcji
[eq:S-TGP-unified-M911-canonical] z metryką [eq:metric-M911-canonical],
człon czasowy w konwencji |g^tt| — decyzja user-gate N2 2026-09-14,
dopiski core rem:W-sign-axiomatic(iv) + rem:psi-EOM-R3-branch-status,
sek08a — TYLKO ODCZYT):

- „M(psi) = sqrt(-g)*K*|g^tt| = psi**6/(3*psi - 4)**2" ⟹
  **M(ψ) = K_geo ψ⁶/(c₀(4−3ψ)²)**; „M'(psi) = 12*psi**5*(psi-2)/(3*psi-4)**3"
  ⟹ **M′(ψ) = 12ψ⁵(2−ψ)/(4−3ψ)³**
- „Keff(psi) = psi**4" ⟹ **𝒦(ψ) = K_geo ψ⁴**, 𝒦′ = 4ψ³ (odczyt B DZIEDZICZONY)
- „Ueff = w*V = psi**4/4 - psi**3/3" ⟹ **𝒰(ψ) = γ(ψ⁴/4 − ψ³/3)**;
  𝒰′ = γψ²(ψ−1); 𝒰″ = γψ(3ψ−2); 𝒰‴ = γ(6ψ−2); 𝒰⁗ = 6γ
- „pi = dL/dpsidot = M(psi) psidot"; „EOM: M psiddot + 1/2 M' psidot^2
  = div(Keff grad psi) - 1/2 Keff' |grad psi|^2 - Ueff'"
- **K_geo = γ = c₀ = 1 [LOCK §1]**. Dziedzina ψ∈(0,4/3); pas
  ψ ≥ 4/3−1e−6 lub ψ ≤ 1e−6 = BREAKDOWN-BOUNDARY (klasyfikacja,
  stop); ZERO podłóg/barier.

Radialnie 3D (LOCK §1): Mψ̈ + ½M′ψ̇² = (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′.

## 2. π-formulacja (JAWNIE; HANDOFF) i dyskretyzacja przestrzenna

Para kanoniczna: gęstość pędu **π = M(ψ)ψ̇**. Równania Hamiltona
(gęstość H = π²/(2M) + ½𝒦|∇ψ|² + 𝒰):

- **ψ̇ = π/M(ψ)**
- **π̇ = π²M′(ψ)/(2M(ψ)²) + (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′ − γ_sp(r)π**

(równoważność z EOM: π̇ = Mψ̈ + M′ψ̇² ⟹ Mψ̈ + ½M′ψ̇² = RHS_sp;
wkład ½M′ψ̇² wchodzi przez człon +π²M′/2M²; tożsamość zweryfikowana
sympy w Phase 1, deskryptywnie. γ_sp — sponge, §4, tylko r≥160.)

**Siatka radialna (adaptacja z `op-metric-pair-M911/Phase3_relax_M911.py`,
tam gradient flow — silnik czasowy tu nowy):** r_i=(i+½)h, i=0..N−1,
N=R/h; regularność w r=0 przez zerowy strumień (siatka przesunięta,
ψ′(0)=0 realizowane brakiem strumienia przez r=0); brzeg zewnętrzny
r=R: zerowy strumień (Neumann) — dopuszczalne, bo obiekt żyje w r≲40,
a warstwa absorpcyjna [160,200] tłumi (LOCK §1).

**Dyskretna energia (jednostki z czynnikiem 4π; próżnia ODJĘTA —
E[ψ≡1,π≡0]=0, gęstość ≥0 bo 𝒰 ma minimum globalne w ψ=1 na (0,4/3)):**

E = 4π[ Σ_i h r_i²( π_i²/(2M(ψ_i)) + 𝒰(ψ_i) − 𝒰(1) )
      + ½ Σ_{i=0}^{N−2} h r_{i+½}² 𝒦(ψ̄_{i+½}) ((ψ_{i+1}−ψ_i)/h)² ],
ψ̄_{i+½} = (ψ_i+ψ_{i+1})/2, r_{i+½} = (i+1)h.

**Siła przestrzenna dyskretna (wariacyjna, struktura flux+quad
identyczna jak w FlowRadial.rhs poprzednika):**
F_i(ψ) = −(1/(h r_i²)) ∂E_sp/∂ψ_i, gdzie E_sp = człony gradientowe
+ potencjalne. Gwarantuje spójność siły z dyskretną energią (gate
dryfu mierzy dokładnie tę parę).

**E_core (LOCK §1):** suma powyższych członów dla r_i ≤ 80 (człony
gradientowe: r_{i+½} ≤ 80). **ψ(0,t) ≔ ψ(r₀=h/2, t)** (proxy centrum,
błąd O(h²) przy ψ′(0)=0) **[INPUT-MD]**.

## 3. Integrator (FROZEN; LOCK §1: 2. rząd, para (ψ,π))

H nieseparowalny (π²/2M(ψ)) ⟹ **uogólniony Störmer–Verlet**
(symetryczny, 2. rzędu; Hairer–Lubich–Wanner GNI, rozdz. II —
kompozycja dwóch sprzężonych metod Eulera symplektycznego; wybór
uzasadniony SYMETRIĄ czasową ⟹ brak sekularnego dryfu energii):

1. π^{n+½} = π^n + (dt/2)[F(ψ^n) + G(ψ^n, π^{n+½})]  (implicit w π)
2. ψ^{n+1} = ψ^n + (dt/2) π^{n+½}[1/M(ψ^n) + 1/M(ψ^{n+1})]  (implicit w ψ)
3. π^{n+1} = π^{n+½} + (dt/2)[F(ψ^{n+1}) + G(ψ^{n+1}, π^{n+½})]  (explicit)

gdzie F(ψ) = siła przestrzenna dyskretna (§2) − 𝒰′ już w F,
G(ψ,π) = π²M′/(2M²) − γ_sp π. Kroki implicit: iteracja punktu
stałego, tolerancja ‖Δ‖∞ ≤ 1e−14·max(1,‖·‖∞), max 100 iteracji
(brak zbieżności ⟹ krok liczony jak niefinityczny → BREAKDOWN)
**[INPUT-MD]**. dt=0.005 [LOCK]; kontrola dt/2 przy zdarzeniach (§6).

Uwaga deskryptywna (nie zmienia LOCKa): schemat jawny przestrzennie —
lokalna prędkość c(ψ)=(4−3ψ)/ψ rośnie przy ψ→0; przekroczenie CFL
przy głębokich wychyleniach objawia się niefinitycznością ⟹
klasyfikacja **BREAKDOWN** (deskryptywna, jak u poprzednika M911),
z kontrolą dt/2.

## 4. Sponge (FROZEN; LOCK §1: gładka, r∈[160,200])

γ_sp(r) = γ₀·S(x), x = (r−160)/40 obcięte do [0,1],
S(x) = 6x⁵−15x⁴+10x³ (smootherstep, C²; γ_sp(160)=0 gładko,
γ_sp(200)=γ₀). **γ₀ = 1.0 [INPUT-MD]**. Sponge działa wyłącznie na π
(człon −γ_sp π w π̇). Zakaz modyfikacji po pierwszym biegu (LOCK §3b).
Sponge ON we wszystkich biegach produkcyjnych Phase 3 i w P2a;
OFF w testach zamkniętego pudła (P2b dyspersja, P2c dryf energii —
pomiar zachowania energii wymaga układu zamkniętego) **[INPUT-MD]**.

## 5. Detektor oscylonu (FROZEN; LOCK §1) i klasyfikacja biegów

Próbkowanie: ψ(0,t), E_core(t) co **dt_out=0.1**; pełny profil co 50
j.cz.; T₀=2π; t_transient=50; E_ref ≔ E_core(50).

- **Okno podtrzymania:** t_hold = największe t takie, że
  E_core(t′) ≥ 0.5·E_ref dla WSZYSTKICH próbek t′∈[50,t].
- **KANDYDAT:** (t_hold−50) ≥ 100·T₀ = 628.318... ORAZ liczba przejść
  ψ(0,t) przez 1 (zmiany znaku ψ(0,t)−1 po próbkach) na [50,t_hold]
  ≥ 50.
- **Czas życia:** τ = t_hold − 50; jeśli t_hold = t_max=1000 ⟹ τ
  cenzurowane („≥950").
- **POTWIERDZENIE (OSCILLON):** kandydat na h=0.05 ORAZ h=0.025 ORAZ
  w biegu dt/2 (dt=0.0025, na h=0.05 **[INPUT-MD]**); zgodność czasów
  życia: |τ_a−τ_b|/max ≤ 10% dla par (h,h/2) i (dt,dt/2); jeśli
  jedna wartość cenzurowana — zgodność ⟺ druga ≥ 0.9·950=855;
  obie cenzurowane ⟹ zgodne **[INPUT-MD]**.
- **RADIATED:** kandydat NIE wystąpił ORAZ E_core(t_max) ≤ 0.05·E_ref
  **[INPUT-MD — operacjonalizacja „E_core→0 przed progiem"]**;
  zbieżnie = obie siatki dają RADIATED.
- **BREAKDOWN-BOUNDARY:** max ψ ≥ 4/3−1e−6 lub min ψ ≤ 1e−6
  (stop, klasyfikacja; kontrola dt/2). **BREAKDOWN:** niefinityczność
  (deskryptywna). **INCONCLUSIVE:** cała reszta (w tym rozjazd siatek).
- **Zdarzenie (⟹ kontrola dt/2):** kandydat na dowolnej siatce lub
  BREAKDOWN* na dowolnej siatce **[INPUT-MD]**.
- **Pomiar ω (kandydaci):** FFT ψ(0,t)−⟨ψ(0,t)⟩ z oknem Hanna na
  odcinku stabilnym t∈[max(300, t_hold−600), t_hold], interpolacja
  paraboliczna piku; raport pik dominujący + harmoniki (piki o mocy
  ≥1% dominującego) **[INPUT-MD — okno]**.

## 6. Starty (LOCK §1, deterministyczne; π₀≡0)

(i) ψ = 1 + a·exp(−r²/2σ²), a∈{−0.3,−0.15,+0.15,+0.25}, σ∈{3,6} — 8;
(ii) ψ = 1 + a·(sin r)/r·exp(−r²/2σ_w²), a∈{−0.2,+0.2}, σ_w=15 — 2
(sin(r)/r przez np.sinc(r/π), regularne w r=0);
(iii) kontrola: ψ≡1 (zero alarmów detektora).
Siatki h∈{0.05,0.025}; R=200; dt=0.005; t_max=1000.

## 7. Definicje gate'ów Phase 2 (FROZEN)

- **P2a (próżnia):** konfiguracja produkcyjna (sponge ON), ψ≡1, π≡0,
  100·T₀, obie siatki; gate: ‖ψ−1‖∞ ≤ 1e−10 przez cały bieg; detektor
  zero alarmów.
- **P2b (dyspersja):** pudło zamknięte (sponge OFF), R=200, h=0.05
  (primary; h=0.025 deskryptywnie), start ψ = 1+1e−3·exp(−r²/(2·2²)),
  π₀=0, t=400, zapis u(r,t)=r·(ψ−1) co dt_out=0.2; FFT 2D z oknem
  Hanna (r i t); dla binów k najbliższych {0.6, 1.0, 1.4} (gate; biny
  {1.8, 2.2} deskryptywnie) pik ω>0 z interpolacją paraboliczną;
  **gate: |ω_meas−√(k_bin²+1)|/√(k_bin²+1) ≤ 1% na WSZYSTKICH
  3 binach gate'owych** (litera LOCKa „≥3 modach") **[INPUT-MD —
  wybór binów i pulsu]**.
- **P2c-energia:** pudło zamknięte (sponge OFF), obie siatki, start
  ψ = 1+1e−3·exp(−r²/(2·3²)), t=100·T₀; E(t) pełne (próżnia odjęta);
  **dryf ≔ |⟨E⟩_{t∈[90T₀,100T₀]} − ⟨E⟩_{t∈[0,10T₀]}| / ⟨E⟩_{[0,10T₀]}
  ≤ 1e−6** (średnie po oknach eliminują ograniczoną oscylację energii
  właściwą metodom symetrycznym; dryf = zmiana sekularna — to ona jest
  przedmiotem gate'u LOCKa) **[INPUT-MD — operacjonalizacja]**;
  deskryptywnie także max|E−E(0)|/E(0).
- **P2c-sponge (odbicie):** bieg A: R=200 sponge ON; bieg B
  (referencja): R=400 sponge OFF (ściana w 400 nie zdąży odbić do
  r≤120 przed t=250); start wspólny ψ = 1+1e−3·exp(−(r−100)²/(2·5²));
  h=0.05, t=250; w zmiennej u=r(ψ−1):
  **gate: max_{r≤120,t≤250}|u_A−u_B| / max_{r∈[120,160],t≤250}|u_B|
  ≤ 1e−3** (różnica biegów = czysty artefakt sponge; mianownik =
  amplituda padająca na wejściu warstwy) **[INPUT-MD]**.
- FAIL któregokolwiek ⟹ STOP (LOCK §2).

## 8. Phase 1 (sympy; zero ewolucji)

- P1a: linearyzacja stacjonarna ψ=1+e^{−iωt}f(r) w pełnym EOM ⟹
  ∇²f = −κ²f, κ² = (M(1)ω²−𝒰″(1))/𝒦(1); wyprowadzenie sympy, klasy
  ω<1 / ω>1; mapowanie na linearyzację R3: κ=1 ⟺ ω²=2.
- P1b: Lindstedt–Poincaré O(a²) dla modu jednorodnego rdzenia
  (M(ψ)ψ̈+½M′ψ̇²+𝒰′=0, ψ=1+u, u=a·u₁+a²u₂+a³u₃, ω=1+a²ω₂; warunek
  braku sekularności w O(a³)); ω₂ symbolicznie (wchodzą 𝒰‴(1),𝒰⁗(1),
  M′(1) i — jeśli wyjdzie z rachunku — M″(1); raport form);
  **PREDYKCJA pre-rejestrowana: znak ω₂ (ω₂<0 ⟺ miękka nieliniowość
  ⟺ warunek istnienia oscylonu małej amplitudy ω(a)<m=1)** —
  zakaz reinterpretacji po Phase 3 (LOCK §3h).
- P1c (gate): M,M′,𝒦,𝒦′,𝒰,𝒰′,𝒰″,𝒰‴,𝒰⁗ w ψ∈{0.9,1,1.1}: sympy
  (na dokładnej binarnej reprezentacji double wejścia — lekcja
  correction note poprzednika) vs float; próg
  |Δ| ≤ 1e−12·max(1,|wartość|) **[INPUT-MD — normalizacja]**.

## 9. Phase 4 (warunkowo przy Q-E-PASS): rodziny (n, ω, τ)

- **Obwiednia:** A(r) = ⟨|ψ(r,t)−1|⟩_t po oknie stabilnym (okno FFT
  §5), r ≤ 60; wygładzenie średnią ruchomą o szerokości 0.45 j.dł.
  (9 pkt przy h=0.05, 17 przy h=0.025) **[INPUT-MD]**.
- **Liczba węzłów n:** liczba minimów lokalnych A(r) w (0,60)
  spełniających A(r_min) < 0.3·min(sąsiednie maksima lokalne)
  **[INPUT-MD — próg prominencji]**.
- Rodziny: grupowanie po n; w rodzinie raport ω (średnia i rozrzut),
  τ. Werdykt Q-F wg litery LOCKa §2 Phase 4 (próg czasu życia rodzin
  n≥3 = ten sam próg detektora 100·T₀).
- Klasa częstości: raport, czy ω<1 (rdzeń związany) czy dyskretne
  ω>1 (rezonanse) — litera Q-F(b).

## 10. Rejestr WEJŚĆ (flagowane)

[LOCK]: K_geo=γ=c₀=1; rodziny startów (i)/(ii)/(iii) §6; h∈{0.05,
0.025}; R=200; sponge r∈[160,200]; dt=0.005; t_max=1000; t_transient
=50; próg detektora 0.5·E_ref; 100·T₀; ≥50 przejść; E_core r≤80;
zbieżność ≤10%; pas 4/3−1e−6 i 1e−6; dryf ≤1e−6/100T₀; odbicie
≤1e−3; dyspersja ≤1% ≥3 mody; próżnia 1e−10; brak seeda.
[INPUT-MD]: γ₀=1.0 + smootherstep; sponge OFF w P2b/P2c; dt_out=0.1
(P2b: 0.2); puls dyspersyjny a=1e−3 σ=2; puls energetyczny a=1e−3
σ=3; puls odbiciowy a=1e−3 σ=5 @ r=100; biny gate'owe k≈{0.6,1.0,1.4};
tol. punktu stałego 1e−14 (max 100 iter.); RADIATED: E_core(t_max)
≤0.05·E_ref; reguła cenzurowania τ; dt/2 na h=0.05; okno FFT
[max(300,t_hold−600), t_hold]; obwiednia r≤60, wygładzenie 0.45,
prominencja 0.3; ψ(0,t)≔ψ(h/2,t); profil co 50 j.cz.; energia
z czynnikiem 4π i odjętą próżnią.

**FROZEN. Zmiany poniżej po starcie obliczeń = forbidden move.**
