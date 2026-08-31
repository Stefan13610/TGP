# NEEDS — op-substrate-fluctuation-channel-2026-08-23 (user-gated)

> **LOG WYKONANIA (2026-08-31):**
> - **N3 ✅ EXECUTED** (user-gate krok 2) — remark
>   `rem:fluctuation-channel-bridge` w dodatekB (przy cor:entropy-potential):
>   inwentarz kanałów poziomu 0 + połączenie m_eff²=γ(1+T_Γ),
>   γ≈12Λ_eff/Φ₀ ⟹ reżim efektywnie krytyczny wewnątrz horyzontu.
>   Build gate PASS (build_gate_2026-08-31.log, zero nowych undefined refs).
> - **N5 ✅ EXECUTED** (user-gate krok 2) — POST-SCRIPTUM w ANALIZA_N2
>   (op-lattice-bath-runaway) + pkt (iii) remarku `rem:W-sign-axiomatic`
>   w sek08a (rozrzedzenie vs zagęszczenie; odwrotna rola gęstości).
> - **N1 + N2 → realizowane 2026-08-31** jako cykl-następca
>   `op-fluctuation-extended-nbody-2026-08-31` (własny LOCK Phase 0;
>   autoryzacja: user 2026-08-31 „zajmij się
>   research/op-substrate-fluctuation-channel-2026-08-23/").
> - **N4 — OPEN** (bez zmian; czeka na decyzję).

Rdzeń .tex: dopiski addytywne wykonane wyłącznie w ramach user-gate wyżej.

## N1 (research): kanał fluktuacyjny dla obiektów ROZCIĄGŁYCH
Defekt punktowy daje na krytyczności potencjał ~−1/d² (nie Newtonowskie
−1/d). Dla obiektów rozciągłych (sfera zamrożonego pola o promieniu R,
model rdzenia solitonu) wykładnik i prefaktor się zmieniają — policzyć
F_fl(d; R) (macierz kowariancji wielopunktowa, ta sama maszyneria co
Phase 2, det większej macierzy). Pytanie binarne: czy istnieje reżim
R/d, w którym kanał fluktuacyjny daje −1/d.

## N2 (research): N-ciałowość kanału fluktuacyjnego
F_fl dla trzech+ defektów NIE jest sumą par (log-det nie jest addytywny) —
zbadać odchylenie od addytywności parowej vs wynik op-nbody-additivity
(kanał klasyczny: addytywny z poprawką ~exp). Test spójności z bilansem
programu „most do grawitacji".

## N3 (core, user-gate): remark łączący w dodatekB
Dopisek (addytywny) przy prop:continuum-conditions / cor:entropy-potential:
ξ = 1/√γ z γ ≈ 12Λ_eff/Φ₀ ⟹ reżim efektywnie krytyczny wewnątrz horyzontu;
inwentarz kanałów poziomu 0 (ten cykl): jedyny kanał o uniwersalnym znaku
przyciągania = fluktuacyjny, zasięg 2m vs m, krytycznie −1/d² vs −1/d.

## N4 (research): QB poza MFT
Krzywizna w kąpieli z pełną całką jednowęzłową (dokładna 1D kwadratura
zamiast MFT; opcjonalnie mały MC) — czy próg Φ_c/Φ_vac ≈ 0.20–0.33
przeżywa fluktuacje termiczne on-site.

## N5 (doc, user-gate): nota porównawcza do ANALIZA_N2 (op-lattice-bath-runaway)
Dopisek: na poziomie 0 MFT znak tachionowy pochodzi z ROZRZEDZENIA
(spinodala on-site; wiązanie gradientowe zawsze stabilizuje, +8zJs_b⁶),
z odwrotną rolą gęstości niż w hipotezie dwóch sektorów — do uwzględnienia
przy interpretacji wyników op-bath-two-sectors Q2.
