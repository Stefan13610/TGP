# Plan P6 — broader equations dla cuprates, hydrydów i f-metali

**Status:** roboczy; powstał po zakończeniu P5 (ps1–ps8).

## Motywacja — czego nie pokazuje P5

Global fit ps5 dał r(log-log)=0.48 na 19 znanych SC, co jest **zbyt słabe** dla predykcji. Krytyczne outliery (|Δlog₁₀(T_c)| > 0.5):

| Klasa | Przykład | T_obs | T_pred (ps5) | off | Mechanizm problemu |
|-------|----------|-------|--------------|-----|--------------------|
| hydryd | H3S | 203 K | 5.5 K | -1.57 | coupling fonon H-substrat nie ujęty |
| hydryd | LaH10 | 250 K | 250 K* | 0 | hit bounds (sztuczne) |
| cuprate | YBCO | 92 K | 15 K | -0.79 | 2D d-wave, CuO₂ plane, nie 3D soliton lattice |
| cuprate | BiSCCO | 110 K | 15 K | -0.86 | idem |
| dwu-przerwowy | MgB2 | 39 K | 8 K | -0.68 | σ+π dwa kanały |
| f-metal (ps8) | Ce, Yb | <5 K | 1300 K | +2.5 | A_f P4 nie mapuje 4f w metalu |

## Trzy osie poszerzenia

### P6.A — Dwu-przerwowe, anizotropia, d-wave

**Cel:** zamienić skalar `T_c = k_d * J * M` na **tensor**

$$
\hat{T}_c = k_d(z) \cdot \sum_{\alpha,\beta} J_{\alpha\beta}(a, g_0) \cdot M_{\alpha\beta}(\text{struktura}) \cdot \Lambda_E
$$

gdzie α,β ∈ {σ, π, d} indeksują kanały parujące.

**Konkretnie dla MgB2:**
- σ-band: A_sp + A_mu (B-B σ-bond, silne sprzeżenie)
- π-band: A_sp (B-B π-bond, słabsze)
- T_c(MgB2) = max(T_c^σ, T_c^π) lub średnia ważona.

**Dla cuprates (YBCO, BiSCCO):**
- Model 2D (SC→BCC redukcja z=12→4 dla 2D plakaty) + d-wave znak.
- k_d^2D(z=4) = 0.893 (BKT threshold), ale z anomalnym spikeprzy progu.
- Modyfikacja: `M_dwave = cos²(kx a) - cos²(ky a)` — faza topo.

**Testy:**
- ps9: formula dwu-przerwowa dla MgB2; porównaj z danymi d¹σ=1.8 meV, d¹π=7.2 meV.
- ps10: 2D soliton-lattice fit dla YBCO/BiSCCO/Hg-1223 (trzy cuprates).

### P6.B — Hydrydy i coupling fonon-substrat

**Kluczowa fizyka:** H-atomy mają niskie masy → wysokie ω_Debye → mocne coupling do oscylacji substratu Φ.

**Hipoteza:**
$$
\Lambda_E^{(eff)}(\text{hydryd}) = \Lambda_E^{(0)} \cdot \left(\frac{\omega_H}{\omega_0}\right)^\alpha
$$

gdzie ω_H ≈ 200 meV (H-optical mode w H3S/LaH10), ω_0 ≈ 15 meV (fonon Debye dla Al), α do wyznaczenia.

**Dla H3S:** ω_H/ω_0 ≈ 13; jeśli α≈1.5, to Λ_E wzrasta ~50× → T_c(H3S) ~5.5·50 = 275 K ✓ (obs 203 K)

**Dla LaH10:** ω_H/ω_0 ≈ 14; α≈1.5 → T_c(LaH10) ~ 85·2 ~ 170 K (obs 250 K) — wciąż trochę mało, ale rząd OK.

**Testy:**
- ps11: fit α na 5 hydrydach (H3S, LaH10, YH9, CeH10, CaH6).
- ps12: coupling factor = f(ω_H, m_H/m_substrate) — wyprowadzenie z TGP.

### P6.C — f-metale i orbital switching pod cisnieniem

**Problem:** ps8 dla Ce, Yb daje T_c~1300 K przy ambient — absurd. A_f=2.03 z fitu 5c jest zbyt duże dla 4f w metalu.

**Hipoteza:** elektrony 4f w metalu przechodzą w hybrydyzację z pasmem przewodnictwa → efektywnie stają się **d-like** (nie τ-soliton).

Formalnie:
$$
A_{eff}(\text{f-metal}) = (1 - \eta) A_f + \eta A_d,\quad \eta(\rho, P) \to 1 \text{ pod cisnieniem}
$$

gdzie η(ρ, P) to parametr hybrydyzacji.

**Testy:**
- ps13: fit η dla szeregu f-metali (Ce, Yb, Sm, Eu, Gd...).
- ps14: orbital-switching pod cisnieniem (Ce: γ→α faza ~0.7 GPa).

## Co pozostaje w domenie P5

Wariant 5c z current parametrami dobrze opisuje:
- BCS klasyczne (|Δlog|<0.3): Al, Sn, Hg, Pb, Zn, In, Nb, V, Ta, Tc
- A15 z hybrydą łańcuchową (ps7 RMS=0.64 vs 0.86 Bravais): V3Si, Nb3Sn, Nb3Ge, Nb3Al
- niskie-T_c pressure-driven: Ba, K, Cs, Sr

P5 w wersji publikowanej powinien jasno zaznaczyć że dla cuprates/hydrydów/f-metali trzeba rozszerzenia P6.

## Zmienne globalne do kalibracji P6

| Nazwa | ps5 wartość | P6 będzie | 
|-------|-------------|-----------|
| C_0 | 48.8222 | = (przenośne) |
| a*_TGP | 4.088 A | = (przenośne) |
| Λ_E | 0.131 meV | zmienia się per klasa (P6.B) |
| σ_a | 2.59 A | = (przenośne) |
| A_s, A_sp, A_d, A_f | fit per variant | tensorowe A_{αβ} (P6.A) |
| η_orbital | — | nowy (P6.C) |
| α_phonon | — | nowy (P6.B) |

## Kolejność prac

1. ~~ps6-ps8~~ done.
2. **ps9 (MgB2 dwu-przerwowy)** — 1 parametr, 1 materiał, bardzo konkretny test.
3. **ps10 (cuprates 2D)** — redukcja wymiaru; test na 3-5 cuprates.
4. **ps11 (hydrydy coupling α)** — kalibracja na H3S/LaH10, predykcja nowych.
5. **ps13 (f-metale η)** — orbital hybridization.
6. Sformalizowanie w dodatek `dodatekP6_broader_sc_mechanisms.tex`.

## Ryzyka

1. **Przedopasowanie**: dodając kanały, wagi, coupling, łatwo trafić każdy punkt — trzeba trzymać parametr/materiał < 0.5.
2. **Podwójne liczenie**: P5 już ma A_orb, P6.A chce tensor J_αβ — nie można jednoczesnie mieć obu, inaczej degeneracja.
3. **Domena stosowalności**: P6.B zakłada ω_H >> ω_0. Dla organic SC (TMTSF₂X) ω jest inne — może inny reżim.
