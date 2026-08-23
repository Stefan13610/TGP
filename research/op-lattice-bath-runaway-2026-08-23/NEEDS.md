---
title: "NEEDS — op-lattice-bath-runaway (po STOP bramki A na A3)"
date: 2026-08-23
type: needs
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[PhaseA_report.md]]"
  - "[[../../core/formalizm/dodatekH_lancuch_wyprowadzen.tex]]"
  - "[[../../core/_meta_latex/status_map.tex]]"
---

# NEEDS (user-gated)

Bramka A padła na A3 (fazy ogona niereprodukowalne). Propozycje dalszej
drogi, do decyzji użytkownika:

> **UPDATE 2026-08-23 (autoryzacja użytkownika: „zajmij się N1 i N2"):**
> - **N1 → 🟢 WYKONANE (śledztwo):** fazy ZNALEZIONE i zreprodukowane co do
>   setnej stopnia — pochodzą z `_archiwum/.../p131_eta_refinement.py`,
>   układ eq:J-ode (F=1+2α_eff·ln g, Formulacja B), NIE układ potęgowy O-L5.
>   Pełna analiza: [[ANALIZA_N1_pochodzenie-faz_2026-08-23.md]]
>   (+ [[N1_provenance_check.py]]). Korekty rdzenia: nadal USER-GATED.
> - **N2 → 🟢 WYKONANE (dokumentacyjnie):** źródło rozdwojenia znaku W
>   zlokalizowane W SAMYM sek08a (prop:field-eq-from-action vs inline
>   „Nota kanoniczna 2026-04-07", której EL = równanie maszynerii 2);
>   L04 canonicalization rozstrzygnął tylko kinetykę, nigdy znak potencjału.
>   Hipoteza autora (pojedynczy obiekt vs kąpiel) ma dokładny odpowiednik
>   strukturalny — bez derywacji (kandydat do N3). Derywacja jednego znaku
>   z akcji: OTWARTA. Pełna analiza: [[ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]].
>   Korekty rdzenia: nadal USER-GATED.
> - N3, N4: bez zmian (poniżej).

## N1 — Ustalenie pochodzenia faz w dodatekH lin. 1126–1129 (PRIORYTET)

Wartości δ_e=−81.4°, δ_μ=+38.6°, Δ(e→μ)=120.01°≈2π/3 (oraz mapa fazowa
p127–p128: δ_e=−81.1°, δ_μ=+43.8°, A_μ=0.3861 przy g₀=1.4550) nie są
reprodukowalne z udokumentowanego ODE w żadnym z 6 zbadanych wariantów,
a A_μ=0.3861 jest sprzeczne z własnymi skryptami rdzenia
(atail_functional: A(1.455)≈0.594). Potrzebne: odszukanie SKRYPTU/
rachunku, który wygenerował te liczby (być może inna kalibracja g₀,
inna definicja fazy, inny wariant biegnącego α, albo błąd przepisania
do .tex). Do rozstrzygnięcia u źródła zanim „Δ=2π/3" będzie dalej
cytowane. Dotyczy statusu O-K1 („Δ(e→μ)=120° — potwierdzone") —
proponowany dopisek Limitations/flagi [UNVERIFIED-SOURCE] przy tych
liniach (zmiana w core/.tex — wyłącznie za zgodą użytkownika).

## N2 — Rozstrzygnięcie znaku W z akcji TGP (rysa z A5)

A5 zlokalizowało różnicę między maszynerią 2 a „F-A kanoniczną"
z AUDYT_KRYTYCZNY Dod. A: znak członu źródłowego (W̃=g⁷/7−g⁸/8 vs
W=u⁸/8−u⁷/7). Rdzeń definiuje równanie (źródło g²(1−g) po podzieleniu
przez K), ale nie wyprowadza jawnie W z akcji TGP; EL akcji F-S
(W_FS=g³/3−g⁴/4) z K=g⁴ daje trzecie, jeszcze inne równanie
((1−g)/g²). Potrzebna jednoznaczna derywacja: które W wynika z akcji
(L04 canonicalization) — bo od znaku W″(1) zależy sama oscylacyjność
ogona (fundament locku oscylacyjnego z retrospektywy 2026-08-16).

## N3 — Decyzja: Phase 1–2 na fazach ZMIERZONYCH zamiast dokumentowanych

Rachunek centralny (runaway w kąpieli, ω²(n)) pozostaje niepoliczony,
a jego składniki, które przeżyły audyt (ω=1, oscylacyjność, A≈|g₀−1|,
próg 8/5), wystarczają technicznie do jego wykonania z fazami
ZMIERZONYMI w A3 (δ_e=−75.50°, δ_μ=+88.48°, Δ=163.98° — z faktycznego
ODE) zamiast z dokumentacji rdzenia. Wymaga NOWEGO cyklu z nowym
Phase0-lockiem (ten cykl pozostaje CLOSED-GATE-STOP; reguła bramki nie
pozwala kontynuować pod obecnym lockiem). Uwaga: wartość poznawcza
takiego rachunku nie zależy od tego, czy Δ=120° — drabina d* istnieje
dla dowolnych faz.

## N4 — Rysa τ przy g₀=4 (biegnące α_eff)

W niezależnej integracji τ (g₀=4, η_K=181/15) KOLAPSUJE; g₀=4 leży
na granicy istnienia (lim α→0 progu = 4). Spójne z p134e−g (A_τ
monotoniczna, kres skanu), ale wymaga flagi przy „g₀^τ=4 wyprowadzone
z K-samodzielności": obiekt może nie istnieć jako regularne
rozwiązanie. Do weryfikacji u źródła przy okazji N1.
