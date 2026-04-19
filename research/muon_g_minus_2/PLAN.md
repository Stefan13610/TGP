# Anomalny moment magnetyczny mionu $(g-2)_\mu$ z TGP

**Data startu:** 2026-04-20
**Status:** plan / open  **(Flagship falsyfikowalne zadanie)**
**Kategoria:** TGP applications → leptony → poprawki radiacyjne

## 1. Problem fizyczny

Anomalny moment magnetyczny $a_\mu = (g-2)/2$ mionu jest jednym z **najprecyzyjniej
mierzonych i najprecyzyjniej obliczanych** obserwabli w fizyce.

**Eksperyment:**
- BNL E821 (2006): $a_\mu = 116592080(63)\times 10^{-11}$
- Fermilab Run-1+2+3 (2021-2023): $a_\mu = 116592055(24)\times 10^{-11}$
- World average: $a_\mu^\text{exp} = 116592059(22)\times 10^{-11}$

**Teoria (Standard Model):**
- White Paper 2020 (Aoyama et al.): $a_\mu^\text{SM} = 116591810(43)\times 10^{-11}$
  → **rozbieżność $\approx 4.2\sigma$, $\Delta a_\mu \approx 249\times 10^{-11}$**
- Jednak BMW-lattice 2021 i następne 2024-2025 obniżają HVP kontrybucję o ~100×10⁻¹¹,
  **zmniejszając rozbieżność do ~1-2σ**. Jury still out.

**Obecny status (kwiecień 2026):** rozbieżność zawęziła się do $(100\pm 50)\times 10^{-11}$
w zależności od wyboru estymacji hadronowej vacuum polarisation (HVP) — data-driven
(e+e-→hadrons) vs lattice QCD. Ostateczna precyzja oczekiwana po pełnym Fermilab
Run-6 (2026-2027).

## 2. Dlaczego TGP

**Kluczowa hipoteza:** w TGP metryka efektywna
$g_{ij} = e^{+2U/c_0^2}\delta_{ij}$
modyfikuje **propagatory wirtualnych cząstek** w pętlach Feynmana mionu.
Poprawka może wejść w:

- **Schwinger diagram** (1-loop QED, $\alpha/2\pi$) — wewnętrzna pętla fotonu
  w TGP-metryce daje dodatkowy czynnik $\langle e^{+2U/c^2}\rangle\approx 1 + O(U/c^2)$.
- **Hadronic vacuum polarisation** — pętla kwarkowa "czuje" lokalny substrat;
  ten sam mechanizm może wyjaśnić różnicę e+e- vs lattice.
- **"Geometric self-energy"** — poprawka od nieliniowości potencjału
  $V(\Phi) = \beta\Phi^3/3\PhiZero - \gamma\Phi^4/4\PhiZero^2$ sprzęgnięta
  do $\mu$-spin-current.

**Oczekiwany rząd wielkości TGP poprawki:**
$$\Delta a_\mu^\text{TGP} \sim \frac{\alpha}{2\pi}\cdot \frac{m_\mu^2 c_0^2}{\Lambda_E^2}
  \cdot f(\beta,\gamma,\text{topology})$$
gdzie $\Lambda_E$ to skala TGP (mikrofala meV w SC, ale efektywnie inaczej
w QED mionu). Jeśli $f \sim 1$ i $\Lambda_E$ = właściwa skala, można dostać
$\Delta a_\mu^\text{TGP} \sim 100\times 10^{-11}$.

## 3. Cele badawcze

### G1 — 1-loop QED z TGP metryką

Przerobić rachunek Schwingera 1948 w metryce $e^{+2U/c_0^2}\delta_{ij}$
(+$ e^{-2U/c_0^2}dt^2$). W granicy $U\to 0$ odzyskać $\alpha/2\pi$.
W skończonym $U$ (lokalny potencjał efektywnej TGP-pola indukowanej
przez samą masę mionu) przewidywać korektę.

Kluczowa trudność: jaki to jest $U$? Propozycje:
- **Self-consistent**: $U = G m_\mu/r$ w promieniu Comptona $\hbar/m_\mu c$ —
  daje poprawkę $\sim G m_\mu^2/\hbar c \sim 10^{-42}$, zbyt mała.
- **Substrate-local**: skala w której $\Phi$ deformuje się od mionu to nie
  Compton grawitacyjny, lecz **TGP-Compton** $\hbar c/\Lambda_E^\text{TGP}$
  — jeśli $\Lambda_E^\text{TGP} \sim $ MeV lub GeV → $\Delta a_\mu \sim 10^{-9}$.

### G2 — HVP reinterpretacja w TGP

Jeśli TGP modyfikuje propagator kwarkowy w pętli HVP, to rozbieżność
e+e-→hadrons (data-driven) vs lattice może być **efektem środowiskowym**:
e+e- dzieje się w próżni laboratoryjnej (nisko-$\Phi$), lattice jest
obliczana w "idealnym" wakuumie (też nisko-$\Phi$, ale bez energetycznego
tła). Różnica znika na poziomie partonu, ale może ujawnić się na poziomie
hadronowym. **Testowalne** bo obie metody powinny dać ten sam wynik.

### G3 — Light-by-light (HLbL) i TGP

Aoyama et al. HLbL contribution: $92(18)\times 10^{-11}$. Pętla 4-fotonowa
jest czuła na nieliniowość propagatora — idealne miejsce na test nieliniowości
$\Phi^3, \Phi^4$ TGP.

### G4 — Elektron $(g-2)_e$ jako "zero-point"

$(g-2)_e$ jest zmierzone do $10^{-13}$ i zgadza się z SM w $2.5\sigma$
(nieużywając $\alpha$ z Cs), ale używając $\alpha$ z Cs/Rb różni się o $5\sigma$.
Ten problem JEST niezależny od mionu. TGP musi dać spójny wynik:
$a_e, a_\mu, a_\tau$ skalują się z $m_\ell^2/\Lambda_E^2$, więc
$(\Delta a_\mu/\Delta a_e) \approx (m_\mu/m_e)^2 \approx 43000$ w TGP,
vs. "BSM universal" ≈ 1, vs. "mass-dependent BSM" ≈ $(m_\mu/m_e)^2$.
**To test topologii sprzężenia TGP do leptonu**.

### G5 — Tau: przewidywanie

Skoro TGP modyfikuje $a_\mu$, powinno też modyfikować $a_\tau$. Obecny limit:
$-0.052 < a_\tau < 0.013$ (DELPHI 2004). Skala TGP ~$(m_\tau/m_\mu)^2 \approx 280$
razy większa niż $\mu$ — potencjalnie testowalne w Belle-II lub FCC-ee.

## 4. Plan numeryczny

- **ps01_schwinger_tgp_1loop.py** — 1-pętla QED w TGP-metryce; $\Delta a_\mu$
  od samej geometrii.
- **ps02_hvp_tgp_correction.py** — integracja kernelu HVP z $\Phi$-modyfikacją.
- **ps03_hlbl_tgp.py** — HLbL z nieliniowością $\Phi^4$.
- **ps04_lepton_universality_tgp.py** — skan $(a_e, a_\mu, a_\tau)$ jako
  funkcja $\Lambda_E^\text{lepton}$.
- **ps05_summary_muon_g2.py** — końcowe porównanie z FNAL Run-6 + J-PARC.

## 5. Literatura startowa

- Aoyama et al., *The anomalous magnetic moment of the muon in the Standard Model*,
  Phys. Rep. 887, 1 (2020) — Muon g-2 Theory Initiative White Paper
- FNAL Muon g-2: Aguillard et al., PRL 131, 161802 (2023) — Run-1+2+3 combined
- BMW Collab., Nature 593, 51 (2021) — lattice HVP reducing tension
- Hanneke-Fogwell-Gabrielse, PRL 100, 120801 (2008) — electron g-2
- DELPHI: Abdallah et al., Eur. Phys. J. C 35, 159 (2004) — tau g-2 bound

## 6. Relacje z innymi sektorami TGP

- **particle_sector_closure + brannen_sqrt2 + cabibbo_correction**: masy
  i mieszanie leptonów z substratu. $(g-2)$ jest naturalnym dodatkowym
  testem.
- **metric_ansatz**: efektywna metryka wokół punktowego mionu.
- **uv_completion**: zachowanie TGP na skali ~$m_\mu c_0$ — ważne dla
  regularności pętli.
- **qm_foundations**: TGP poprawka do propagatora kwantowego.

## 7. Falsyfikowalność

**To jest jeden z czystszych testów TGP**:
- $(g-2)_\mu$ mierzony na ~$10^{-11}$ — mały margines dla BSM
- Do 2028 Fermilab + J-PARC powinny obniżyć błąd eksperymentalny o ~2×
- Lattice HVP ma ~3 grupy dający spójne wyniki → rozbieżność zostanie
  ustalona w ~1 roku
- Przy rozdzielczości $10^{-11}$, predykcja TGP albo trafia, albo nie

**Warunek na "sukces"**: TGP przewiduje $\Delta a_\mu, \Delta a_e, \Delta a_\tau$
z **jednego** parametru ($\Lambda_E^\text{lepton}$), wszystkie trzy muszą zgadzać się
z eksperymentem w $2\sigma$.

## 8. Otwarte pytania

1. Czym jest "TGP-Compton skala" dla leptonu — $\hbar c/m_\ell c_0^2$,
   czy $\hbar c/\Lambda_E^\text{univ}$?
2. Czy substrat wpływa na pętle QED ALE NIE na mass renormalisation leptonu
   (inaczej dostalibyśmy sprzeczność z precyzyjnymi pomiarami $m_e$)?
3. Czy anomalia e+e- vs lattice HVP jest rzeczywiście "environmental" w TGP,
   czy ortogonalna do sektora TGP?

## 9. Priorytet

Ten sektor ma **najbliższy eksperymentalny horyzont** (Fermilab Run-6 finalne
dane 2027, J-PARC E34 start 2028). Jeśli TGP ma trafny $\Lambda_E^\text{lepton}$,
może być jednym z **kilku testów TGP z precyzją $10^{-10}$ w ciągu 2 lat**.

## 10. Link do rdzenia TGP

Core paper [[papers/core/tgp_core.pdf]] § *Emergent metric* + *Applications*.
Sugerowany companion paper docelowy: `papers/leptons/tgp_leptons.tex`
(po uzyskaniu konkretnych predykcji dla $a_e, a_\mu, a_\tau$).
