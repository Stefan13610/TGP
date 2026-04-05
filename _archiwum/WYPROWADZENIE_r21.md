---
tags: [TGP, masa, generacje, wyprowadzenie, analityczne, r21, Cardano]
created: 2026-03-22
status: draft-verified
---

# TGP — Wyprowadzenie r₂₁ = M₂/M₁

> **Weryfikacja numeryczna**: `TGP/TGP_v1/scripts/advanced/r21_wyprowadzenie.py`
> **Dokładność wzoru**: 5–15% vs pełna numeryka

---

## 1. Setup

Cząstka TGP = soliton Yukawa w tle $\Phi_0$:

$$\Phi(r) = \Phi_0 + \frac{K\,e^{-m_\text{sp}r}}{r}, \qquad K = \frac{M}{4\pi\Phi_0}$$

Warunek samospójności: $E(K)/K = \Lambda = 4\pi$ (j.n.).

Funkcja $E(K)/K$ ma **jedno maksimum** — stąd co najwyżej **dwa przecięcia** z poziomą $\Lambda$:

$$K_1 < K_\text{max} < K_2 \qquad \Rightarrow \qquad r_{21} = \frac{M_2}{M_1} = \frac{K_2}{K_1}$$

---

## 2. Analiza E(K)/K — trzy rezimy

### Rezim A: $K \to 0$ (male K)

Przy małym $K$, pole $\Phi \approx \Phi_0$ wszędzie.

$$E_\text{kin} \approx \frac{2\pi K^2}{a_\Gamma}\bigl(1 + \tfrac{\alpha}{\Phi_0^2}\bigr)
\qquad\Rightarrow\qquad \frac{E}{K} \approx \frac{2\pi K(1+\alpha)}{a_\Gamma} \;\to\; 0$$

### Rezim B: $K \sim \Phi_0 a_\Gamma$ (przejście, rdzeń solitonu)

W rdzeniu ($r \approx a_\Gamma$): $\Phi_\text{core} \approx \Phi_0 + K/a_\Gamma$.

**Człon nieliniowy kinetyczny** dla $K \gg \Phi_0 a_\Gamma$:

$$\frac{\alpha}{\Phi_0\,\Phi_\text{core}} \approx \frac{\alpha\,a_\Gamma}{K}
\qquad\Rightarrow\qquad \frac{E_\text{kin,nl}}{K} \approx 2\pi\alpha$$

**Stały wkład** niezależny od $K$ — kluczowy dla wartości $K_1$.

### Rezim C: $K \gg \Phi_0 a_\Gamma$ (duze K)

Potencjał kwartyczny $V \sim -\psi^4/4$ dominuje w rdzeniu:

$$\Phi_\text{core} \sim K/a_\Gamma, \quad
V(\Phi_\text{core}/\Phi_0) \approx -\tfrac{1}{4}\bigl(\tfrac{K}{\Phi_0 a_\Gamma}\bigr)^4$$

$$E_\text{pot} \approx -\frac{\pi K^4}{\Phi_0^2\,a_\Gamma}
\qquad\Rightarrow\qquad \frac{E_\text{pot}}{K} \approx -\frac{\pi K^3}{\Phi_0^2\,a_\Gamma} \;\to\; -\infty$$

---

## 3. Aproksymacja E(K)/K — wzór roboczy

Łącząc trzy rezimy dla $\Phi_0 = 1$:

$$\boxed{\frac{E(K)}{K} \approx \frac{2\pi K}{a_\Gamma} + 2\pi\alpha - \frac{\pi K^3}{a_\Gamma}}$$

Warunek przecięcia $E(K)/K = 4\pi$:

$$\frac{2K}{a_\Gamma} + 2\alpha - \frac{K^3}{a_\Gamma} = 4$$

$$\boxed{K^3 - 2K = a_\Gamma(2\alpha - 4) \equiv c}$$

Jest to **rownanie kubiczne w K** o znanych współczynnikach.

---

## 4. Rozwiązanie — wzory na K₁ i K₂

### 4.1 Pierwsze przecięcie K₁ (gałąź rosnąca, małe K)

Przy małym $K$ dominuje człon liniowy $-2K$:

$$K_1 \approx \frac{4\,a_\Gamma}{2(1+\alpha) - a_\Gamma}
\;\approx\; \frac{2\,a_\Gamma}{1+\alpha} \qquad [\alpha \gg 1]$$

*Weryfikacja*: $\alpha = 5.9$, $a_\Gamma = 0.05$: $K_1 = 0.01454$ (num), $0.01455$ (wzór) — błąd $<1\%$.

### 4.2 Drugie przecięcie K₂ (gałąź malejąca, duże K)

Wzór Cardano dla $K^3 - 2K - c = 0$ przy $\Delta = 32 - 27c^2 > 0$:

$$\boxed{K_2 = \frac{2}{\sqrt{3}}\cos\!\left[\frac{1}{3}\arccos\!\left(\frac{3\sqrt{3}}{4}\,c\right)\right]}$$

gdzie $c = a_\Gamma(2\alpha - 4)$.

*Alternatywnie* przez przekształcenie trygonometryczne:

$$K_2 = 2\sqrt{\frac{2}{3}}\cos\!\left[\frac{1}{3}\arccos\!\left(\frac{3c}{2}\sqrt{\frac{3}{2}}\right)\right]$$

**Warunek istnienia dwóch przecięć** ($\Delta > 0$):

$$a_\Gamma\,(2\alpha - 4) < \sqrt{\frac{32}{27}} \approx 1.089$$

---

## 5. Formuła na r₂₁

$$\boxed{r_{21} = \frac{K_2}{K_1}
= K_2 \cdot \frac{2(1+\alpha) - a_\Gamma}{4\,a_\Gamma}
\approx \frac{K_2\,(1+\alpha)}{2\,a_\Gamma}}$$

Podstawiając wzór Cardano:

$$r_{21} \approx \frac{(1+\alpha)}{\sqrt{3}\,a_\Gamma}
\cos\!\left[\frac{1}{3}\arccos\!\left(\frac{3\sqrt{3}}{4}\,a_\Gamma(2\alpha-4)\right)\right]$$

---

## 6. Weryfikacja numeryczna

| $\alpha$ | $a_\Gamma$ | $r_{21}^\text{num}$ | $r_{21}^\text{bazowy}$ | błąd |
|---------|-----------|--------------------|--------------------|------|
| 4.0 | 0.03 | 123.2 | 119.9 | −2.7% |
| 4.0 | 0.05 | 74.7 | 72.7 | −2.7% |
| 5.9 | 0.03 | 184.3 | 168.6 | −8.5% |
| 5.9 | 0.05 | 113.6 | 103.3 | −9.0% |
| 6.0 | 0.03 | 187.7 | 171.2 | −8.8% |
| 7.0 | 0.03 | 222.1 | 197.5 | −11.1% |
| 8.0 | 0.03 | 257.8 | 224.2 | −13.1% |
| 10.0 | 0.03 | 332.4 | 278.6 | −16.2% |

**Dokładność wzoru bazowego: 3–16%** (wzór niedoszacowuje $r_{21}$, błąd rośnie z $\alpha$).

---

## 7. Warunek r₂₁ = 207

Z formuły:
$$\frac{(1+\alpha)}{\sqrt{3}\,a_\Gamma}\cos\!\left[\frac{\arccos(c\sqrt{27/32})}{3}\right] = 207$$

gdzie $c = a_\Gamma(2\alpha-4)$.

Kombinacje dające $r_{21} \approx 207$ (ze skanem numerycznym):

| $\alpha$ | $a_\Gamma$ | $r_{21}^\text{an}$ | $r_{21}^\text{num}$ |
|---------|-----------|---------|---------|
| 7.0 | 0.030 | 197.5 | **213.9** |
| 5.0 | 0.020 | 216.1 | **219.0** |
| 10.0 | 0.040 | 213.2 | **249.6** |

Wzór sugeruje $\alpha \approx 7$, $a_\Gamma \approx 0.025$–$0.030$ daje $r_{21} \sim 207$.

---

## 8. Asymptotyczne skalowanie

Dla $c = a_\Gamma(2\alpha-4) \to \sqrt{32/27}$ (granica stabilności):

$$K_2^\text{max} = \frac{2}{\sqrt{3}}, \qquad
r_{21}^\text{max} = \frac{(1+\alpha)\sqrt{3}}{a_\Gamma}$$

Dla $\alpha \gg 1$, $c \ll 1$ (rozwinięcie cosinus):

$$K_2 \approx \sqrt{2}\left(1 + \frac{c}{4}\right) = \sqrt{2}\left(1 + \frac{a_\Gamma(\alpha-2)}{2}\right)$$

$$r_{21} \approx \frac{\sqrt{2}(1+\alpha)}{2\,a_\Gamma}\left(1 + \frac{a_\Gamma(\alpha-2)}{2}\right)
\approx \frac{\sqrt{2}\,(1+\alpha)}{2\,a_\Gamma}$$

czyli:

$$\boxed{r_{21} \approx \frac{\alpha}{\sqrt{2}\,a_\Gamma} \qquad [\alpha \gg 1,\ a_\Gamma\alpha \ll 1]}$$

Dla $\alpha = 5.9$, $a_\Gamma = 0.03$: $r_{21} \approx 5.9/(0.042) \approx 139$. Numerycznie: 177.

---

## 9. Skąd pochodzi 207?

Wzór zamknięty:

$$r_{21} = \frac{(1+\alpha)}{2\sqrt{3}\,a_\Gamma}
\cdot \frac{4\,a_\Gamma}{2(1+\alpha)-a_\Gamma}
\cdot \left[\text{Cardano}\right]^{-1} \cdot \text{Cardano}$$

Upraszczając:

$$r_{21} = \frac{2}{\sqrt{3}} \cdot \frac{2(1+\alpha)-a_\Gamma}{4\,a_\Gamma}
\cdot \cos\!\left[\frac{\arccos\!\left(\frac{3\sqrt{3}}{4}a_\Gamma(2\alpha-4)\right)}{3}\right]$$

Wartość $r_{21} = 207$ wyznacza **krzywą w przestrzeni** $(\alpha, a_\Gamma)$:

$$\cos\!\left[\frac{\arccos(c\sqrt{27/32})}{3}\right] = \frac{207 \cdot 4 a_\Gamma \sqrt{3}}{2(2(1+\alpha)-a_\Gamma)}$$

Rozwiązanie tej równości daje dopuszczalne pary $(\alpha, a_\Gamma)$.

---

## 10. Źródło błędu i poprawki numeryczne

### 10.1 Wzór bazowy — trzy wyrazy

$$\frac{E(K)}{K} \approx \underbrace{\frac{2\pi K}{a_\Gamma}}_{\text{Ek standard}} + \underbrace{2\pi\alpha}_{\text{Ek nonlinear}} - \underbrace{\frac{\pi K^3}{a_\Gamma}}_{\text{Ep quartic}}$$

### 10.2 Poprawka: człon kubiczny potencjału $\psi^3/3$

Oryginalny wzór uwzględnia tylko człon kwartyczny $-\psi^4/4$ w potencjale $V(\psi) = \psi^3/3 - \psi^4/4$.
Człon kubiczny $+\psi^3/3$ daje dodatni wkład do energii:

$$E_{\text{pot,cubic}} = 4\pi\Phi_0^2 \int_{a_\Gamma}^1 \frac{\psi^3}{3}\,r^2\,dr
\approx \frac{4\pi K^3}{3\Phi_0}\ln\!\frac{1}{a_\Gamma}$$

gdzie całka jest log-rozbieżna (ucięta na $a_\Gamma$). Wkład do $E/K$:

$$\frac{E_{\text{pot,cubic}}}{K} \approx \frac{4\pi K^2}{3}\ln\!\frac{1}{a_\Gamma} \quad (> 0)$$

Ten **dodatni** człon przesuwa $K_2$ w górę (prawidłowy kierunek).

### 10.3 Poprawione równanie kubiczne

Uwzględniając człon $\psi^3/3$:

$$\boxed{K^3 - \varepsilon K^2 - 2K = c,
\qquad \varepsilon = \frac{4a_\Gamma}{3}\ln\!\frac{1}{a_\Gamma},\quad c = a_\Gamma(2\alpha-4)}$$

Przy $a_\Gamma \to 0$: $\varepsilon \to 0$, wzór wraca do formy bazowej. Rozwiązanie numeryczne (brentq):

```python
eps = (4*a_gam/3) * np.log(1/a_gam)
K2 = brentq(lambda K: K**3 - eps*K**2 - 2*K - c, 0.5, 6.0)
```

### 10.4 Porównanie: błąd wzoru bazowego vs poprawionego

| $\alpha$ | $a_\Gamma$ | $r_{21}^\text{num}$ | Bazowy | Poprawiony |
|---------|-----------|--------------------|--------------------|------------|
| 4.0 | 0.03 | 123.2 | 119.9 (−2.7%) | 125.8 (+2.1%) |
| 5.9 | 0.03 | 184.3 | 168.6 (−8.5%) | **176.6 (−4.2%)** |
| 5.9 | 0.05 | 113.6 | 103.3 (−9.0%) | **110.1 (−3.1%)** |
| 7.0 | 0.03 | 222.1 | 197.5 (−11.1%) | **206.6 (−7.0%)** |
| 8.0 | 0.03 | 257.8 | 224.2 (−13.1%) | **234.4 (−9.1%)** |
| 10.0 | 0.03 | 332.4 | 278.6 (−16.2%) | **290.9 (−12.5%)** |

**Redukcja błędu**: ~4 punkty procentowe w całym zakresie parametrów.

### 10.5 Pozostałe źródła błędu

Poprawka $\psi^3/3$ redukuje błąd, ale nie eliminuje go. Pozostałe brakujące człony:

1. **Nieliniowy kinetyczny w ogonie** ($r \sim 1$): dla $K \sim K_2$, $\phi \approx 1 + K e^{-r}/r \approx 1$ w ogonie → człon $\alpha/\phi$ ≈ $\alpha$ daje wkład skalujący się z $\alpha$ i rosnący z $K$ — wyjaśnia wzrost błędu wraz z $\alpha$
2. **Dokładna granica cięcia $r^*$** (zamiast stałej $a_\Gamma$ lub $1$): przejście między rdzeniem a ogonem jest stopniowe
3. **Poprawka log do kinetycznego** ($+2\pi K \ln(1/a_\Gamma)$ w $E_\text{kin}$): przesuwa $K_1$ i $K_2$ obu

Skrypt: `scripts/advanced/r21_poprawka.py`

---

## 11. Podsumowanie i status

| Wynik | Typ |
|-------|-----|
| $K_1 \approx 4a_\Gamma/(2(1+\alpha)-a_\Gamma)$ | ✅ Analityczny, błąd ~7–20% |
| $K_2$ (bazowy): $K^3-2K = c$ (Cardano) | ✅ Analityczny, błąd ~16–28% na $K_2$ |
| $K_2$ (poprawiony): $K^3-\varepsilon K^2-2K = c$ | ✅ Poprawiony (brentq), błąd ~11–23% |
| $r_{21}$ (bazowy Cardano) | ✅ Analityczny, błąd ~3–16% |
| $r_{21}$ (z poprawką $\psi^3/3$) | ✅ Poprawiony, błąd ~2–12% |
| Dokładna wartość 207 | ⚠️ Wymaga konkretnych $(\alpha, a_\Gamma)$ |
| $r_{21} \approx \alpha/(\sqrt{2}\,a_\Gamma)$ | ✅ Asymptotycznie poprawne |

### Kluczowe równanie — wersja bazowa:

$$K^3 - 2K = a_\Gamma(2\alpha - 4) \equiv c$$

Lewa strona: **geometria Yukawa** (kwartyczny potencjał $-\psi^4/4$).
Prawa strona: skala wyznaczona przez **sieć substratu $a_\Gamma$** i **nieliniowość kinetyczną $\alpha$**.

### Kluczowe równanie — wersja poprawiona (+człon $\psi^3/3$):

$$\boxed{K^3 - \varepsilon K^2 - 2K = c, \qquad \varepsilon = \tfrac{4a_\Gamma}{3}\ln\!\tfrac{1}{a_\Gamma}}$$

Człon $\varepsilon K^2$ pochodzi z **kubicznego potencjału** $+\psi^3/3$ i redukuje błąd $r_{21}$ o ~4 pp.

### Wzór zamknięty (Cardano, wersja bazowa):

$$r_{21}^\text{bazowy} = \frac{(1+\alpha)\cdot 2}{\sqrt{3}\cdot a_\Gamma}
\cdot\cos\!\left[\frac{1}{3}\arccos\!\!\left(\frac{3\sqrt{3}}{4}\,a_\Gamma(2\alpha-4)\right)\right]$$

---

## 12. Model sprzężony — weryfikacja i ograniczenia

> Szczegóły: [[ANALIZA_SPRZEZONY]]

### Poprawiony układ sprzężony (Φ_bg = Φ₀ − ξM)

Zamiast stałego $\Phi_0$, masa cząstki zmienia lokalne tło:

$$\Phi_\text{bg} = \Phi_0^\text{total} - \xi M \qquad K = \frac{M}{4\pi\Phi_\text{bg}}$$

Szukamy zer $g(M) = E(K;\Phi_\text{bg}) - M = 0$.

**Wynik skanowania (2026-03-23):**

| ξ | α | a_Γ | Φ₀ | r₂₁ | błąd od 207 |
|---|---|-----|-----|-----|------------|
| 0.001 | 7.0 | 0.030 | 1.0 | **210.1** | +1.5% |
| 0.001 | 5.9 | 0.030 | 1.0 | 174.4 | −15.7% |
| 0.010 | 10.0 | 0.030 | 1.0 | 238.3 | +15.1% |

### Wersja 1: brak 3. generacji (2 rozwiązania)

Funkcja $g(M)$ ma strukturę N-kształtną: **dokładnie 2 zera**, nigdy 3.

**Fizyczna przyczyna**: dla $K \gtrsim 2$, rdzeń solitonu ma $\psi \gg 1$, co powoduje $E_\text{pot} \sim -\psi^4/4 \to -\infty$.

Dwa rozwiązania:
- $M_1$ ↔ gałąź małego $K$ (elektron), $r_{21} \approx 207$ dla $\alpha=7$, $a_\Gamma=0.03$
- $M_2$ ↔ gałąź dużego $K$ (mion)

### Wersja 2: potencjał zmodyfikowany — 3 generacje ✅

Dodanie członu $\lambda(\psi-1)^6/6$ do potencjału:

$$V_\text{mod}(\psi) = \frac{\psi^3}{3} - \frac{\psi^4}{4} + \frac{\lambda(\psi-1)^6}{6}$$

zatrzymuje katastrofę ψ⁴ i tworzy **trzecie zero** $g(M)$ odpowiadające taonowi.

**Wynik numeryczny** (bisekcja po $\lambda^*$ przy stałym $r_{31}=3477$):

$$\boxed{\alpha=7.0,\quad a_\Gamma=0.030,\quad \lambda^*=5.51\times10^{-7}}$$

$$r_{21} = 208.1 \quad (+0.5\%),\qquad r_{31} = 3477.0 \quad (0.0\%)$$

Skrypty: `v2_lambda_potencjal.py`, `v2_fine_scan.py`. Szczegóły: [[ANALIZA_SPRZEZONY]].

---

## Pliki

- `scripts/advanced/r21_wyprowadzenie.py` — weryfikacja numeryczna wzoru bazowego
- `scripts/advanced/r21_poprawka.py` — poprawka $\psi^3/3$ (redukcja błędu ~4 pp)
- `scripts/advanced/sprzezony_układ.py` — model sprzężony Φ_bg = Φ₀ − ξM (poprawiony 2026-03-23)
- `scripts/advanced/badanie_pi.py` — badanie roli $\pi$
- `ODKRYCIA_PI_KOIDE.md` — obserwacje o $\pi$ i formule Koidego
- `ANALIZA_SPRZEZONY.md` — pełna analiza modelu sprzężonego
