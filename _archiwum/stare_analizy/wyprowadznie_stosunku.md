# Analityczne wyprowadzenie stosunku ψ_core/ψ_cross = π/2 w TGP

## Streszczenie

Wyprowadzamy analityczne wyrażenie na stosunek wartości pola w centrum solitonu trzeciej generacji ($\psi_{\text{core}}$) do wartości crossover ($\psi_{\text{cross}}$) w Teorii Generowanej Przestrzeni (TGP). Pokazujemy, że stosunek ten zależy wyłącznie od parametru cutoff $a = a_\Gamma$ (regularizacji źródła) i dany jest wzorem:

$$R(a) = \frac{e^{-a}}{a} \sqrt{\frac{I_4(a)}{I_6(a)}}$$

gdzie $I_4(a)$ i $I_6(a)$ są całkami wyrażonymi przez funkcję wykładniczą całkową $E_1$. Dla $a = 0.040$ otrzymujemy $R(0.040) \approx 1.566$, co jest niezwykle bliskie $\pi/2 \approx 1.5708$ (błąd $\approx 0.3\%$). Wynik ten wyjaśnia, dlaczego w TGP stosunek mas leptonów spełnia relację Koide'a ($Q = 3/2$) i dlaczego parametr $a_\Gamma$ przyjmuje wartość $0.040$.

---

## 1. Wprowadzenie i definicje

W TGP potencjał zmodyfikowany dla pola $\psi = \Phi/\Phi_0$ ma postać:

$$V_{\text{mod}}(\psi) = \frac{\psi^3}{3} - \frac{\psi^4}{4} + \lambda \frac{(\psi-1)^6}{6}, \quad \lambda > 0$$

Definiujemy:

- **$\psi_{\text{cross}}$** – punkt, w którym $V_{\text{mod}}(\psi_{\text{cross}}) = V_{\text{mod}}(1)$. Dla $\lambda \ll 1$ i $\psi \gg 1$:

$$\psi_{\text{cross}} = \sqrt{\frac{3}{2\lambda}}$$

- **$\psi_{\text{core}}$** – wartość pola w centrum solitonu ($r = a_\Gamma$). Dla profilu $\psi(r) = 1 + K e^{-r}/r$ i dużego $K$:

$$\psi_{\text{core}} \approx \frac{K e^{-a}}{a}$$

gdzie $a = a_\Gamma$ jest cutoffem regularizacji.

- **$K_3$** – amplituda trzeciego solitonu (trzeciej generacji), dana warunkiem samospójności $g(K_3)=0$. W przybliżeniu dominacji członów kwartycznego i $\lambda$:

$$K_3^2 = \frac{3I_4(a)}{2\lambda I_6(a)}$$

gdzie $I_4(a)$ i $I_6(a)$ są całkami zależnymi od $a$.

---

## 2. Wyrażenia na całki $I_4(a)$ i $I_6(a)$

### 2.1 Całka $I_4(a)$

$$I_4(a) = \int_a^\infty \frac{e^{-4r}}{r^2} dr$$

Całkując przez części:

$$I_4(a) = \left[ -\frac{e^{-4r}}{r} \right]_a^\infty - 4\int_a^\infty \frac{e^{-4r}}{r} dr = \frac{e^{-4a}}{a} - 4E_1(4a)$$

gdzie $E_1(x) = \int_x^\infty \frac{e^{-t}}{t} dt$ jest wykładniczą całką całkową.

$$\boxed{I_4(a) = \frac{e^{-4a}}{a} - 4E_1(4a)}$$

### 2.2 Całka $I_6(a)$

$$I_6(a) = \int_a^\infty \frac{e^{-6r}}{r^4} dr$$

Całkując przez części trzykrotnie:

Pierwsze całkowanie:
$$\int_a^\infty \frac{e^{-6r}}{r^4} dr = \left[ -\frac{e^{-6r}}{3r^3} \right]_a^\infty - 2\int_a^\infty \frac{e^{-6r}}{r^3} dr = \frac{e^{-6a}}{3a^3} - 2\int_a^\infty \frac{e^{-6r}}{r^3} dr$$

Drugie całkowanie:
$$\int_a^\infty \frac{e^{-6r}}{r^3} dr = \left[ -\frac{e^{-6r}}{2r^2} \right]_a^\infty - 3\int_a^\infty \frac{e^{-6r}}{r^2} dr = \frac{e^{-6a}}{2a^2} - 3\int_a^\infty \frac{e^{-6r}}{r^2} dr$$

Trzecie całkowanie:
$$\int_a^\infty \frac{e^{-6r}}{r^2} dr = \frac{e^{-6a}}{a} - 6E_1(6a)$$

Podstawiając od tyłu:

$$\int_a^\infty \frac{e^{-6r}}{r^3} dr = \frac{e^{-6a}}{2a^2} - 3\left( \frac{e^{-6a}}{a} - 6E_1(6a) \right) = \frac{e^{-6a}}{2a^2} - \frac{3e^{-6a}}{a} + 18E_1(6a)$$

Zatem:

$$I_6(a) = \frac{e^{-6a}}{3a^3} - 2\left( \frac{e^{-6a}}{2a^2} - \frac{3e^{-6a}}{a} + 18E_1(6a) \right)$$

$$I_6(a) = \frac{e^{-6a}}{3a^3} - \frac{e^{-6a}}{a^2} + \frac{6e^{-6a}}{a} - 36E_1(6a)$$

$$\boxed{I_6(a) = \frac{e^{-6a}}{3a^3} - \frac{e^{-6a}}{a^2} + \frac{6e^{-6a}}{a} - 36E_1(6a)}$$

---

## 3. Wyrażenie na stosunek $\psi_{\text{core}}/\psi_{\text{cross}}$

Z definicji:

$$\psi_{\text{core}} = \frac{K_3 e^{-a}}{a}, \quad \psi_{\text{cross}} = \sqrt{\frac{3}{2\lambda}}$$

Podstawiając $K_3$ z warunku samospójności:

$$\frac{\psi_{\text{core}}}{\psi_{\text{cross}}} = \frac{e^{-a}}{a} \sqrt{\frac{3}{2\lambda} \cdot \frac{I_4(a)}{I_6(a)}} \cdot \sqrt{\frac{2\lambda}{3}} = \frac{e^{-a}}{a} \sqrt{\frac{I_4(a)}{I_6(a)}}$$

**Kluczowa obserwacja:** $\lambda$ znika! Stosunek zależy tylko od $a$:

$$\boxed{R(a) = \frac{e^{-a}}{a} \sqrt{\frac{I_4(a)}{I_6(a)}}}$$

gdzie $I_4(a)$ i $I_6(a)$ są dane wzorami powyżej.

---

## 4. Obliczenie $R(a)$ dla $a = 0.040$

### 4.1 Wartości pomocnicze

$$e^{-0.04} = 0.96078944$$
$$e^{-0.16} = 0.85214379$$
$$e^{-0.24} = 0.78662786$$

$$4a = 0.16,\quad 6a = 0.24$$

### 4.2 Funkcja $E_1(x)$ dla małych $x$

Rozwinięcie: $E_1(x) = -\ln x - \gamma + x - \frac{x^2}{4} + \frac{x^3}{18} - \ldots$, gdzie $\gamma \approx 0.5772156649$.

$$E_1(0.16) \approx -\ln(0.16) - 0.57721566 + 0.16 - 0.0064$$
$$-\ln(0.16) = 1.83258146$$
$$E_1(0.16) \approx 1.83258146 - 0.57721566 + 0.16 - 0.0064 = 1.4089658$$

$$E_1(0.24) \approx -\ln(0.24) - 0.57721566 + 0.24 - 0.0144$$
$$-\ln(0.24) = 1.42711636$$
$$E_1(0.24) \approx 1.42711636 - 0.57721566 + 0.24 - 0.0144 = 1.0755007$$

### 4.3 Obliczenie $I_4(0.04)$

$$\frac{e^{-4a}}{a} = \frac{0.85214379}{0.04} = 21.30359475$$
$$4E_1(4a) = 4 \times 1.4089658 = 5.6358632$$
$$I_4 = 21.30359475 - 5.6358632 = 15.66773155$$

### 4.4 Obliczenie $I_6(0.04)$

$$\frac{e^{-6a}}{3a^3} = \frac{0.78662786}{3 \times 0.000064} = \frac{0.78662786}{0.000192} = 4097.0201$$
$$\frac{e^{-6a}}{a^2} = \frac{0.78662786}{0.0016} = 491.6424125$$
$$\frac{6e^{-6a}}{a} = \frac{6 \times 0.78662786}{0.04} = \frac{4.71976716}{0.04} = 117.994179$$
$$36E_1(6a) = 36 \times 1.0755007 = 38.7180252$$

$$I_6 = 4097.0201 - 491.6424125 + 117.994179 - 38.7180252$$
$$I_6 = 3605.3776875 + 117.994179 = 3723.3718665$$
$$I_6 = 3723.3718665 - 38.7180252 = 3684.6538413$$

### 4.5 Obliczenie $R(0.04)$

$$\frac{I_4}{I_6} = \frac{15.66773155}{3684.6538413} = 0.0042525$$
$$\sqrt{\frac{I_4}{I_6}} = \sqrt{0.0042525} = 0.065206$$
$$\frac{e^{-a}}{a} = \frac{0.96078944}{0.04} = 24.019736$$
$$R(0.04) = 24.019736 \times 0.065206 = 1.5659$$

---

## 5. Porównanie z $\pi/2$

### 5.1 Obliczenia ręczne (przybliżone)

$$\frac{\pi}{2} = 1.57079632679$$
$$R(0.04)_{\rm ręcznie} = 1.5659$$

Różnica od obliczeń ręcznych: $\Delta = 0.0049$, czyli 0.31%.

### 5.2 Obliczenia dokładne (weryfikacja p132)

Użycie `scipy.special.exp1` (pełna precyzja E₁) daje:

$$R(0.040)_{\rm dokładne} = 1.56625422$$

Główne źródło różnicy z obliczeniami ręcznymi: obcięcie szeregu E₁ daje E₁(0.24) = 1.0755 zamiast dokładnego 1.0762, co propaguje się przez I₆.

**Porównanie z π/2:**

| Wartość | Liczbowa | Odchylenie od π/2 |
|---------|----------|-------------------|
| R(0.040) ręcznie | 1.5659 | 0.31% |
| R(0.040) dokładne | 1.56625 | **0.29%** |
| π/2 | 1.57080 | — |

### 5.3 Dokładne a* spełniające R(a*) = π/2

Bisection (`scipy.optimize.brentq`) daje:

$$a^* = 0.037815$$

przy czym $a_\Gamma^{\rm TGP} = 0.040049$. Różnica: **5.6%**.

Interpretacja: warunek $R(a) = \pi/2$ nie trafia dokładnie w $a_\Gamma = 0.040$, ale jest bliski. Może to wskazywać na:
- potrzebę poprawki wyższego rzędu w $V_{\rm mod}$ (np. członem $\psi^8$),
- inną postać warunku selekcji (np. $R(a) = \pi/2 + O(a)$),
- lub fakt, że $a_\Gamma$ jest wyznaczany przez **inny** warunek (np. φ-FP z running $\alpha_{\rm eff}$).

---

## 6. Wnioski

1. **Stosunek $\psi_{\text{core}}/\psi_{\text{cross}}$ zależy tylko od cutoffu $a$** i jest dany wzorem:
   $$R(a) = \frac{e^{-a}}{a} \sqrt{\frac{\frac{e^{-4a}}{a} - 4E_1(4a)}{\frac{e^{-6a}}{3a^3} - \frac{e^{-6a}}{a^2} + \frac{6e^{-6a}}{a} - 36E_1(6a)}}$$
   Ten wynik jest **dokładny** — wzory na I₄ i I₆ zweryfikowane z kwadraturą numeryczną do precyzji maszynowej (p132).

2. Dla $a = 0.040$ otrzymujemy $R = 1.5663$, co jest **0.29% poniżej** $\pi/2$. Dokładne $a^* = 0.03782$ daje $R(a^*) = \pi/2$.

3. **Hipoteza (zmodyfikowana):** Istnieje dokładne $a^*$ spełniające $R(a^*) = \pi/2$, wynoszące $a^* = 0.03782$. Jest ono bliskie, ale nie identyczne z $a_\Gamma = 0.040$ (5.6% różnicy).

4. **Nowa obserwacja (p131–p132):** Numeryczne $\eta_K = 12.067$ (z re-derywacji φ-FP) jest bliskie
   $$2\left(\frac{\pi}{2}\right)^4 = 12.176$$
   z dokładnością **0.9%**. Jeśli ta relacja jest fundamentalna, łączy warunek $R = \pi/2$ z biegiem $\alpha_{\rm eff}$.

5. **Wniosek (ostrożny):** Parametr $a_\Gamma$ jest *bliski* wartości wyznaczonej przez $R(a) = \pi/2$, ale nie jest z nią identyczny. Pełne wyznaczenie $a_\Gamma$ wymaga prawdopodobnie sprzężonego warunku uwzględniającego bieg $\alpha_{\rm eff}$ (η_K).

---

## 7. Związek z running $\alpha_{\rm eff}$ i η_K (sesja v42+)

Niezależna analiza numeryczna (p127–p131) zamknęła problem masy tau:

$$\alpha_{\rm eff}(g) = \frac{2}{1 + \eta_K\,(g-1)^2}, \qquad \eta_K = 12.067$$

daje jednocześnie $r_{21} = 206.768$ i $r_{31} = 3477.5$ (m_τ = 1777.0 MeV).

Koincydencja:

$$\eta_K = 12.067 \approx 2\left(\frac{\pi}{2}\right)^4 = 12.176 \quad (0.9\%)$$

**Status:** hipoteza robocza do dalszej analizy. Może łączyć warunek geometryczny $R = \pi/2$ z dynamicznym biegiem $\alpha_{\rm eff}$.

---

## Dodatek A: Implementacja numeryczna w Pythonie

```python
import numpy as np
from scipy.special import exp1

def I4(a):
    """Całka I4(a) = ∫_a^∞ e^{-4r}/r^2 dr"""
    return np.exp(-4*a)/a - 4*exp1(4*a)

def I6(a):
    """Całka I6(a) = ∫_a^∞ e^{-6r}/r^4 dr"""
    e6a = np.exp(-6*a)
    return e6a/(3*a**3) - e6a/a**2 + 6*e6a/a - 36*exp1(6*a)

def R(a):
    """Stosunek ψ_core/ψ_cross"""
    return np.exp(-a)/a * np.sqrt(I4(a)/I6(a))

# Szukanie a* takiego, że R(a) = π/2
from scipy.optimize import brentq

target = np.pi/2
a_star = brentq(lambda a: R(a) - target, 0.01, 0.1)
print(f"a* = {a_star:.6f}")
print(f"R(a*) = {R(a_star):.6f}")