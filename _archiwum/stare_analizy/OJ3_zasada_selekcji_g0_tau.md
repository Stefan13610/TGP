# O-J3: zasada selekcji dla \( g_0^\tau \)

## Status
**✅ ZAMKNIĘTY NUMERYCZNIE — ex113_oj3_tau_koide.py 8/8 PASS (2026-04-02)**

### Wynik końcowy

| Wielkość | Wartość | Źródło |
|---|---|---|
| g₀^e | 1.24915 | φ-FP (ex106) |
| g₀^μ | 2.02117 = φ·g₀^e | O-J2 ✓ |
| **g₀^τ** | **3.18912 ≈ φ²·g₀^e** | Koide+A_tail (ex113) |
| r₂₁ | 206.768 (PDG: 206.768, δ=0.0001%) | A_tail^4 |
| r₃₁ | 3477.44 (PDG: 3477.48, δ=0.0012%) | Koide+A_tail |

**Zasada selekcji φ-drabiny:**
- elektron: g₀^e (bazowa)
- mion: g₀^μ = φ·g₀^e (O-J2, φ-FP)
- **tau: g₀^τ ≈ φ²·g₀^e z korekcją ε_φ²≈2.5%**

Stara hipoteza g₀^τ ≈ 2·g₀^e (OJ3.md) dotyczyła innych g₀ (przed φ-FP) i jest zastąpiona przez φ²-selekcję.

---

## Status pierwotny
**Hipoteza robocza (2026-03-30) — zastąpiona przez wynik ex113**

Nie udało się jeszcze wyprowadzić w pełni analitycznie sztywnej wartości
\[
g_0^\tau = 1.7294
\]
wyłącznie z obecnych aksjomatów rdzenia TGP. Można jednak sformułować silną zasadę selekcji wiodącej:

\[
g_0^\tau \approx 2 g_0^e,
\]

a wartość dokładną \( g_0^\tau = 1.7294 \) interpretować jako niewielką korektę dynamiczną pochodzącą z pełnego ODE substratowego oraz warunku globalnego domknięcia trypletu leptonowego.

---

## Dane wejściowe

Z wcześniejszych wyników przyjmujemy:
\[
g_0^e = 0.8694,
\qquad
g_0^\mu = 1.4068,
\qquad
g_0^\tau = 1.7294.
\]

Ponadto z O-J2 mamy:
\[
g_0^\mu = \varphi\, g_0^e.
\]

Dla sektora \(\tau\) otrzymujemy:
\[
\frac{g_0^\tau}{g_0^e}
=
\frac{1.7294}{0.8694}
\approx 1.9892.
\]

Jest to wartość bardzo bliska:
\[
2.
\]

---

## Zasada selekcji wiodącej

Naturalna propozycja brzmi:
\[
\boxed{
g_0^\tau \approx 2\,g_0^e
}
\]

albo dokładniej:
\[
\boxed{
g_0^\tau = (2-\varepsilon_\tau)\, g_0^e,
\qquad
\varepsilon_\tau \ll 1.
}
\]

Dla wartości numerycznych:
\[
\varepsilon_\tau
=
2-\frac{g_0^\tau}{g_0^e}
=
2-1.9892
\approx 0.0108.
\]

Zatem odchylenie od idealnej zasady \(2 g_0^e\) wynosi około **1.1%** na poziomie ilorazu.

---

## Interpretacja fizyczna

W obecnej strukturze TGP widać dwa różne typy selekcji:

1. **Selekcja złota** dla mionu:
   \[
   g_0^\mu = \varphi g_0^e.
   \]

2. **Selekcja harmoniczna** dla tau:
   \[
   g_0^\tau \approx 2 g_0^e.
   \]

Sugeruje to, że:
- elektron odpowiada poziomowi bazowemu,
- mion jest stanem wyróżnionym przez regułę \(\varphi\),
- tau odpowiada drugiemu poziomowi harmonicznemu / radialnemu tej samej rodziny substratowej.

W tym obrazie tau nie jest kolejnym krokiem „złotej drabiny”, lecz stanem należącym do innej reguły selekcyjnej.

---

## Powiązanie z Koidem

Jeżeli:
- \(r_{21}\) jest już ustalone przez mechanizm \(\varphi\)-FP,
- a relacja Koidego jest traktowana jako warunek globalnego domknięcia trypletu \((e,\mu,\tau)\),

to \(r_{31}\) przestaje być wolnym parametrem.

Z relacji Koidego:
\[
\frac{m_e+m_\mu+m_\tau}{(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2}
=
\frac{2}{3}
\]
oraz z ustalonego \(r_{21}\) dostaje się fizyczne rozwiązanie:
\[
r_{31} \approx 3477.44.
\]

Jeśli następnie masa efektywna w TGP jest zadana przez mapę monotoniczną
\[
M(g_0)
=
\left(
\frac{A_{\rm tail}(g_0)}{A_{\rm tail}(g_0^e)}
\right)^4,
\]
to warunek
\[
M(g_0^\tau)=r_{31}
\]
wybiera jednoznacznie wartość
\[
g_0^\tau \approx 1.7294.
\]

---

## Wniosek roboczy

Najmocniejsza obecnie forma O-J3 nie brzmi:

\[
g_0^\tau = 1.7294
\quad\text{(wyprowadzone dokładnie z aksjomatów),}
\]

lecz raczej:

\[
\boxed{
\text{O-J3': }\;
g_0^\tau = (2-\varepsilon_\tau)\,g_0^e,
\qquad
\varepsilon_\tau \sim 10^{-2},
}
\]
gdzie:
- część główna \(2g_0^e\) jest zasadą selekcji harmonicznej,
- mała poprawka \(\varepsilon_\tau\) pochodzi z pełnej nieliniowej dynamiki substratowej i warunku globalnego domknięcia leptonów.

---

## Interpretacja końcowa

Najkrócej:
- **mion**: selekcja złota,
- **tau**: selekcja harmoniczna,
- **dokładna liczba \(1.7294\)**: renormalizacja poziomu \(2g_0^e\) przez pełną dynamikę substratu.

---

## Status epistemiczny

To sformułowanie należy obecnie traktować jako:

- **mocną hipotezę strukturalną** dla sektora \(\tau\),
- **nie pełne twierdzenie analityczne**.

Do pełnego wyprowadzenia brakuje jeszcze jednego z dwóch elementów:

1. albo reguły topologiczno-radialnej dającej wprost
   \[
   g_0^{(n)} \approx n g_0^e,
   \]
2. albo lokalnej analitycznej postaci funkcji \(A_{\rm tail}(g_0)\) na gałęzi \(\tau\), pozwalającej wyliczyć poprawkę \(\varepsilon_\tau\) bez odwołań czysto numerycznych.
