# gs9 — Program badawczy: efektywna redukcja wymiarowa w TGP

## Kontekst

Zbadano 5 mechanizmów mikro (gs7a-gs7f) i sprzężenie dysformalne (gs8).
Wszystkie perturbacyjne mechanizmy OBALONO. Ale analiza gs8 doprowadziła
do nowego kierunku: **efektywne przejście wymiarowe 3D → 2D**.

### Kluczowe odkrycia (gs8)

1. Model z efektywną redukcją do 2D naturalnie daje:
   - F ~ 1/r² blisko (zachowanie 3D) → Newton
   - F ~ 1/r daleko (efektywne 2D) → płaska RC
   - v⁴ = GM·a₀ automatycznie (Tully-Fisher)
2. Parametr przejścia κ = a₀/(cH₀) = 1/(2π) — naturalna wartość
3. Model DGP-like: g_obs = g_N + √(g_N·a₀) — zgodny z MOND w limitach
4. Różnica od MOND: ~27% przy g_bar ≈ a₀ — testowalny z SPARC

### Co obserwujemy — co to znaczy

**TGP wykazuje sygnały efektywnej redukcji do 2D.**

To NIE musi oznaczać, że przestrzeń jest dosłownie 2D (holografia, braneworld).
Możliwe interpretacje — które z nich wyłoni się z równań:

**Interpretacja A: Wymiar efektywny dla dynamiki informacji**
2D jest wymiarem w którym działa właściwa miara, geometria lub przepływ
stopni swobody. Obserwowalna przestrzeń pozostaje 3D, ale korelacje,
propagacja informacji grawitacyjnej lub entropia mają charakter 2D
na dużych skalach.

**Interpretacja B: Sektor matematyczny**
2D jest specjalnym sektorem w którym rdzeń TGP prostuje się i pokazuje
swoją właściwą strukturę. Podobnie jak równanie solitonowe jest
najczyściej widoczne w 1D, a fizycznie realizowane w 3D — tak dynamika
grawitacji na dużych skalach „preferuje" 2D opis.

**Interpretacja C: Fizyczna membrana/brana**
Substrat TGP jest dosłownie 2D obiektem w wyższym wymiarze (holografia,
braneworld). Najbardziej radykalna opcja — nie zakładamy jej z góry.

**Podejście**: Niech równania pokażą, który obraz jest właściwy.
Nie fiksujemy się na żadnej interpretacji — badamy matematykę
przejścia 3D→2D i patrzymy co wyłania się naturalnie.

## Plan badawczy

### Etap 1: Struktura matematyczna przejścia (gs9a)

**Cel**: Zbadać MATEMATYCZNIE jak równanie TGP zachowuje się
w różnych wymiarach i czy przejście 3D→2D pojawia się naturalnie.

**Co policzyć**:
1. Równanie solitonowe w d wymiarach:
   g'' + g'²/g + (d-1)g'/r + g = 1
   Jak zmienia się rozwiązanie z d? Ogon, masa, energia?
2. Funkcja Greena liniowego TGP (δ'' + (d-1)δ'/r + δ = source) w 2D i 3D:
   - 3D: sin(r)/r (oscylujący, zanika jak 1/r)
   - 2D: J₀(r) (oscylujący, zanika jak 1/√r — wolniej!)
   Co się dzieje gdy d zmienia się płynnie od 3 do 2?
3. Efektywny wymiar: jeśli zdefiniujemy d_eff(r) z zachowania potencjału,
   jak d_eff zależy od odległości? Czy jest sygnał d_eff → 2 przy r → ∞?
4. Lagrangian TGP: L = |∇g|²/g + (g-1)²
   W jakim wymiarze d ma najczystszą strukturę (symetrie, całkowalność)?

**Kryteria sukcesu**:
- [ ] Zachowanie rozwiązania jako funkcja d zbadane numerycznie
- [ ] Zidentyfikowany mechanizm preferujący d=2 na dużych skalach
- [ ] Lub: wykazane, że żaden taki mechanizm nie istnieje → odrzucenie

### Etap 2: Efektywny propagator z przejściem wymiarowym (gs9b)

**Cel**: Skonstruować propagator grawitacji który zachowuje się
jak 3D blisko i 2D daleko, i porównać z obserwacjami.

**Co policzyć**:
1. Propagator z „efektywnym wymiarem" d(r):
   G(r) rozwiązanie ∇²_d G + G = δ(r) z d = d(r)
2. Najprostszy model: d(r) = 2 + 1/(1 + r/r_c) (przejście 3→2)
3. Wynikowa siła F(r) i krzywa rotacji v²(r)
4. RAR: g_obs(g_bar) — dokładna forma interpolacji
5. Porównanie z MOND, DGP, i danymi

**Uwaga**: Nie zakładamy DLACZEGO d→2. Badamy CZY model z d→2
pasuje do obserwacji lepiej niż MOND czy CDM.

**Kryteria sukcesu**:
- [ ] Gładkie przejście dające płaską RC
- [ ] BTFR z wykładnikiem 4 odtworzony
- [ ] RAR porównywalna z MOND lub lepsza
- [ ] r_c powiązany z a₀: r_c = √(GM/a₀)

### Etap 3: Porównanie z danymi SPARC (gs9c)

**Cel**: Ilościowe porównanie z obserwacjami.

**Co policzyć**:
1. RAR: g_obs vs g_bar — χ² vs MOND
2. Krzywe rotacji wybranych galaktyk
3. Freeman limit z modelu
4. Scatter w RAR
5. Klastry: czy model daje więcej „phantom DM" niż MOND?

**Kryteria sukcesu**:
- [ ] χ² ≤ MOND na danych SPARC
- [ ] Freeman limit w zakresie 100-200 M☉/pc²
- [ ] Predykcje dla klastrów

### Etap 4: Dlaczego 2D? — szukanie mechanizmu (gs9d)

**Cel**: Dopiero PO potwierdzeniu że przejście 3D→2D działa fenomenologicznie,
szukać DLACZEGO równania TGP preferują 2D.

**Możliwe kierunki** (nie przesądzamy):
- Dynamika informacji/korelacji na substracie
- Holograficzna interpretacja (jeśli wyłoni się z matematyki)
- Matematyczny sektor czysty w 2D (jak całkowalność w 1+1D)
- Efektywna teoria pola: RG flow ku d=2 w podczerwieni
- Napięcie powierzchniowe substratu → Λ

**Kryterium**: Mechanizm musi wynikać z równań, nie z założeń.

### Etap 5: Unikalne predykcje (gs9e)

**Predykcje do identyfikacji**:
1. Interpolacja RAR: konkretna forma f(g_bar/a₀) ≠ MOND
2. a₀(z): ewoluuje z H(z) czy stałe?
3. Efekt Zewnętrznego Pola (EFE)
4. Klastry i Bullet Cluster
5. Ultradiffuse galaxies
6. Testy Układu Słonecznego

## Kryteria odrzucenia

Model jest **odrzucony** jeśli:
- RAR gorsza niż MOND o >3σ na SPARC
- Nie daje BTFR z wykładnikiem 4
- Sprzeczność z testami Układu Słonecznego
- Nie daje Freeman limit w zakresie 100-200 M☉/pc²
- Przejście 3D→2D wymaga fine-tuningu

## Kluczowa zasada

> Nie zakładamy holografii, braneworldów ani membran.
> Obserwujemy sygnał matematyczny (efektywne 2D na dużych skalach)
> i badamy go rzetelnie. Niech równania pokażą interpretację.

## Powiązania

- [[gs8_disformal_coupling.py]] — analiza przejścia 3D→2D
- [[gs7f_comparison.py]] — porównanie 5 mechanizmów
- [[gs1_flat_well_model.py]] — fenomenologia do reprodukcji
- [[ct6_mechanism_diagnosis.py]] — problem dwóch skal
