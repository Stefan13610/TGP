# RIBBONS_group_compare.py
# Etap: substrat -> wstazki (arena Phi-2, nieabelowa).
# Cel: policzyc U ZRODLA strukture grupowa dwoch kandydatow na target M = SO(3)/H,
#      gdzie wstazki = linie dysklinacji klasyfikowane przez pi_1(M) = H~ (binarna).
#   - Q12 = Dic3 (binarna dwuscienna, rzad 12)  = pi_1(SO(3)/D3)
#   - 2T  = SL(2,3) (binarna tetraedralna, rzad 24) = pi_1(SO(3)/T)
#
# Dla kazdej grupy raportujemy (bez fitowania, czysta teoria grup):
#   1. rzad, histogram rzedow elementow
#   2. centrum (element -1 rzedu 2 = generator spinu-1/2, FR/Berry)
#   3. klasy sprzezonosci = TYPY WSTAZEK (topologiczne ladunki linii)
#   4. podgrupa komutatorow G' i abelianizacja G/G' = ile sektorow ladunku widac ABELOWO
#   5. TEST KOLORU: czy Z3 (trialnosc) przezywa abelianizacje (obserwowalne),
#      czy tylko siedzi w strukturze klas (jak obalony no-go writhe->Z)
#
# Zadnych parametrow numerycznych, zadnego pudla, zadnego fitu. Sama struktura.

import itertools, cmath

# ---------- generyczny silnik: zamknij grupe z generatorow, policz strukture ----------

def close_group(gens, mul, ident, key):
    elems = {key(ident): ident}
    frontier = [ident]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                p = mul(g, h)
                k = key(p)
                if k not in elems:
                    elems[k] = p; new.append(p)
        frontier = new
    return list(elems.values())

def order_of(g, mul, ident, key):
    p = g; n = 1
    while key(p) != key(ident):
        p = mul(p, g); n += 1
        if n > 1000: raise RuntimeError("order blew up")
    return n

def inverse(g, elems, mul, ident, key):
    for h in elems:
        if key(mul(g, h)) == key(ident): return h
    raise RuntimeError("no inverse")

def conj_classes(elems, mul, ident, key):
    kmap = {key(e): e for e in elems}
    seen = set(); classes = []
    for e in elems:
        if key(e) in seen: continue
        cls = set()
        for g in elems:
            gi = inverse(g, elems, mul, ident, key)
            c = mul(mul(g, e), gi)
            cls.add(key(c))
        classes.append([kmap[k] for k in cls])
        seen |= cls
    return classes

def commutator_subgroup(elems, mul, ident, key):
    comms = []
    for a in elems:
        ai = inverse(a, elems, mul, ident, key)
        for b in elems:
            bi = inverse(b, elems, mul, ident, key)
            comms.append(mul(mul(a, b), mul(ai, bi)))
    # domknij podgrupe generowana przez komutatory
    return close_group(comms, mul, ident, key)

def analyze(name, elems, mul, ident, key, note=""):
    n = len(elems)
    orders = {}
    for e in elems:
        o = order_of(e, mul, ident, key)
        orders[o] = orders.get(o, 0) + 1
    # centrum
    center = []
    for z in elems:
        if all(key(mul(z, g)) == key(mul(g, z)) for g in elems):
            center.append(z)
    # klasy sprzezonosci
    classes = conj_classes(elems, mul, ident, key)
    class_info = sorted([(len(c), order_of(c[0], mul, ident, key)) for c in classes])
    # abelianizacja
    Gp = commutator_subgroup(elems, mul, ident, key)
    ab_order = n // len(Gp)
    # rzedy elementow w abelianizacji (przez rzedy warstw) -> struktura
    # warstwy: klucz warstwy = zbior kluczy g*Gp
    Gp_keys = set(key(x) for x in Gp)
    def coset_key(g):
        return frozenset(key(mul(g, x)) for x in Gp)
    cosets = {}
    for g in elems:
        ck = coset_key(g)
        cosets.setdefault(ck, g)
    # rzad warstwy w grupie ilorazowej
    coset_orders = {}
    reps = list(cosets.values())
    for r in reps:
        # rzad r modulo Gp
        p = r; m = 1
        while coset_key(p) != coset_key(ident):
            p = mul(p, r); m += 1
        coset_orders[m] = coset_orders.get(m, 0) + 1

    print("="*70)
    print(f"{name}   |G| = {n}   {note}")
    print("-"*70)
    print(f"  histogram rzedow elementow: {dict(sorted(orders.items()))}")
    print(f"  centrum |Z| = {len(center)}  (rzedy: {sorted(order_of(z,mul,ident,key) for z in center)})")
    print(f"  liczba klas sprzezonosci = {len(classes)}  (= liczba typow wstazek)")
    print(f"  klasy (rozmiar, rzad_reprezentanta): {class_info}")
    has_order3_class = any(o == 3 for _, o in class_info)
    print(f"  czy istnieje klasa rzedu 3 (kandydat-kolor w strukturze klas)? {has_order3_class}")
    print(f"  |G'| (podgrupa komutatorow) = {len(Gp)}")
    print(f"  |G/G'| (abelianizacja) = {ab_order}   struktura(rzedy warstw): {dict(sorted(coset_orders.items()))}")
    ab_has_3 = (ab_order % 3 == 0)
    print(f"  >>> TEST KOLORU: czy Z3 PRZEZYWA abelianizacje (obserwowalne abelowo)? {ab_has_3}")
    if ab_has_3:
        print(f"      => trialnosc widoczna na poziomie ABELOWYM ({ab_order} sektorow, 3 | {ab_order}).")
    else:
        print(f"      => ab.= {ab_order}, brak Z3: kolor GINIE pod abelianizacja (jak no-go writhe->Z).")
    print()
    return {"name": name, "n": n, "classes": len(classes), "abelianization": ab_order,
            "color_survives_abelianization": ab_has_3, "order3_class": has_order3_class}

# ---------- Grupa 1: Q12 = Dic3, faithful 2x2 complex rep ----------
# a = diag(z, z^-1), z = exp(i*pi/3) (6-ty pierwiastek); b = [[0,1],[-1,0]]
# a^6=I, b^2=-I=a^3, b a b^-1 = a^-1
z = cmath.exp(1j*cmath.pi/3)
def matmul(A, B):
    return ((A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]),
            (A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]))
def matkey(A):  # zaokraglony klucz
    return tuple(round(x.real,6)+0.0 for row in A for x in row) + \
           tuple(round(x.imag,6)+0.0 for row in A for x in row)
I2 = ((1+0j,0j),(0j,1+0j))
a = ((z,0j),(0j,1/z))
b = ((0j,1+0j),(-1+0j,0j))
Q12 = close_group([a,b], matmul, I2, matkey)
r1 = analyze("Q12 = Dic3  (target M = SO(3)/D3)", Q12, matmul, I2, matkey,
             note="minimalny nieabelowy z Z3")

# ---------- Grupa 2: 2T = SL(2,3), 2x2 nad GF(3), det=1 ----------
def m3mul(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) % 3 for j in range(2)) for i in range(2))
def m3key(A): return (A[0][0],A[0][1],A[1][0],A[1][1])
I3 = ((1,0),(0,1))
# generatory SL(2,3)
g1 = ((1,1),(0,1))
g2 = ((1,0),(1,1))
T2 = close_group([g1,g2], m3mul, I3, m3key)
r2 = analyze("2T = SL(2,3)  (target M = SO(3)/T)", T2, m3mul, I3, m3key,
             note="binarna tetraedralna, bogata trialnosc")

print("="*70)
print("PODSUMOWANIE POROWNANIA (substrat->wstazki):")
for r in (r1, r2):
    print(f"  {r['name'][:24]:24s} | klasy={r['classes']:2d} | ab=Z{r['abelianization']:<2d} | "
          f"kolor przezywa abelianizacje: {r['color_survives_abelianization']}")
