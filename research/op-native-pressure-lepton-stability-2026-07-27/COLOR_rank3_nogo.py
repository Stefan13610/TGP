"""
POPRAWKA testu [3] z COLOR_z3_sources_check.py.
Tam uzylem IZOTROPOWEGO C(a) -> zerowal sie tez rank-2, wiec test nic nie rozstrzygal.

Wlasciwy test. Jedyny wiez na korelator wiazaniowy to TRANSLACYJNA NIEZMIENNICZOSC:
    C(-a) = <s_i s_{i-a}> = <s_{i+a} s_i> = <s_i s_{i+a}> = C(a)      [exact]
Zaden dodatkowy postulat o symetrii sieci nie jest potrzebny.

Pytanie: przy DOWOLNYM anizotropowym C spelniajacym C(-a)=C(a):
    rank-2  M_ij  = sum_a C(a) a_i a_j          -- czy niezerowy?
    rank-3  T_ijk = sum_a C(a) a_i a_j a_k      -- czy niezerowy?
"""
import itertools, numpy as np
rng = np.random.default_rng(7)

def shells(kind):
    if kind == "6":  return [v for v in itertools.product([-1,0,1],repeat=3) if sum(map(abs,v))==1]
    if kind == "18": return [v for v in itertools.product([-1,0,1],repeat=3) if 1<=sum(map(abs,v))<=2]
    if kind == "26": return [v for v in itertools.product([-1,0,1],repeat=3) if v!=(0,0,0)]
    if kind == "124":return [v for v in itertools.product([-2,-1,0,1,2],repeat=3) if v!=(0,0,0)]

print("="*74)
print("NO-GO: rank-3 z korelatora 2-punktowego (dowolne anizotropowe C, C(-a)=C(a))")
print("="*74)
print(f"{'sasiedztwo':<14}{'#wiazan':>9}{'max|M_ij^TF| (rank2)':>24}{'max|T_ijk| (rank3)':>22}")

for kind in ["6","18","26","124"]:
    nb = shells(kind)
    # losowe ANIZOTROPOWE C, ale wymuszamy C(-a)=C(a) (translacyjna niezmienniczosc)
    C = {}
    for a in nb:
        na = tuple(-x for x in a)
        if na in C: C[a] = C[na]
        else:       C[a] = rng.normal()
    M = np.zeros((3,3)); T = np.zeros((3,3,3))
    for a,c in C.items():
        av = np.array(a,dtype=float)
        M += c*np.outer(av,av)
        T += c*np.einsum('i,j,k->ijk',av,av,av)
    Mtf = M - np.trace(M)/3*np.eye(3)
    print(f"{kind:<14}{len(nb):>9}{np.abs(Mtf).max():>24.6f}{np.abs(T).max():>22.3e}")

print("\nWniosek: rank-2 (sigma_ab) jest GENERYCZNIE NIEZEROWY -- nematyk istnieje.")
print("         rank-3 znika TOZSAMOSCIOWO, dla kazdego sasiedztwa i kazdego C.")
print("         Powod: T_ijk jest NIEPARZYSTY w a, a C jest PARZYSTY => sumowanie po")
print("         parach (a,-a) kasuje sie skladnik po skladniku. To nie jest artefakt siatki.")

print("\n" + "="*74)
print("Czy 3-punktowy korelator ratuje sytuacje? (Z_2: s -> -s)")
print("="*74)
print("  <s s_+a s_+b> ma 3 czynniki s  =>  NIEPARZYSTY  =>  = 0 w fazie Z_2-symetrycznej.")
print("  W fazie zlamanej (v=<s> != 0, T<T_c) rozklad s = v + ds daje:")
print("    <s s_+a s_+b> = v^3 + v*[<ds ds> x 3 pary] + <ds ds ds>_c")
print("  - czlon v^3            : staly, bezkierunkowy (rank 0)")
print("  - czlon v*<ds ds>      : REDUKOWALNY do sigma_ab (nic nowego)")
print("  - czlon <ds ds ds>_c   : wymaga sprzezenia KUBICZNEGO w H_Gamma (lam3 s^3),")
print("                           ktore jest ZABRONIONE przez Z_2.  => = 0")
print("\n  ==> Zaden nowy niezalezny parametr porzadku rangi 3 nie powstaje.")

print("\n" + "="*74)
print("Jedyna ocalala droga: korelator 4-punktowy <s s_+a s_+b s_+c> (parzysty w s)")
print("="*74)
print("  Z_2-parzysty (4 czynniki)  => NIE znika.  Trzy niezalezne przesuniecia a,b,c")
print("  => moze nosic strukture rank-3 w indeksach przestrzennych.")
print("  ALE: (i) nie wystepuje w rdzeniu TGP jako obserwabla poziomu 0,")
print("       (ii) w teorii Gaussa/MFT rozklada sie na pary (tw. Wicka) => redukowalny,")
print("       (iii) czesc spojna <....>_c jest niezerowa tylko poza MFT (lam0 s^4).")
print("  Status: OTWARTE, ale to bylby NOWY POSTULAT o obserwabli, nie wniosek.")
