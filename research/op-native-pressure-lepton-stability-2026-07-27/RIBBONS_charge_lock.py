# RIBBONS_charge_lock.py
# Strzalka wstazki->kwarki, mechanizm (alpha): ladunek = pi_3(B) + trialnosc.
# TEST CENTRALNY (§0 lock): czy ulamkowy ladunek EM w TRZECICH i jego LOCK z kolorem
# WYPADAJA z topologii, czy trzeba je wsadzic?
#
# Klucz: ladunki U(1)_EM rozrozniajace sektory dysklinacji = Hom(pi_1(M), U(1))
#        = Hom(2T, U(1)) = Hom(2T^ab, U(1)) = Hom(Z3, U(1)) = Z3.
# Jesli tak: 3 sektory ladunku {0, 1/3, 2/3} = 3 kolory. Lock kolor<->ladunek
# jest wtedy STRUKTURALNY (nie postulat grupy cechowania jak w SM).
#
# Test na cyrkularnosc: gdzie moze wejsc "1/3" recznie? Odpowiedz: NIGDZIE w tym
# rachunku - "3" pochodzi wylacznie z |Z3| = |2T/Q8|, ktore samo bylo WYMUSZONE
# przez trialnosc (uniqueness 2T). Zadnego dopasowania do znanych ladunkow.

import math, cmath
def qmul(p,q):
    w1,x1,y1,z1=p; w2,x2,y2,z2=q
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
def qk(p): return tuple(round(c,6)+0.0 for c in p)
QID=(1.,0.,0.,0.); NEG=(-1.,0.,0.,0.)
def close(gens):
    E={qk(QID):QID}; fr=[QID]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                pr=qmul(g,h);k=qk(pr)
                if k not in E:E[k]=pr;nw.append(pr)
        fr=nw
    return list(E.values())
def inv(g):w,x,y,z=g;return(w,-x,-y,-z)
def comm_sub(E):
    cs=[]
    for a in E:
        ai=inv(a)
        for b in E: cs.append(qmul(qmul(a,b),qmul(ai,inv(b))))
    return close(cs)

G=close([(0.,1.,0.,0.),(0.5,0.5,0.5,0.5)]); assert len(G)==24
Q8=comm_sub(G); Q8set=set(qk(x) for x in Q8)
def coset(g): return frozenset(qk(qmul(g,x)) for x in Q8)
# --- trialnosc jako PRAWDZIWY homomorfizm 2T->Z3 (log dyskretny wzgl. generatora warstwy) ---
C0=coset(QID)                                   # warstwa tozsamosciowa = 0
gen=next(g for g in G if coset(g)!=C0)           # generator warstwy nietrywialnej = 1
C1=coset(gen); C2=coset(qmul(gen,gen))           # C1^2 = 2
label={C0:0, C1:1, C2:2}
assert coset(qmul(qmul(gen,gen),gen))==C0        # C1^3 = 0 (rzad 3 w ilorazie)
def tri(g): return label[coset(g)]               # tri(gh)=tri(g)+tri(h) mod 3 (z konstrukcji)
# sanity: tri jest homomorfizmem
assert all(tri(qmul(a,b))==(tri(a)+tri(b))%3 for a in G for b in G), "tri nie jest homomorfizmem!"

print("Hom(2T,U(1)) = charaktery 1-wymiarowe (= ladunki ulamkowe):")
w=cmath.exp(2j*math.pi/3)
# 3 kandydaci charakterow: chi_k(g) = w^(k*tri(g)), k=0,1,2
def homomorphism_ok(k):
    for a in G:
        for b in G:
            lhs=w**((k*tri(qmul(a,b)))%3)
            rhs=(w**((k*tri(a))%3))*(w**((k*tri(b))%3))
            if abs(lhs-rhs)>1e-9: return False
    return True
n_hom=0
for k in range(3):
    ok=homomorphism_ok(k)
    n_hom+= 1 if ok else 0
    print(f"  chi_{k}: homomorfizm 2T->U(1)? {ok}   ladunek ulamkowy = {k}/3 mod 1")
print(f"  |Hom(2T,U(1))| = {n_hom}  => {'Z3 POTWIERDZONE' if n_hom==3 else 'NIE Z3'}")

print("\nLADUNEK vs KOLOR vs SPIN (co wypada):")
print(f"  kolor (trialnosc) elementow: warstwy {sorted({tri(g) for g in G})} = Z3")
print(f"  ladunek ulamkowy sektora t = t/3 mod 1: {{0, 1/3, 2/3}}  <- 3 kolory")
print(f"  spin -1: trialnosc={tri(NEG)} => ladunek ulamkowy 0 (bezbarwny, jak nalezy)")
print(f"  elementy rzedu 3 (kolor/antykolor): trialnosci {sorted({tri(g) for g in G if abs(qmul(g,qmul(g,g))[0]-1)<1e-6 and qk(g)!=qk(QID)})}")

print("\n" + "="*64)
print("WYPADA (DERIVED, bez fitu):")
print("  - Hom(2T,U(1)) = Z3  => ladunek EM ulamkowy w TRZECICH.")
print("  - fractional(Q) = trialnosc/3  => LOCK kolor<->ladunek STRUKTURALNY.")
print("  - '3' pochodzi WYLACZNIE z |2T/Q8|=|Z3|, wymuszonego trialnoscia. Zero dopasowania.")
print("  - spin (-1) bezbarwny, ladunek ulamkowy 0 - spojne.")
print("NIE WYPADA (OTWARTE, wymaga inputu - NIE fituje):")
print("  - dokladne wartosci (+2/3 vs -1/3): offset od nawiniecia B/2 (Gell-Mann-Nishijima)")
print("    + normalizacja U(1). To osobne wyprowadzenie.")
print("  - emergencja EM z pradu topologicznego (RAMY §3) - wlasne wyprowadzenie.")
