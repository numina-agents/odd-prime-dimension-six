"""Cyclic-restriction regression checks for both explicit dimension-p^6 models.

No imports from the other checks. Evaluate actual cocycles on (g,g^j,g),
not just the claimed quadratic restriction formula. All arithmetic is exact.
The paper supplies the proof for arbitrary primes; these are small controls.
"""
from itertools import product
from collections import Counter
from math import gcd, lcm


def nonsquare(p):
    return next(a for a in range(1, p) if pow(a, (p-1)//2, p) == p-1)


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def phi(a, b, c, p):
    return (a[2]*b[0]*c[1] + pow(2,-1,p)*a[0]**2*b[1]*c[1]
            + a[0]*b[0]*b[1]*c[1]) % p


def chi(a, b, c, p):
    if p == 3:
        return (phi(a,b,c,p) + a[0]*((b[0]+c[0])//p)
                + a[1]*((b[1]+c[1])//p)) % p
    e1,e2,e3=det(a,b),det(b,c),det(a,c)
    return (pow(2,-1,p)*a[2]*e2
            + pow(12,-1,p)*(2*e1*e2+e1*e3)) % p


def hp_mul(a,b,p):
    correction = a[0]*b[1] if p == 3 else pow(2,-1,p)*det(a,b)
    return ((a[0]+b[0])%p,(a[1]+b[1])%p,(a[2]+b[2]+correction)%p)


def forms(p,e):
    return [(0,0,0),(1,0,0),(e,0,0)]+[(1,0,d) for d in range(1,p)]


for p in (3,5,7,11):
    e=nonsquare(p)
    restrictions=[]
    for g in product(range(p),repeat=3):
        h=(0,0,0); hsum=0; carries=[0,0,0]
        for _ in range(p):
            hsum += chi(g,h,g,p)
            carries[0] += g[0]*((h[0]+g[0])//p)
            carries[1] += g[0]*((h[1]+g[1])//p)
            carries[2] += g[1]*((h[1]+g[1])//p)
            h=hp_mul(h,g,p)
        assert h==(0,0,0)
        hsum%=p; carries=[x%p for x in carries]
        assert carries==[g[0]**2%p,g[0]*g[1]%p,g[1]**2%p]
        assert hsum==(int(g[:2]!=(0,0)) if p==3 else 0)
        restrictions.append((hsum,carries))
    counts=Counter()
    for c in (1,e):
        for q in forms(p,e):
            nonzero=any((c*h+sum(a*b for a,b in zip(q,cs)))%p
                        for h,cs in restrictions)
            centre_exponent=p*p if nonzero else p
            expected=p if p>3 and q==(0,0,0) else p*p
            # B has exponent p, and C_0 contains a pointed factor of exponent p.
            assert centre_exponent==expected,(p,c,q,centre_exponent,expected)
            counts[centre_exponent]+=1
    print(f'p={p}: product normal forms by exact exponent {dict(sorted(counts.items()))}')


def gamma_mul(a,b,p):
    return ((a[0]+b[0])%p,(a[1]+b[1])%p,
            (a[2]+b[2]+a[0]*b[1])%p,
            (a[3]+b[3]+(a[0]+b[0])//p)%p)


def prime_omega(a,b,c,p,cc,v,s):
    rb=b[0]+p*b[3];rc=c[0]+p*c[3]
    return (cc*phi(a,b,c,p)+v*a[0]*((rb+rc)//(p*p))
            +s*a[1]*((b[1]+c[1])//p))%p


for p in (3,5):
    e=nonsquare(p)
    for cc in (1,e):
        v=(-e*pow(cc,-1,p))%p
        for s in range(p):
            exponent=1
            for g in product(range(p),repeat=4):
                order=1 if g==(0,0,0,0) else (p*p if g[0] else p)
                h=(0,0,0,0);total=0
                for _ in range(order):
                    total+=prime_omega(g,h,g,p,cc,v,s)
                    h=gamma_mul(h,g,p)
                assert h==(0,0,0,0)
                cocycle_order=p if total%p else 1
                exponent=lcm(exponent,order*cocycle_order)
            assert exponent==p**3,(p,cc,s,exponent)
    print(f'p={p}: all {2*p} prime normal forms have exponent {p**3}')
print('Cyclic exponent regressions pass; p=5,Q=0 has exponent 5, not 25.')
