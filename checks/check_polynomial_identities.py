"""Independent characteristic-zero polynomial checks of the two cocycles.

Vanishing over Q implies the asserted reductions in every allowed characteristic.
This does not test the classification or replace its equivalence argument.
"""
from sympy import Rational, expand, symbols

g, h, k, l = [symbols(f'{name}x {name}y {name}z') for name in 'ghkl']

def det(a, b):
    return a[0] * b[1] - a[1] * b[0]

def original(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2]+a[0]*b[1])

def bch(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2]+det(a,b)/2)

def phi(a, b, c):
    return a[2]*b[0]*c[1] + Rational(1,2)*a[0]**2*b[1]*c[1] + a[0]*b[0]*b[1]*c[1]

def chi(a, b, c):
    return Rational(1,2)*a[2]*det(b,c) + Rational(1,12)*(2*det(a,b)*det(b,c)+det(a,b)*det(a,c))

def differential(f, mul):
    return expand(f(h,k,l)-f(mul(g,h),k,l)+f(g,mul(h,k),l)-f(g,h,mul(k,l))+f(g,h,k))

assert differential(phi, original) == 0
assert differential(chi, bch) == 0

a = symbols('a')
central = (0, 0, a)
for name, f, mul in [('phi', phi, original), ('chi', chi, bch)]:
    slant = f(central,g,h)+f(g,h,central)-f(g,central,h)
    dtau = -a*g[2]-a*h[2]+a*mul(g,h)[2]
    assert expand(slant-dtau) == 0
    print(name + ': cocycle and central-slant identities vanish over Q')
