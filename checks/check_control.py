"""Small exact check of the Heisenberg cocycle used as a p=3 control.

This is a consistency check, not a classification or data-generation job.
The general cocycle identity is proved by the cup-cochain calculation in notes.
"""
from itertools import product

p = 3
inv2 = pow(2, -1, p)
elts = list(product(range(p), repeat=3))

def mul(g, h):
    x, y, z = g
    u, v, w = h
    return ((x+u) % p, (y+v) % p, (z+w+x*v) % p)

def phi(g, h, k):
    x, y, z = g
    u, v, w = h
    a, b, c = k
    return (z*u*b + inv2*x*x*v*b + x*u*v*b) % p

mult = [[elts.index(mul(g,h)) for h in elts] for g in elts]
f = [[[phi(g,h,k) for k in elts] for h in elts] for g in elts]
for i,j,k,l in product(range(len(elts)), repeat=4):
    residue = (f[j][k][l] - f[mult[i][j]][k][l]
               + f[i][mult[j][k]][l] - f[i][j][mult[k][l]]
               + f[i][j][k]) % p
    assert residue == 0, (i,j,k,l,residue)

# Central slants: alpha_a(g,h)=phi(a,g,h)+phi(g,h,a)-phi(g,a,h).
# The trivializing exponent is t_a(g)=-a*z_g.
for a,g,h in product(range(p), elts, elts):
    z_a = (0,0,a)
    alpha = (phi(z_a,g,h)+phi(g,h,z_a)-phi(g,z_a,h)) % p
    tdelta = (-a*g[2]-a*h[2]+a*mul(g,h)[2]) % p
    assert alpha == tdelta

# The tensor correction on central degrees a,b is identically zero.
for a,b,g in product(range(p), range(p), elts):
    za,zb=(0,0,a),(0,0,b)
    correction=(phi(za,zb,g)+phi(g,za,zb)-phi(za,g,zb))%p
    assert correction==0

print('p=3: all 27^4 cocycle equations, central slants and tensor corrections pass')
print('central metric q(a) = exp(-2*pi*i*a^2/3); complement q(a,b) = exp(-2*pi*i*(a^2+b^2)/3)')
