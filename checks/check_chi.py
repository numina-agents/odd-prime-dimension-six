"""Small independent identities for the invariant Heisenberg class."""
from itertools import product
import numpy as np

# p=3: cyclic restrictions of the corrected cohomology class.
p=3
G=list(product(range(p),repeat=2))
f=lambda x,y:(-x*x*y*y+x*x+y*y)%p
assert all(f(x,y)==int((x,y)!=(0,0)) for x,y in G)
for a,b,c,d in product(range(p),repeat=4):
    det=(a*d-b*c)%p
    if not det:continue
    for x,y in G:
        assert f((a*x+b*y)%p,(c*x+d*y)%p)==det*det*f(x,y)%p
print('p=3: corrected cyclic restrictions are GL2-invariant on all 9 vectors')

# p=5: the universal Pluecker primitive, without the Heisenberg z coordinate.
p=5
G=np.array(list(product(range(p),repeat=2)),dtype=np.int64)
N=len(G)
E=(G[:,0,None]*G[None,:,1]-G[:,1,None]*G[None,:,0])%p
idx={tuple(v):i for i,v in enumerate(G)}
add=np.array([[idx[tuple((v+w)%p)] for w in G] for v in G])
F1=E[:,:,None]*E[None,:,:]
F2=E[:,:,None]*E[:,None,:]
T=(2*F1+F2)*pow(12,-1,p)%p
j=np.arange(N)[:,None,None]
k=np.arange(N)[None,:,None]
l=np.arange(N)[None,None,:]
for i in range(N):
    dt=(T[j,k,l]-T[add[i,:,None,None],k,l]+T[i,add[:,:,None],l]
        -T[i,j,add[None,:,:]]+T[i,j,k])%p
    target=E[i,:,None,None]*E[None,:,:]*pow(4,-1,p)%p
    assert np.array_equal(dt,target)
print('p=5: all 25^4 identities for the invariant Pluecker primitive pass')
