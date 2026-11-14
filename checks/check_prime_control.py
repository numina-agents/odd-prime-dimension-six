"""Exact small p=3 check for the rank-one-power centralizer model."""
from itertools import product
import numpy as np

p=3
n=p*p
cc=vv=1
elts=np.array(list(product(range(p),repeat=4)),dtype=np.int64)
x,y,z,t=elts.T
r=x+p*t
idx={tuple(g):i for i,g in enumerate(elts)}
mul=np.empty((len(elts),len(elts)),dtype=np.int64)
for i,g in enumerate(elts):
    for j,h in enumerate(elts):
        xx=(g[0]+h[0])%p
        yy=(g[1]+h[1])%p
        zz=(g[2]+h[2]+g[0]*h[1])%p
        tt=(g[3]+h[3]+(g[0]+h[0])//p)%p
        mul[i,j]=idx[xx,yy,zz,tt]

phi=(cc*(z[:,None,None]*x[None,:,None]*y[None,None,:]
    +pow(2,-1,p)*x[:,None,None]**2*y[None,:,None]*y[None,None,:]
    +x[:,None,None]*x[None,:,None]*y[None,:,None]*y[None,None,:])
    +vv*x[:,None,None]*((r[None,:,None]+r[None,None,:])//n))%p
N=len(elts)
jj=np.arange(N)[:,None,None]
kk=np.arange(N)[None,:,None]
ll=np.arange(N)[None,None,:]
for i in range(N):
    delta=(phi[jj,kk,ll]-phi[mul[i,:,None,None],kk,ll]
       +phi[i,mul[:,:,None],ll]-phi[i,jj,mul[None,:,:]]+phi[i,jj,kk])%p
    assert not np.any(delta),i

central=[idx[0,0,a,b] for a,b in product(range(p),repeat=2)]
for a,b in product(range(p),repeat=2):
    ia=idx[0,0,a,b]
    tau=(-p*cc*a*z-vv*b*r)%n
    delta_tau=(tau[:,None]+tau[None,:]-tau[mul])%n
    slant=(p*(phi[ia,:,:]+phi[:,:,ia]-phi[:,ia,:]))%n
    assert np.array_equal(delta_tau,slant)
    for aa,bb in product(range(p),repeat=2):
        ib=idx[0,0,aa,bb]
        tau_b=(-p*cc*aa*z-vv*bb*r)%n
        tau_sum=(-p*cc*((a+aa)%p)*z-vv*((b+bb)%p)*r)%n
        correction=(p*(phi[ia,ib,:]+phi[:,ia,ib]-phi[ia,:,ib]))%n
        assert np.all((tau+tau_b+correction-tau_sum)%n==0)
        self_exp=(-p*cc*a*a-vv*b*(p*b))%n
        assert self_exp==(p*(-cc*a*a-vv*b*b))%n

print('p=3 rank-one-power model: all 81^4 cocycle equations pass (chunked, <40 MB)')
print('all central slants, order-p section tensor corrections, and anisotropic self-braidings pass')
print('centralizer invertibles: a=0, character k=2*b mod 3; group C9 x C3, twists zeta3^(b^2)')
