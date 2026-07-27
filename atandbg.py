import numpy as np
from scipy.optimize import linprog
t = np.linspace(0.0,1.0,6001); st=np.sqrt(t)
f = np.where(t==0,1.0,np.arctan(np.where(t==0,1,st))/np.where(t==0,1,st))
ship=np.array([1.1271711055988247,0.2849778513254418,0.008830042167832291,
               1.4605043e0,5.718157e-1,5.0166193e-2])
P=np.stack([t,t**2,t**3],axis=1)
def build(eps):
    A1=np.hstack([-P,(f*(1-eps))[:,None]*P]); b1=1.0-f*(1-eps)
    A2=np.hstack([P,-(f*(1+eps))[:,None]*P]); b2=f*(1+eps)-1.0
    return np.vstack([A1,A2]),np.concatenate([b1,b2])
for eps in (2e-9,1.95e-9,1.9e-9,1e-9):
    Aub,bub=build(eps)
    viol=np.max(Aub@ship-bub)
    r=linprog(np.zeros(6),A_ub=Aub,b_ub=bub,bounds=[(None,None)]*6,method='highs')
    print("eps %.2e  shipped max violation %.3e  lp status %d %s"%(eps,viol,r.status,r.message[:40]))
