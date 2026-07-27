import numpy as np, sys
from scipy.optimize import minimize
t=np.linspace(0.0,1.0,20001); st=np.sqrt(t)
f=np.where(t==0,1.0,np.arctan(np.where(t==0,1,st))/np.where(t==0,1,st))
T=np.stack([t,t**2,t**3],axis=1)
ship=np.array([1.1271711055988247,0.2849778513254418,0.008830042167832291,
               1.4605043e0,5.718157e-1,5.0166193e-2])
def err(v):
    N=1.0+T@v[:3]; D=1.0+T@v[3:]
    if np.any(D<=0): return 1.0
    return np.max(np.abs(N/D-f)/f)
print("shipped: %.4e (%.4f ulp)"%(err(ship),err(ship)/2**-24)); sys.stdout.flush()
best=ship.copy(); bv=err(ship)
scale=np.abs(ship)*1e-4
for it in range(12):
    r=minimize(lambda z: err(best+z*scale), np.zeros(6), method='Nelder-Mead',
               options={'maxiter':40000,'maxfev':40000,'xatol':1e-10,'fatol':1e-18})
    cand=best+r.x*scale
    if err(cand)<bv: bv=err(cand); best=cand
    scale=scale*0.5
print("refit  : %.4e (%.4f ulp)"%(bv,bv/2**-24))
print("headroom: %.2fx"%(err(ship)/bv))
b32=np.float32(best)
N=1.0+T@np.float64(b32[:3]); D=1.0+T@np.float64(b32[3:])
e32=np.max(np.abs(N/D-f)/f)
print("f32-rounded refit: %.4e (%.4f ulp)"%(e32,e32/2**-24))
print("coeffs a0,a1,a2,b0,b1,b2 =",[float(x) for x in best])
