import numpy as np
from scipy.optimize import linprog
import sys

t = np.linspace(0.0, 1.0, 6001)
st = np.sqrt(t)
f = np.where(t==0, 1.0, np.arctan(np.where(t==0,1,st))/np.where(t==0,1,st))

a=[0.008830042167832291,0.2849778513254418,1.1271711055988247]  # a2,a1,a0
b=[5.0166193e-2,5.718157e-1,1.4605043e0]                        # b2,b1,b0
# code: numer/x = ((a2*t+a1)*t+a0)*t+1 -> N = 1 + a0 t + a1 t^2 + a2 t^3
A=[a[2],a[1],a[0]]; B=[b[2],b[1],b[0]]
def poly(c,t): return 1.0 + c[0]*t + c[1]*t**2 + c[2]*t**3
cur = np.max(np.abs(poly(A,t)/poly(B,t) - f)/f)
print("shipped idealized max rel err: %.4e  (%.3f ulp-equiv)"%(cur,cur/2**-24)); sys.stdout.flush()

def feasible(eps):
    # vars: a0,a1,a2,b0,b1,b2
    P = np.stack([t,t**2,t**3],axis=1)
    n=len(t)
    A1 = np.hstack([-P, (f*(1-eps))[:,None]*P]); b1 = 1.0 - f*(1-eps)
    A2 = np.hstack([ P, -(f*(1+eps))[:,None]*P]); b2 = f*(1+eps) - 1.0
    Aub=np.vstack([A1,A2]); bub=np.concatenate([b1,b2])
    r=linprog(np.zeros(6),A_ub=Aub,b_ub=bub,bounds=[(None,None)]*6,method='highs')
    return r.status==0, (r.x if r.status==0 else None)
lo,hi=1e-10,1e-3; best=None
for _ in range(60):
    mid=np.sqrt(lo*hi)
    ok,x=feasible(mid)
    if ok: hi=mid; best=x
    else: lo=mid
    if hi/lo<1.0001: break
print("optimal [3/3] idealized max rel err: %.4e  (%.3f ulp-equiv)"%(hi,hi/2**-24))
print("headroom factor: %.2fx"%(cur/hi))
if best is not None:
    A2c=[best[0],best[1],best[2]]; B2c=[best[3],best[4],best[5]]
    print("N coeffs (a0,a1,a2):",A2c); print("D coeffs (b0,b1,b2):",B2c)
    e=np.max(np.abs(poly(A2c,t)/poly(B2c,t)-f)/f); print("verify %.4e"%e)
    A32=np.float32(A2c); B32=np.float32(B2c)
    e32=np.max(np.abs(poly(np.float64(A32),t)/poly(np.float64(B32),t)-f)/f)
    print("f32-rounded: %.4e (%.3f ulp-equiv)"%(e32,e32/2**-24))
