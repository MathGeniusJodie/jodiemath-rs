import sympy as sp

x, y = sp.symbols('x y', positive=True)

inv = y**3   # inverse of g(x) = x**(1/3)
f = 1 / (inv - x)

fp = [f]
for _ in range(6):
    fp.append(sp.diff(fp[-1], y))

newton       = sp.simplify(y + fp[0] / fp[1])
halley       = sp.simplify(y + 2 * fp[1] / fp[2])
householder3 = sp.simplify(y + 3 * fp[2] / fp[3])
householder4 = sp.simplify(y + 4 * fp[3] / fp[4])
householder5 = sp.simplify(y + 5 * fp[4] / fp[5])
householder6 = sp.simplify(y + 6 * fp[5] / fp[6])
