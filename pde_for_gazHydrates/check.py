import sympy



t, s = sympy.symbols('t, s')
a = sympy.symbols('a', real=True, positive=True)

f = sympy.exp(-a*t)
f = sympy.laplace
print(f)