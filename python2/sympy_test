from sympy import *
from sympy import init_printing
import numpy
#from sympy.utilities.autowrap import autowrap
#from sympy.utilities.autowrap import ufuncify
x,y,z = symbols("x y z")
Eq(x+1,4)

a = (x+1)**2
b = x**2 + 2*x + 1
print(simplify(a-b))

c = x**2-2*x+1
print(simplify(a-c))

expr = cos(x)+1
expr2 = expr.subs(x,y)
print(expr2)

print(expr2.evalf(subs={y: 0}))

expr3 = sin(x)/x
f = lambdify(x, expr3)
print(f(3.14))
f2 = lambdify(y, expr2)
print(f2(3.14))

f3 = lambdify(y, expr2, 'numpy')

data = numpy.linspace(1,10,100)
print(f2(data))

# ufuncify does not work - do not know why
# f4 = ufuncify([x], expr, backend='f2py')
#print(f4(data))

# sympy functions
init_printing()
k, m, n = symbols('k m n')
facexp = factorial(n)
print(facexp)
print(binomial(n,k))

print(diff(y*exp(x**2), x,y))
