from sympy import *

# experimental data
x11, x12 = symbols('x11, x12', positive=True)
v11, v12 = symbols('v11, v12', positive=True)
ac1, ac2 = symbols('ac1, ac2', positivie=True)

# parameters
K1ac, k1cat = symbols('K1ac, k1cat', positive=True)
flux_11 = v11 - k1cat*ac1/(K1ac + ac1)
flux_12 = v12 - k1cat*ac2/(K1ac + ac2)

# nonlinear solution for parameters
system_flux = [flux_11, flux_12]
parameters = [K1ac, k1cat]
variables = [ac1, v11, ac2, v12]
solutions = nonlinsolve(system_flux, parameters)
solutions = list(solutions)
# print(len(solutions[0]))

# function of expressions from above symbolic solutions
func_k1ac = lambdify(variables,solutions[0][0],"numpy")
func_k1cat = lambdify(variables,solutions[0][1],"numpy")
print(func_k1ac(10, 1, 2, .1))
