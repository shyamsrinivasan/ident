from sympy import *

# experimental data
x11, x12 = symbols('x11, x12', positive=True)
v21, v22 = symbols('v21, v22', positive=True)

# parameters
K2pep, V2max = symbols('K2pep, V2max', positive=True)
flux_21 = v21 - V2max * x11 / (K2pep + x11)
flux_22 = v22 - V2max * x12 / (K2pep + x12)

# nonlinear solution for parameters
system_flux = [flux_21, flux_22]
parameters = [K2pep, V2max]
variables = [x11, v21, x12, v22]
solutions = nonlinsolve(system_flux, parameters)
solutions = list(solutions)
# print(len(solutions[0]))

# function of expressions from above symbolic solutions
func_k2pep = lambdify(variables,solutions[0][0],"numpy")
func_v2max = lambdify(variables,solutions[0][1],"numpy")
print(func_k2pep(10, 1, 2, .1))