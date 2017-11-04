from sympy import *

# experimental data
x11, x12, x13, x14, x21, x22, x23, x24 = symbols('x11, x12, x13, x14, x21, x22, x23, x24', positive=True)
v31, v32, v33, v34 = symbols('v31, v32, v33, v34')

# parameters
V3max, K3fdp, K3pep, L3fdp = symbols('V3max, K3fdp, K3pep, L3fdp', positive=True)
# ratio1, ratio2, ratio3, ratio4 = symbols('ratio1, ratio2, ratio3, ratio4',positive=True)

# flux equations
reg1_term = (1 + x11/K3pep)**(-4)
reg2_term = (1 + x12/K3pep)**(-4)
reg3_term = (1 + x13/K3pep)**(-4)
reg4_term = (1 + x14/K3pep)**(-4)
flux_31 = v31 - (V3max * (1 + x21/K3fdp) * (x21/K3fdp)**3)/((1 + x21/K3fdp)**4+L3fdp*reg1_term)
flux_32 = v32 - (V3max * (1 + x22/K3fdp) * (x22/K3fdp)**3)/((1 + x22/K3fdp)**4+L3fdp*reg2_term)
flux_33 = v33 - (V3max * (1 + x23/K3fdp) * (x23/K3fdp)**3)/((1 + x23/K3fdp)**4+L3fdp*reg3_term)
flux_34 = v34 - (V3max * (1 + x24/K3fdp) * (x24/K3fdp)**3)/((1 + x24/K3fdp)**4+L3fdp*reg4_term)

# nonlinear solution for parameters
system_flux = [flux_31, flux_32, flux_33]
parameters = [V3max, K3fdp, K3pep]
variables = [x11, x21, v31, x12, x22, v32, x13, x23, v33]
solutions = nonlinsolve(system_flux, parameters)
solutions = list(solutions)
# print(len(solutions[0]))

# function of expressions from above symbolic solutions
func_v3max = lambdify(variables,solutions[0][0],"numpy")
func_k3fdp = lambdify(variables,solutions[0][1],"numpy")
func_k3pep = lambdify(variables,solutions[0][2],"numpy")
print(func_v3max(10, 1, 2, .1))
