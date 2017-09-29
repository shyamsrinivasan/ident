from sympy import *

# experimental data
x11, x12, x13, x21, x22, x23 = symbols('x11, x12, x13, x21, x22, x23', positive=True)
v31, v32, v33 = symbols('v31, v32, v33')

# parameters
V3max, K3fdp, K3pep = symbols('V3max, K3fdp, K3pep', positive=True)
# ratio1, ratio2, ratio3, ratio4 = symbols('ratio1, ratio2, ratio3, ratio4',positive=True)

# flux equations
# regulation_inhibition = 1/(1 + pep_sat) # for future reference
reg1_term = 1/(1+K3pep/x11)
reg2_term = 1/(1+K3pep/x12)
reg3_term = 1/(1+K3pep/x13)
fdp_sat_1 = x21/K3fdp
fdp_sat_2 = x22/K3fdp
fdp_sat_3 = x23/K3fdp
flux_31 = v31 - (reg1_term*V3max*fdp_sat_1/(1+fdp_sat_1))
flux_32 = v32 - (reg2_term*V3max*fdp_sat_2/(1+fdp_sat_2))
flux_33 = v33 - (reg3_term*V3max*fdp_sat_3/(1+fdp_sat_3))

# nonlinear solution for parameters
system_flux = [flux_31, flux_32] # , flux_33]
parameters = [K3fdp, K3pep] # , K3fdp, K3pep]
variables = [x11, x21, v31, x12, x22, v32, x13, x23, v33]
solutions = nonlinsolve(system_flux, parameters)
solutions = list(solutions)
# print(len(solutions[0]))
print(solutions)

# function of expressions from above symbolic solutions
# func_v3max = lambdify(variables,solutions[0][0],"numpy")
# func_k3fdp = lambdify(variables,solutions[0][1],"numpy")
# func_k3pep = lambdify(variables,solutions[0][2],"numpy")
# print(func_v3max(10, 1, 2, .1))
