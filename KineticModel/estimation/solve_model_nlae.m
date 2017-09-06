function ymodel = solve_model_nlae(nlaefh,nc,npert,xinit,p_all,options)

[~,flux,oderhs,~,x,p,dfx] = nlaefh(nc,npert);
ode_dfx_fun = casadi.Function('ode_dfx_fun',{x,p},{oderhs,dfx});

nlae = @(xstate)ode_dfx_fun(xstate,p_all);
options = optimoptions('fsolve',options,'SpecifyObjectiveGradient','on');
[ymodel,fval,exitflag] = fsolve(nlae,xinit,options);
