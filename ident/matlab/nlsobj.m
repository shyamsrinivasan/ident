function [objfx,Dobjfx,D2objfx] = nlsobj(nvar,ntime)

yexp = casadi.SX('yexp',nvar,ntime);
ymodel = casadi.SX('ymodel',nvar,ntime);
yexp_var = casadi.SX('yexp_var',nvar,ntime);

obj = sum(sum(((yexp-ymodel)./yexp_var).^2,2)); 
% ymodel has to be expressed as casadi.SX expression here to use parameters
% theta here

objfx = casadi.Function('objfx',{yexp,ymodel,yexp_var},{obj});

Dobjfx = casadi.Function('Dobjfx',{yexp,ymodel,yexp_var},{jacobian()});


