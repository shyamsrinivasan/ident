% or load pre-calculated data
load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noisy_model/pdata_aug24');

[obj,var,par,x,p,vareps,p_usl,ac,flux,wts] =...
kotte_pest_allf_obj(exp_sol{1}.xss,exp_sol{1}.fss,3,6,4);
cons =...
kotte_pest_allf_cons(exp_sol{1}.xss,exp_sol{1}.fss,x,p,p_usl,ac,flux,vareps);

estprob = struct('x',var,'f',obj,'g',cons);
sol = estprob('x0',x0,'p',p,'lbx',lb,'ubx',ub);


x = casadi.SX.sym('x',1);
y = casadi.SX.sym('y',1);
nlp = struct('x',[x;y],'f',x/y);
options.warn_initial_bounds = 1;
NLP = casadi.nlpsol('solver','ipopt',nlp,options);
sol = NLP('x0',[.1 .1],'lbx',[.001 .001],'ubx',[10 10]);
% int = struct('x',x,'ode',exp(-x));
% opt.grid = 0:.1:20;
% F = casadi.integrator('F','cvodes',int,opts);




