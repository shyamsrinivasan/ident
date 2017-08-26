% or load pre-calculated data
load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noisy_model/pdata_aug24');

nc = 3;
nf = 6;
npert = 4;
np = 13;
nvar = nc*npert+np+nf*npert+2;

xexp = exp_sol{1}.xss;
vexp = exp_sol{1}.fss;
xexp = reshape(xexp,[nc*npert,1]);
vexp = reshape(vexp,[nf*npert,1]);

% sym obj
[obj,var,par,x,p,flux,vareps,p_usl,ac,wts] =...
kotte_pest_allf_obj(xexp,vexp,nc,nf,npert);
% sym cons
cons =...
kotte_pest_allf_cons(xexp,vexp,nc,nf,npert,x,p,flux,vareps,p_usl,ac);
ncons = size(cons,1);

% get parameters
p0_obj = optimdata.odep(1:13)';
scale = ones(np,1);
scale(3) = 1e-6;

weigths = [1000;1;1;1000];
p0_const = optimdata.odep(14:end)';
par_val = [p0_const;weigths];

% get initial value
x0 = [xexp;p0_obj.*scale;vexp;0;0];
% set bounds for x0
lb = zeros(nvar,1);
ub = zeros(nvar,1);
lb(1:nc*npert) = 1e-7;
ub(1:nc*npert) = 300;
lb(nc*npert+1:nc*npert+np) = .008;
ub(nc*npert+1:nc*npert+np) = 3;
lb(nc*npert+np+1:nc*npert+np+nf*npert) = 0;
ub(nc*npert+np+1:nc*npert+np+nf*npert) = 100;
lb(nvar-2+1:nvar) = 0;
ub(nvar-2+1:nvar) = 1;
% set bounds for cons
lbg = zeros(ncons,1);
ubg = zeros(ncons,1);


estprob = struct('x',var,'f',obj,'g',cons);
sol = estprob('x0',x0,'p',par_val,'lbx',lb,'ubx',ub,'lbg',lbg,'ubg',ubg);


x = casadi.SX.sym('x',1);
y = casadi.SX.sym('y',1);
nlp = struct('x',[x;y],'f',x/y);
options.warn_initial_bounds = 1;
NLP = casadi.nlpsol('solver','ipopt',nlp,options);
sol = NLP('x0',[.1 .1],'lbx',[.001 .001],'ubx',[10 10]);
% int = struct('x',x,'ode',exp(-x));
% opt.grid = 0:.1:20;
% F = casadi.integrator('F','cvodes',int,opts);




