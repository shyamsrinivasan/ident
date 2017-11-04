% parameter optimization based on sysid example in casadi examples pack
% get non-noisy experimental data
% generate experimental data - get initial ss
tspan = 0:0.1:100;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise_dyn(opts,xss,fss,odep_bkp,pt_sol_id);
close all

p_truth = opts.odep(2:end)';
p_guess = opts.odep(2:end)';

% cas model
x = casadi.SX.sym('x',3,1);
p = casadi.SX.sym('p',16,1);

% K1ac = p(1);    % or 0.02
K3fdp = p(1);
L3fdp = p(2);
K3pep = p(3);
K2pep = p(4);
vemax = p(5);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(6);        % or 0.45
ne = p(7);             % or 2
d = p(8);
V4max = p(9);
k1cat = p(10);   
V3max = p(11);    
V2max = p(12);   
ac = p(16);
K1ac = .1;

flux1 = k1cat.*x(3).*ac./(ac+K1ac);
flux2 = vemax.*(1-1./(1+(KeFDP./x(2)).^ne));
ratio = 1+x(2)./K3fdp;
flux3 = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1)./K3pep).^(-4));
flux4 = V2max.*x(1)./(x(1)+K2pep);
flux5 = V4max.*x(1);
flux6 = d.*x(3);

oderhs = [flux1-flux4-flux5;flux4-flux3;flux2-flux6];
ode = casadi.Function('ode',{x,p},{oderhs});

% build custom RK4 integrator
dt = .1;
k1 = ode(x,p);
k2 = ode(x+dt/2.0*k1,p);
k3 = ode(x+dt/2.0*k2,p);
k4 = ode(x+dt*k3,p);

states_final = x+dt/6.0*(k1+2*k2+2*k3+k4);

one_step = casadi.Function('one_step',{x,p},{states_final});
xstate = one_step(x,p);
one_sample = casadi.Function('one_sample',{x,p},{xstate});
one_sample = one_sample.expand();

all_samples = one_sample.mapaccum('all_samples',3000);

% run to get time course solution
pert_p = opts.odep(2:end)';
pert_p(10) = 2;
get_pert_data = all_samples(no_noise_sol(4).xss,repmat(pert_p,1,3000));
% all_samples(opts.x0,repmat(opts.odep(2:end)',1,3000))
scale = ones(16,1);
scale(1) = .01;
scale(2) = 1e6;
x_sym = all_samples(no_noise_sol(4).xss,repmat(p.*scale,1,3000));

% obj = sum(sum((get_pert_data - x_sym).^2,2));
variability = (get_pert_data - x_sym);
obj = .5*dot(variability,variability);

nlp = struct('x',p,'f',obj);
solver = casadi.nlpsol('solver','ipopt',nlp);
sol_opt = solver('x0',opts.odep(2:end)');

% all step integrator
problem = struct('x',x,'ode',oderhs,'p',p);
solver_opts = opts.solver_opts;
% solver_opts.tf = 300;
solver_opts.output_t0 = 1;
solver_opts.print_stats = 1;
solver_opts.grid = tspan;
% options = struct('output_t0',1,'print_stats',1,'grid',tspan);
Fint = casadi.integrator('Fint','cvodes',problem,solver_opts);
sol = Fint('x0',x,'p',p);
yout = full(sol.xf);

all_steps = casadi.Function('all_steps',{x,p},{yout});


