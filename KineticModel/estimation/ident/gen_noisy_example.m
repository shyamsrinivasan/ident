function ystate_noisy = gen_noisy_example()
x = casadi.SX.sym('x',6,1);
% y = casadi.SX.sym('y',2,1);
p = casadi.SX.sym('p',9,1);

% fluxes
v1 = p(1).*x(1).*x(2);
v2 = p(1).*p(2).*x(3);
v3 = p(3).*p(4);
v4 = p(3).*x(2);
v5 = p(5).*x(3);
v6 = p(6).*x(4);
v7 = p(7).*x(4);
v8 = p(8).*x(4);

% odes
% Epo
dx1 = -v1+v2+v6;
% EpoR
dx2 = -v1+v2+v3-v4+v6;
% Epo_EpoR
dx3 = v1-v2-v5;
% Epo_EpoR_i
dx4 = v5-v6-v7-v8;
% dEpo_i
dx5 = v7;
% dEpo_e
dx6 = v8;

y = [p(9).*(x(1)+x(6));p(9).*x(3)];

tspan = 0:.1:750;
odeoutput = casadi.Function('odef',{x,p},{y});

dx = [dx1;dx2;dx3;dx4;dx5;dx6];
ode = struct('x',x,'ode',dx,'p',p);
solver_opts.output_t0 = 1;
solver_opts.print_stats = 1;
solver_opts.grid = tspan;
p0_log = [-4.091;2.583;-1.758;2.821;-1.177;-2.447;-2.73;-1.884;2.21];
p0 = 10.^(p0_log);
x0 = [10^(3.322);10^(2.821);0;0;0;0];
int = casadi.integrator('int','cvodes',ode,solver_opts);
sol = int('x0',x0,'p',p0);
yout = full(sol.xf);
figure
plot(tspan,yout);

% if nsmp>1    
    pd = makedist('Uniform','lower',-.05,'upper',.05);    
% else
%     pd = makedist('Uniform','lower',-.05,'upper',.05);    
% end
met_noise = random(pd,6,1);
noisy_yexp = yout.*(1+met_noise);

ystate_noisy = zeros(2,length(tspan));
for i = 1:length(tspan)
    ystate_noisy(:,i) = full(odeoutput(noisy_yexp(:,i),p0));
end
