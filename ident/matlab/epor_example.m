% redo the EPoR ligand-binding and ligand trafficking example from Raue et al., 2010

ystate_noisy = gen_noisy_example();

x = casadi.SX.sym('x',6,7501);
% y = casadi.SX.sym('y',2,1);
p = casadi.SX.sym('p',9,1);

% fluxes
v1 = p(1).*x(1,:).*x(2,:);
v2 = p(1).*p(2).*x(3,:);
v3 = p(3).*p(4);
v4 = p(3).*x(2,:);
v5 = p(5).*x(3,:);
v6 = p(6).*x(4,:);
v7 = p(7).*x(4,:);
v8 = p(8).*x(4,:);

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

dx = [dx1;dx2;dx3;dx4;dx5;dx6];
y = [p(9).*(x(1,:)+x(6,:));p(9).*x(3,:)];

xode = casadi.SX.sym('xode',6,1);
% fluxes
v1ode = p(1).*xode(1,:).*xode(2,:);
v2ode = p(1).*p(2).*xode(3,:);
v3ode = p(3).*p(4);
v4ode = p(3).*xode(2,:);
v5ode = p(5).*xode(3,:);
v6ode = p(6).*xode(4,:);
v7ode = p(7).*xode(4,:);
v8ode = p(8).*xode(4,:);
% odes
% Epo
dx1ode = -v1ode+v2ode+v6ode;
% EpoR
dx2ode = -v1ode+v2ode+v3ode-v4ode+v6ode;
% Epo_EpoR
dx3ode = v1ode-v2ode-v5ode;
% Epo_EpoR_i
dx4ode = v5ode-v6ode-v7ode-v8ode;
% dEpo_i
dx5ode = v7ode;
% dEpo_e
dx6ode = v8ode;
dxode = [dx1ode;dx2ode;dx3ode;dx4ode;dx5ode;dx6ode];

tspan = 0:.1:750;
ode = struct('x',xode,'ode',dxode,'p',p);
solver_opts.output_t0 = 1;
% solver_opts.print_stats = 1;
solver_opts.grid = tspan;
int = casadi.integrator('int','cvodes',ode,solver_opts);

CVODES_int = casadi.Function('CVODES_int',{xode,p},{int});



odeoutput = casadi.Function('odef',{x,p},{y});

ystate_sym = odeoutput(x,p);
xerror = (ystate_noisy-ystate_sym);
obj = .5*dot(xerror,xerror);

nlp = struct('x',p,'f',obj);
solver = casadi.nlpsol('solver','ipopt',nlp);





% p0_log = [-4.091;2.583;-1.758;2.821;-1.177;-2.447;-2.73;-1.884;2.21];
% p0 = 10.^(p0_log);
% x0 = [10^(3.322);10^(2.821);0;0;0;0];

% sol = int('x0',x0,'p',p0);
% yout = full(sol.xf);
% ystate_noisy = zeros(2,length(tspan));
for i = 1:length(tspan)
    ystate_noisy(:,i) = full(odeoutput(yout(:,i),p));
end
% figure
% plot(tspan,yout);

xerror = casadi.DM(ystate_noisy) - yout

odeoutput(x,p)












