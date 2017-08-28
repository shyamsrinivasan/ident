% estimationn script - simultaneous ss estimation of fluxes, concenrations
% and parameters
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
kotte_pest_allf_typeb_obj(xexp,vexp,nc,nf,npert);
% sym cons
cons =...
kotte_pest_allf_typeb_cons(xexp,vexp,nc,nf,npert,x,p,flux,vareps,p_usl,ac);
ncons = size(cons,1);
consfun = casadi.Function('consfun',{var,par},{cons});

% get parameters
p0_obj = odep_bkp(1:13)';
scale = ones(np,1);
scale(3) = 1e-6;

weigths = [1;1;1;100];
p0_const = odep_bkp(14:end)';
par_val = [p0_const;weigths];

% get initial value
xval = [xexp;p0_obj.*scale;vexp];
x0 = [xval;0;0];
cons_check = full(consfun(x0,par_val));

% set bounds for x0
lb = zeros(nvar,1);
ub = zeros(nvar,1);
% conc
lb(1:nc*npert) = 1e-7;
ub(1:nc*npert) = 300;
% parameters
[p_bounds_lb,p_bounds_ub] = p_bounds();
lb(nc*npert+1:nc*npert+np) = p_bounds_lb;
ub(nc*npert+1:nc*npert+np) = p_bounds_ub;
% fluxes
lb(nc*npert+np+1:nc*npert+np+nf*npert) = 0;
ub(nc*npert+np+1:nc*npert+np+nf*npert) = 20;
% uncertainty
lb(nvar-2+1:nvar-1) = 0;
ub(nvar-2+1:nvar-1) = 1;
lb(nvar-1+1:nvar) = 0;
ub(nvar-1+1:nvar) = 1;

% set bounds for cons
lbg = ones(ncons,1);
ubg = ones(ncons,1);
% ss cons bounds (=) type a 
% lbg(1:nc*npert) = 0;
% ubg(1:nc*npert) = 0;
% % nl flux cons (=)
% lbg(nc*npert+1:nc*npert+nf*npert) = 0;
% ubg(nc*npert+1:nc*npert+nf*npert) = 0;
% % nl noisy conc cons (<=)
% lbg(nc*npert+nf*npert+1:nc*npert+nf*npert+2*nc*npert) = -Inf;
% ubg(nc*npert+nf*npert+1:nc*npert+nf*npert+2*nc*npert) = 0;
% % nl noisy flux cons (<=)
% lbg(nc*npert+nf*npert+2*nc*npert+1:nc*npert+nf*npert+2*nc*npert+2*nf*npert) = -Inf;
% ubg(nc*npert+nf*npert+2*nc*npert+1:nc*npert+nf*npert+2*nc*npert+2*nf*npert) = 0;

% set flux cons bounds type b
lbg(1:nf*npert) = 0;
ubg(1:nf*npert) = 0;
% nl noisy conc cons (<=)
lbg(nf*npert+1:nf*npert+2*nc*npert) = -Inf;
ubg(nf*npert+1:nf*npert+2*nc*npert) = 0;
% nl noisy flux cons (<=)
lbg(nf*npert+2*nc*npert+1:nf*npert+2*nc*npert+2*nf*npert) = -1;
ubg(nf*npert+2*nc*npert+1:nf*npert+2*nc*npert+2*nf*npert) = 0;

estprob = struct('x',var,'f',obj,'g',cons,'p',par);
options.warn_initial_bounds = 1;
options.ipopt.fixed_variable_treatment = 'make_constraint';
solver = casadi.nlpsol('solver','ipopt',estprob,options);
sol = solver('x0',x0,'p',par_val,'lbx',lb,'ubx',ub,'lbg',lbg,'ubg',ubg);
optsol_sim.xval = full(sol.x);
optsol_sim.fval = full(sol.f);
optsol_sim.exitflag = 1;

%% data processing and figures
data = struct('nc',nc,'nf',nf,'nflx',nf,'nvar',nvar,'npert',npert,...
              'np',np,'type',2,'odep',odep_bkp,'p_id',1:13,...
              'pscale',scale,'flxid',1);
% exp_data = exp_sol{1};
opts.tspan = 0:.1:600;
[proc_data,exp_data] = recalcss(optsol_sim,noisy_sol{1},[],data,opts);

%% figures
exp_xss = cat(2,noisy_sol{1}.xss);
pep = [proc_data.opt_xss(1,:)' proc_data.calc_xss(1,:)' exp_xss(1,:)'];
fdp = [proc_data.opt_xss(2,:)' proc_data.calc_xss(2,:)' exp_xss(2,:)'];
enz = [proc_data.opt_xss(3,:)' proc_data.calc_xss(3,:)' exp_xss(3,:)'];

exp_fss = cat(2,noisy_sol{1}.fss);
f1 = [proc_data.opt_fss(1,:)' proc_data.calc_fss(1,:)' exp_fss(1,:)'];
f3 = [proc_data.opt_fss(3,:)' proc_data.calc_fss(3,:)' exp_fss(3,:)'];
f4 = [proc_data.opt_fss(4,:)' proc_data.calc_fss(4,:)' exp_fss(4,:)'];
f5 = [proc_data.opt_fss(5,:)' proc_data.calc_fss(5,:)' exp_fss(5,:)'];

% bar plot
hf1 = figure;
bar(pep);
ylabel('pep a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');
hf2=figure;
bar(fdp);
ylabel('fdp a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');
hf3=figure;
bar(enz);
ylabel('E a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');

hf4 = figure;
bar(f1);
ylabel('v1 a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');
hf5=figure;
bar(f3);
ylabel('v2 a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');
hf6=figure;
bar(f4);
ylabel('v3 a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');
bar(f5);
ylabel('v4 a.u.');
legend('Optimal Value','Calculated Value','Noisy Data');

% scatter plots
hf7 = figure;
hold on
plot(pep(:,end),pep(:,end),'LineWidth',1.5,'Color','k');
plot(pep(:,end),pep(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(pep(:,end),pep(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental pep a.u.');
ylabel('Optimal/ODE Sim pep a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

hf8 = figure;
hold on
plot(fdp(:,end),fdp(:,end),'LineWidth',1.5,'Color','k');
plot(fdp(:,end),fdp(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(fdp(:,end),fdp(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental fdp a.u.');
ylabel('Optimal/ODE Sim fdp a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

hf9=figure;
hold on
plot(enz(:,end),enz(:,end),'LineWidth',1.5,'Color','k');
plot(enz(:,end),enz(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(enz(:,end),enz(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental E a.u.');
ylabel('Optimal/ODE Sim E a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

hf10=figure;
hold on
plot(f1(:,end),f1(:,end),'LineWidth',1.5,'Color','k');
plot(f1(:,end),f1(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(f1(:,end),f1(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental v1 a.u.');
ylabel('Optimal/ODE Sim v1 a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

hf11=figure;
hold on
plot(f3(:,end),f3(:,end),'LineWidth',1.5,'Color','k');
plot(f3(:,end),f3(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(f3(:,end),f3(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental v2 a.u.');
ylabel('Optimal/ODE Sim v2 a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

hf12=figure;
hold on
plot(f4(:,end),f4(:,end),'LineWidth',1.5,'Color','k');
plot(f4(:,end),f4(:,1),'LineStyle','none','Marker','.','MarkerSize',20,'Color','b');
plot(f4(:,end),f4(:,2),'LineStyle','none','Marker','.','MarkerSize',12,'Color','r');
xlabel('Experimental v3 a.u.');
ylabel('Optimal/ODE Sim v3 a.u.');
legend('Noisy Data Line','Optimal Value','Calculated Value');

%% save figure files
dir = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\estimation\simest\typeb\';
set(0,'CurrentFigure',hf1);
fname = 'simest_noisy_pep_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf1);
set(0,'CurrentFigure',hf2);
fname = 'simest_noisy_fdp_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf2);
set(0,'CurrentFigure',hf3);
fname = 'simest_noisy_enz_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf3);
set(0,'CurrentFigure',hf4);
fname = 'simest_noisy_v1_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf4);
set(0,'CurrentFigure',hf5);
fname = 'simest_noisy_v2_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf5);
set(0,'CurrentFigure',hf6);
fname = 'simest_noisy_v3_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf6);
set(0,'CurrentFigure',hf7);
fname = 'simest_noisy_pep_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf7);
set(0,'CurrentFigure',hf8);
fname = 'simest_noisy_fdp_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf8);
set(0,'CurrentFigure',hf9);
fname = 'simest_noisy_enz_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf9);
set(0,'CurrentFigure',hf10);
fname = 'simest_noisy_v1_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf10);
set(0,'CurrentFigure',hf11);
fname = 'simest_noisy_v2_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf11);
set(0,'CurrentFigure',hf12);
fname = 'simest_noisy_v3_scatter_aug28';
print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
print([dir fname],'-dpng','-r200');
close(hf12);










% hfcv = compare_vals(proc_data,noisy_sol{1},[],data,1);

% xopt = full(sol.x);


% x = casadi.SX.sym('x',1);
% y = casadi.SX.sym('y',1);
% nlp = struct('x',[x;y],'f',x/y);
% options.warn_initial_bounds = 1;
% NLP = casadi.nlpsol('solver','ipopt',nlp,options);
% sol = NLP('x0',[.1 .1],'lbx',[.001 .001],'ubx',[10 10]);
% int = struct('x',x,'ode',exp(-x));
% opt.grid = 0:.1:20;
% F = casadi.integrator('F','cvodes',int,opts);




