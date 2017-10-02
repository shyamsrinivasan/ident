% load noisy data
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end

if status == 1    
    load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noisy_model/pdata_oct2');
elseif status == 2    
    load('~/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noisy_model/pdata_oct2');    
end

% % collect only needed perturbations for analysis
% avail_pert = size(noisy_sol{1},2);
% use_pert = 1;
% npert = length(use_pert);
% [exp_select_sol,noisy_select_sol] = parseperturbations(noisy_sol{1},use_pert);
% 
% % use wt as initial value for all perturbations
% xinit = repmat(noisy_xss(:,1),npe1rt,1);
% yinit = repmat(noisy_fss(:,1),npert,1);
% get rhs fun (ode is casadi fun)
[ode,flux,D2FX,oderhs,x,p_var,p_other,acetate] = kotteCASident(3);

npts = 10;
% RK4 integrator
dt = .1;
k1 = ode(x,p_var,p_other,acetate);
k2 = ode(x+dt/2.0*k1,p_var,p_other,acetate);
k3 = ode(x+dt/2.0*k2,p_var,p_other,acetate);
k4 = ode(x+dt*k3,p_var,p_other,acetate);
xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = flux(xfinal,p_var,p_other,acetate);

% Create a function that simulates one step propagation in a sample
xstate_one_step =...
casadi.Function('xstate_one_step',{x,p_var,p_other,acetate},{xfinal,yfinal});

xstate = x;
for i=1:3000
    [xstate,ystate] = xstate_one_step(xstate,p_var,p_other,acetate);
end

% Create a function that simulates all step propagation on a sample
xstate_onepoint =...
casadi.Function('xstate_onepoint',{x,p_var,p_other,acetate},{xstate,ystate});
xstate_onepoint = xstate_onepoint.expand();

xdyn_fun = xstate_onepoint.mapaccum('all_samples',npts);

% Choose an excitation signal
acetate_data = [.1;.2;.4;.6;.8;1;1.2;1.4;1.6;2];
% acetate_data = 2;

xinit = repmat(noisy_xss(:,1),1,1);
scale = ones(9,1);
scale(3) = 1e6;
[x_measured,y_measured] =...
xdyn_fun(xinit,...
         repmat([odep_bkp(1:5)';odep_bkp(10:13)']./scale,1,npts),...
         repmat(odep_bkp(6:9)',1,npts),...
         acetate_data);

y_data = y_measured;

[x_sym,y_sym] =...
xdyn_fun(xinit,...
         repmat(p_var,1,npts),...
         repmat(odep_bkp(6:9)',1,npts),...
         acetate_data);
     
e = y_data([1 3 4 5],:) - y_sym([1 3 4 5],:);
nlp = struct('x', p_var, 'f', 0.5*dot(e,e));
solver = casadi.nlpsol('solver','ipopt', nlp);

p0 = ones(9,1);
p0(3) = 1e6;
sol = solver('x0',p0./scale);
opt_x = full(sol.x.*scale);


