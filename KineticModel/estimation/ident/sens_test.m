% test sensitivity calculation w/ state equations using direct method
% load original kotte model
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end
if status == 1
    load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
elseif status == 2
    load('/home/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/Kotte2014/model/kotte_model.mat');    
end

pconv = [.1,.3,0,pvec(9)]; % extra parameters for CK 'K1pep','K2fdp','rhoA' and 'acetate'
p = [pvec,pconv];
p(9) = [];
ival = [M;zeros(3*17,1)];
clear pvec

[FXaug,FX] = kotteCASwSENS(3,17);

% solve deterministic (no noise) model using casadi(cvodes)
odep = p;
fh = @()kotteCASwSENS(3,17); % @kotte_CAS
solver_opts = struct('abstol',1e-6,'reltol',1e-6);
opts = struct('tspan',0:.1:300,'x0',ival,'solver_opts',solver_opts,'odep',odep);
[xdyn,fdyn,xss1,fss1] = solveODE_cas(fh,opts);
figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout','vEout');