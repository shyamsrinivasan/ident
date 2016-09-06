% function to solve ODE to obtain concentration profiles using SUNDIALS
function [outsol,outss] = callODEsolver(model,pvec,initval,solverP,inputsol)
if nargin < 5
    inputsol = struct([]);
end
nvar = length(initval);

% Set Scalar Absolute Tolerance
AbsTol = zeros(nvar,1);
AbsTol(AbsTol==0) = solverP.MabsTol;
data = struct();

% solution data strcuture
outsol = struct();
if isempty(inputsol)
    % time for initial data point
    t0 = 0.0;
    tspan = solverP.tout:0.1:solverP.tmax;    
    outsol.t = zeros(length(tspan)+1,1);
    outsol.y = zeros(nvar,length(tspan)+1);
    outsol.flux = zeros(size(model.S,2),length(tspan)+1);    
else
    tavail = inputsol.t;    
    % time for initial data point
    t0 = tavail(end);
    tspan = tavail(end)+solverP.tout:0.1:solverP.tmax;
    % initial value
    initval = inputsol.y(:,end);    
    outsol.t = inputsol.t;
    outsol.y = inputsol.y;
    outsol.flux = inputsol.flux;    
end
% Normalize t vectors
% t = mu*t;
% mu = model.Vss(model.bmrxn);
% tout = mu*tout;
% t0 = mu*t0;

tstart = tic;

% Initialize CVode
% options = CVodeSetOptions('UserData',data,...                                                    
%                           'RelTol',solverP.RelTol,...
%                           'AbsTol',AbsTol,...
%                           'MaxNumSteps',solverP.MaxIter);

% Using SUNDIALS provided monitor function  
% mondata.mode = 'graphical';
% mondata.updt = 100;
% mondata.skip = 10;
% mondata.initialized = false;
% options = CVodeSetOptions(options,'MonitorFn',@MonitorCVodeGraph,...
%                           'MonitorData',mondata); 

% ODE Function Called through Anonymous function  
% ecoli model
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,data,model,pvec);

% toy model
callODEmodel = @(t,Y)ToyODEmodel(t,Y,data,model,pvec);
% CVodeInit(callODEmodel,'Adams','Newton',t0,initval,options);

% Solve ODE
t = t0;

% store initial values
outsol.t(1) = t;    
outsol.y(:,1) = initval;
outsol.flux(:,1) = iflux(model,pvec,[initval.*model.imc;model.PM]);

fprintf('Initial time:\t\t\t %4.3g\n',t0);
fprintf('Final time:\t\t\t\t %4.3g\n',tspan(end));
fprintf('Total simulaion time:\t %4.3g\n',tspan(end)-t0);

% get CVode stats
% si = CVodeGetStats;

% [status,t,yout] = CVode(tspan,'Normal');  
% opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
opts = odeset('RelTol',solverP.RelTol,'AbsTol',AbsTol);
[t,yout] = ode15s(callODEmodel,tspan,initval,opts);

% CVodeFree;

telaps = toc(tstart);
fprintf('Model integration time:\t %4.3g\n',telaps);

% collect solution
t = columnVector(t);
outsol.t(2:end) = t';
try
    outsol.y(:,2:end) = yout;
catch
    outsol.y(:,2:end) = yout';
end
try
    allmets = [yout.*repmat(model.imc',length(t),1),...
               repmat(model.PM',length(t),1)];
    allmets = allmets';
catch
    allmets = [yout'.*repmat(model.imc',1,length(t));...
               repmat(model.PM,1,length(t))];
end
% calculate time course fluxes
outsol.flux(:,2:end) = iflux(model,pvec,allmets);
fprintf('Flux calculation time:\t %4.3g\n',toc(tstart)-telaps);

% plot solution for diagnostics
timecourseplots(outsol.t,outsol.y,1,{'pep[c]','fdp[c]','bm[c]','pyr[c]'},model);
timecourseplots(outsol.t,outsol.flux,2,[1 3 4 5],model);

% Get final Steady State
% Obtaining final SS values for initial problem 
f_flag = 0;
i = size(outsol.y,2);
while ~f_flag
    if any(outsol.y(1:nvar,i)<-solverP.MabsTol)
        outss.t = outsol.t(i-1);
        outss.y = outsol.y(:,i-1);        
        outss.flux = outsol.flux(:,i-1);
        f_flag = 1;  
    else
        outss.t = outsol.t(i);
        outss.y = outsol.y(:,i);        
        outss.flux = outsol.flux(:,i);
        f_flag = 1;
    end
    i = i-1;
end
if ~f_flag
    outss.t = outsol.t(end);
    outss.y = outsol.y(:,end); 
    outss.flux = outsol.flux(:,end);
end

fprintf('Completion time:\t\t %4.3g\n',toc(tstart));
return;

