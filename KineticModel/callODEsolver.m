% function to solve ODE to obtain concentration profiles using SUNDIALS
function [Sol,finalSS,status] = callODEsolver(model,pvec,initval,solverP,sol)
if nargin < 5
    sol = struct([]);
end
nvar = length(initval);

% Set Scalar Absolute Tolerance
AbsTol = zeros(nvar,1);
AbsTol(AbsTol==0) = solverP.MabsTol;
data = struct();
% options = odeset('RelTol',solverP.RelTol,'AbsTol',AbsTol);

% Time vector
Sol = struct();
if isempty(sol)
    % time for initial data point
    t0 = 0.0;
%     tout = solverP.tmax;
%     tout = solverP.tout:solverP.tmax/(solverP.MaxDataPoints-1):solverP.tmax;
    tout = solverP.tout:0.1:solverP.tmax;
    
    Sol.t = [];%zeros(length(tout)+1,1);
    Sol.y = zeros(nvar,0);
    Sol.flux = zeros(size(model.S,2),0);
    itime = 0;
else
    tavail = sol.t;    
    % time for initial data point
    t0 = tavail(end);
%     tout = solverP.tmax;
%     tout = tavail(end)+solverP.tout:solverP.tmax/(solverP.MaxDataPoints-1):solverP.tmax;
    tout = tavail(end)+solverP.tout:0.1:solverP.tmax;
    % initial value
    initval = sol.y(:,end);    
    Sol.t = sol.t;%zeros(length(tout),1)];
    Sol.y = sol.y;% zeros(length(initval),length(tout))];
    Sol.flux = sol.flux;% zeros(size(model.S,2),length(tout))];
    itime = length(sol.t);
end
% Normalize t vectors
% t = mu*t;
% mu = model.Vss(model.bmrxn);
% tout = mu*tout;
% t0 = mu*t0;
    
%% %Initialize CVode
totaltstart = tic;
options = CVodeSetOptions('UserData',data,...                                                    
                          'RelTol',solverP.RelTol,...
                          'AbsTol',AbsTol,...
                          'MaxNumSteps',solverP.MaxIter);

% Using SUNDIALS provided monitor function  
mondata.mode = 'graphical';
mondata.updt = 100;
mondata.skip = 10;
% mondata.initialized = false;
% options = CVodeSetOptions(options,'MonitorFn',@MonitorCVodeGraph,...
%                           'MonitorData',mondata); 
% ODE Function Called through Anonymous function  
% ecoli model
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,data,model,pvec);

% toy model
callODEmodel = @(t,Y,data)ToyODEmodel(t,Y,data,model,pvec);
CVodeInit(callODEmodel,'Adams','Newton',t0,initval,options);
% [t,dy] = ode15s(callODEmodel,tout,initval,options);
%% %Solve ODE

t = t0;
itime = 1;
%store initial values
Sol.t = [Sol.t;t];    
Sol.y = [Sol.y initval];
Sol.flux = [Sol.flux iflux(model,pvec,[initval.*model.imc;model.Pimc])];

tstep = 0;
fprintf('Initial time %4.3g\n',t0);
fprintf('Final time %4.3g\n',tout(end));
fprintf('Total simulaion time: %4.3g\n',tout(end)-t0);
% while t < tout        
    itstart = tic;
%     [status,t,dY] = CVode(tout,'OneStep');   
%     si = CVodeGetStats;
%     fprintf('%4.3g\t\n',si.tcur);
%     [status,t,dY] = CVode(tout,'OneStep');  
    [status,t,dY] = CVode(tout,'Normal');  
%     opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
%     [t,dY] = ode15s(callODEmodel,tout,initval,opts);
%     tstep = tstep+1;
%     if ~rem(tstep,100)
%         si = CVodeGetStats;
%         fprintf('%6.5g\t %d\t %d\t %d\t %d\t\n',si.tcur,si.nst,si.nfe,si.netf,si.ncfn);                  
        Sol.t = [Sol.t;t'];
        Sol.y = [Sol.y dY];
%         Sol.flux = [Sol.flux flux];
%     end
%     [status,t,dY] = CVode(tout,'Normal');   
    %Plot Solution
    flux = iflux(model,pvec,[dY.*repmat(model.imc,1,length(tout));...
                             repmat(model.Pimc,1,length(tout))]);
    Sol.flux = [Sol.flux flux];
%     plotflux_timecourse(flux,t,model)
%     plotconc_timecourse(dY,t,model)
    % Collect Solution
    
    
%     %Solution.ys = [Solution.ys, YS];
    
    itfinish = toc(itstart);
%     itime = itime+1;
%     itime = [itime;itfinish];
% end
% EWT = CVodeGet('ErrorWeigths');
% CVodeFree;
totaltfinish = toc(totaltstart);  

% Get final Steady State
% Obtaining final SS values for initial problem 
f_flag = 0;
i = size(Sol.y,2);
while ~f_flag
% for i = size(Solution.initSS.y,2):-1:1
    if any(Sol.y(1:nvar,i)<-solverP.MabsTol)
        finalSS.t = Sol.t(i-1);
        finalSS.y = Sol.y(:,i-1);        
        finalSS.flux = Sol.flux(:,i-1);
        f_flag = 1;  
    else
        finalSS.t = Sol.t(i);
        finalSS.y = Sol.y(:,i);        
        finalSS.flux = Sol.flux(:,i);
        f_flag = 1;
    end
    i = i-1;
end
if ~f_flag
    finalSS.t = Sol.t(end);
    finalSS.y = Sol.y(:,end); 
    finalSS.flux = Sol.flux(:,end);
end

% plot
% rmvlst = model.mets(model.nint_metab+1:model.nt_metab);
% rmvlst = {rmvlst{:},'h2o[c]','h[c]'}; 
% var.mets = removemets(model.mets,rmvlst);
var.mets = {'pep[c]','pyr[c]','atp[c]','adp[c]','g6p[c]',...
            'q8[c]','q8h2[c]','nad[c]','nadh[c]'};
timecourseplots(Sol.t,Sol.y,'c',model,var);
var1.rxns = model.rxns([1 2 3 4]);
timecourseplots(Sol.t,Sol.flux,'f',model,var1);
% plot for debugging purposes only
return;

