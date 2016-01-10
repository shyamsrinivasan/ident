%function to solve ODE to obtain concentration profiles using SUNDIALS
function [Sol,finalSS,status] = callODEsolver(model,pvec,initval,solverP,sol)
if nargin < 5
    sol = struct([]);
end
%Initial Value
% nint_metab = model.nint_metab;
% next_metab = model.next_metab;
% nvar = model.nt_metab;%All Metabolites
% data.flux = zeros(size(model.S,2),1);
% %Setting initial values for ODE
% initval = zeros(nvar,1);
% if ~isemptyr(initialSol)
%     initval = initialSol.y(:,end);
% else
% %     bm_ind = model.bmrxn;
%     Mbio = strcmpi('biomass',model.mets);
% %     Mbio = model.S(:,bm_ind)>0;
% %     nint_metab = nint_metab-length(find(Mbio));
%     initval(1:(nint_metab-1)) = variable.MC(1:(nint_metab-1));
%     initval(Mbio) = 0;
%     initval(nint_metab+1:nint_metab+next_metab) = ...
%     assign_extconc(batch.init{1},batch.init{2},model);
% end
% [initval,initflux] = initializeConcentration(model,pmeter,variable,batch.init{1},...
%                                              batch.init{2},1,initialSol);
% initval = initConcentration(model,batch,variable,1,batch.init{1},...
%                             batch.init{2},initialSol);
% initval = sample_metabolites(model,initval)*1e-3;                    
% % initval = zeros(model.nt_metab,1);  
% initval(initval==0) = 0.001;
% initflux = initFlux(model,pmeter,initval,1);
% 
% data.flux = initflux;                                        
% initval(nint_metab+1:nt_metab) = model.M*variable.MC;

%Set ScalAbsTol
AbsTol = zeros(model.nt_metab,1);
AbsTol(AbsTol==0) = solverP.MabsTol;
data = struct();
% options = odeset('RelTol',solverP.RelTol,'AbsTol',AbsTol);

%Time vector
Sol = struct();
if isempty(sol)
    %time for initial data point
    t0 = 0.0;
%     tout = solverP.tmax;
    tout = solverP.tout:solverP.tmax/(solverP.MaxDataPoints-1):solverP.tmax;
    
    Sol.t = [];%zeros(length(tout)+1,1);
    Sol.y = zeros(length(initval),0);
    Sol.flux = zeros(size(model.S,2),0);
    itime = 0;
else
    tavail = sol.t;    
    %time for initial data point
    t0 = tavail(end);
%     tout = solverP.tmax;
    tout = tavail(end)+solverP.tout:solverP.tmax/(solverP.MaxDataPoints-1):solverP.tmax;
    
    %initial value
    initval = sol.y(:,end);    
    Sol.t = sol.t;%zeros(length(tout),1)];
    Sol.y = sol.y;% zeros(length(initval),length(tout))];
    Sol.flux = sol.flux;% zeros(size(model.S,2),length(tout))];
    itime = length(sol.t);
end
%Normalize t vectors
%t = mu*t;
% mu = model.Vss(model.bmrxn);
% tout = mu*tout;
% t0 = mu*t0;
    
%% %Initialize CVode
totaltstart = tic;
options = CVodeSetOptions('UserData',data,...                                                    
                          'RelTol',solverP.RelTol,...
                          'AbsTol',AbsTol,...
                          'MaxNumSteps',solverP.MaxIter);

%Using SUNDIALS provided monitor function  
mondata.mode = 'text';
mondata.update = 100;
mondata.skip = 10;
% mondata.initialized = false;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,...
                          'MonitorData',mondata); 
%ODE Function Called through Anonymous function                      
callODEmodel = @(t,Y,data)ODEmodel(t,Y,data,model,pvec);
CVodeInit(callODEmodel,'BDF','Newton',t0,initval,options);
% [t,dy] = ode15s(callODEmodel,tout,initval,options);
%% %Solve ODE

t = t0;
itime = 1;
%store initial values
Sol.t = [Sol.t;t];    
Sol.y = [Sol.y initval];
Sol.flux = [Sol.flux iflux(model,pvec,initval.*model.imc)];

tstep = 0;
fprintf('Initial time %4.3g\n',t0);
fprintf('Final time %4.3g\n',tout(end));
fprintf('Total simulaion time: %4.3g\n',tout(end)-t0);
% while t < tout        
    itstart = tic;
%     [status,t,dY] = CVode(tout,'OneStep');  
    [status,t,dY] = CVode(tout,'Normal');  
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
    flux = iflux(model,pvec,dY.*repmat(model.imc,1,length(tout)));
    Sol.flux = [Sol.flux flux];
%     plotflux_timecourse(flux,t,model)
%     plotconc_timecourse(dY,t,model)
    %Collect Solution
    
    
%     %Solution.ys = [Solution.ys, YS];
    
    itfinish = toc(itstart);
%     itime = itime+1;
%     itime = [itime;itfinish];
% end
% EWT = CVodeGet('ErrorWeigths');
CVodeFree;
totaltfinish = toc(totaltstart);  

%Get final Steady State
%Obtaining final SS values for initial problem 
f_flag = 0;
i = size(Sol.y,2);
while ~f_flag
% for i = size(Solution.initSS.y,2):-1:1
    if any(Sol.y(1:model.nint_metab,i)<-solverP.MabsTol)
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

%plot for debugging purposes only
figure
plot(Sol.t,Sol.y(1:48,:));
return;

