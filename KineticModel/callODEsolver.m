%function to solve ODE to obtain concentration profiles using SUNDIALS
function [Sol,finalSS,status] = callODEsolver(model,pvec,initval,solverP)%Initial Value
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
%time for initial data point
t0 = 0.0;
%Set ScalAbsTol
AbsTol = zeros(model.nt_metab,1);
AbsTol(AbsTol==0) = solverP.MabsTol;
data = struct();
% options = odeset('RelTol',solverP.RelTol,'AbsTol',AbsTol);

    
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
%Time vector
tout = solverP.tout:solverP.tmax/(solverP.MaxDataPoints-1):solverP.tmax;
% [t,dy] = ode15s(callODEmodel,tout,initval,options);
%% %Solve ODE
Sol = struct();
Sol.t = zeros(length(tout),1);
Sol.y = zeros(length(initval),length(tout));
Sol.flux = zeros(size(model.S,2),length(tout));
itime = [];
for ktime = 1:length(tout)
    itstart = tic;
    [status,t,dY] = CVode(tout(ktime),'Normal');   
    %Plot Solution
    flux = iflux(model,pvec,dY);
%     plotflux_timecourse(flux,t,model)
%     plotconc_timecourse(dY,t,model)
    %Collect Solution
    Sol.t(ktime) = tout(ktime);
    Sol.y(:,ktime) = dY;
    Sol.flux(:,ktime) = flux;
    %Solution.ys = [Solution.ys, YS];
    itfinish = toc(itstart);
    itime = [itime;itfinish];
end
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
        finalSS.y = Sol.y(:,i-1);
        finalSS.t = Sol.t(i-1);
        f_flag = 1;  
    else
        finalSS.y = Sol.y(:,i);
        finalSS.t = Sol.t(i);
        f_flag = 1;
    end
    i = i-1;
end
if ~f_flag
    finalSS.t = Sol.t(end);
    finalSS.y = Sol.y(:,end); 
end
return;

