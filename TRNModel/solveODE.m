% function [Solution,status,indx,nfcall,totaltfinish,itime] =...
%   test5(trnmodel,defparval,ngene,tnreg,ngap,InConc,tmax,scalAbsTol,initval,tout)
%**************************************************************************
%Wrapper to solve ODE model using SUNDIALS
%December 2013 - version 4.0
%March 24 2014
%Changed from non-stiff to a stiff solver to test for v6 of
%singlepromoteractivity.m
%April 01 2014 - Renamed to solveODE from test5
%**************************************************************************
function [Sol,status,indx,nfcall,tfinish] =...
    solveODE(initval,data,model,FBAmodel,Options)   

indx = 0;
nvar = data.nvar;
ng = data.ng;
if nargin < 5
    Options = struct();
    %Options.tmax = 8000;    
    Options.scalAbsTol = 1e-7;
    Options.RelTol = 1e-8;
    Options.MaxIter = 1000;
    Options.MaxDataPoints = 200;
elseif nargin >= 5
    if ~isfield(Options,'tout')
        Options.tout = 0.01;
    end
    if ~isfield(Options,'tmax')
        Options.tmax = 10;%h
    end
    if ~isfield(Options,'scalAbsTol')
        Options.scalAbsTol = 1e-7;
    end
    if ~isfield(Options,'RelTol')
        Options.RelTol = 1e-8;
    end
    if ~isfield(Options,'MaxIter')
        Options.MaxIter = 1000;
    end  
    if ~isfield(Options,'MaxDataPoints')
        Options.MaxDataPoints = 200;
    end
end

%Set ScalAbsTol
t0 = 0.0;
AbsTol = zeros(nvar,1);
AbsTol(1:ng(1)) = Options.RabsTol;
AbsTol(ng(1)+1:ng(1)+ng(2)+ng(3)) = Options.PabsTol;
AbsTol(ng(1)+ng(2)+ng(3)+1:end) = Options.MabsTol;
AbsTol(AbsTol==0) = Options.MabsTol;

%Initialize CVode
tstart = tic;
options = CVodeSetOptions('UserData',data,...                                                    
                          'RelTol',Options.RelTol,...
                          'AbsTol',AbsTol,...
                          'MaxNumSteps',Options.MaxIter);                             
% mondata.sol = false;
% mondata.select = find(indx);
mondata.mode = 'text';
mondata.updt = 100;
mondata.skip = 10;

%Using SUNDIALS provided monitor function  
% options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,...
%                           'MonitorData',mondata);                      
%ODE Function Called through Anonymous function                      
callmatbalance = @(t,Y,data)integrated_ODEmodel(t,Y,data,model,FBAmodel); 
%Stiff Solver and Method
status = CVodeInit(callmatbalance,'BDF','Newton',t0,initval,options);
%FSA Initialization
nvar = model.nvar;
Ns = size(model.allpar,1);
pScale = zeros(Ns,1);
pScale(1) = 1e9;
pScale(2) = 1e9;
pScale(3:end) = 1e9;
% FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
%                                  'ParamField','par',...
%                                  'ParamList',[1:Ns],...
%                                  'ParamScales',pScale,...
%                                  'DQType','Forward');  

% YS0 = zeros(nvar,Ns);
% YS0(1,1) = 1;
% YS0(1,2) = 1;
% for i=2:ng(1)
%     YS0(i,i+1) = 1;
% end
% CVodeSensInit(Ns,[],YS0,FSAoptions); 

Sol.t = [];
Sol.y = [];
Sol.ys = [];
nfcall = 0;

starttime = Options.tout;
endtime = Options.tmax;
maxdata = Options.MaxDataPoints;

itime = [];
tout = starttime:endtime/(maxdata-1):endtime;
fprintf('\nODE Solver Statistics\n');
fprintf('time\t\t\tnst\t\t\tnfe\t\t\tncfn\t\tqcur\t\thcur\n');
oldtime = 1;
for ktime = 1:length(tout)
    itstart = tic;
    [status,t,dY] = CVode(tout(ktime),'Normal');
    if status < 0
        tfinish = toc(itstart);
        return
    end
    %Print Sover Statistics
    if ktime == 1 || ktime >= oldtime+50
        stats = CVodeGetStats;
        fprintf('%2.3g\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%2.3g\n',...
        stats.tcur,stats.nst,stats.nfe,stats.ncfn,stats.qcur,stats.hcur);
        oldtime = ktime;
    end
    %Collect Solution
    Sol.t = [Sol.t;tout(ktime)];
    Sol.y = [Sol.y, dY];
%     Sol.ys = [Sol.ys, sY];        
end
stats = CVodeGetStats;
fprintf('%2.3g\t\t%d\t\t%d\t\t%d\t\t%d\t\t%2.3g\n\n',...
         stats.tcur,stats.nst,stats.nfe,stats.ncfn,stats.qcur,stats.hcur);
EWT = CVodeGet('ErrorWeights');
CVodeFree;
tfinish = toc(tstart);                          
end
