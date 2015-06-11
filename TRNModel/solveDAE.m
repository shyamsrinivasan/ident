%Wrapper to solve DAE model using IDA SUNDIALS
function [Sol,status,indx,nfcall,tfinish] =...
    solveDAE(y0,yp0,data,model,FBAmodel,Options) 

indx = 0;
nvar = data.nvar;
ng = data.ng;
varid = data.varind;
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

%Initialize IDA
tstart = tic;
options = IDASetOptions('UserData',data,...                                                    
                        'RelTol',Options.RelTol,...
                        'AbsTol',AbsTol,...
                        'VariableTypes',varid,...%0 - algebraic variables
                        'LinearSolver','TFQMR',...
                        'MaxNumSteps',Options.MaxIter);  
mondata.mode = 'text';
mondata.updt = 100;
mondata.skip = 10;    
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);
%ODE Function Called through Anonymous function                      
DAEmodel = @(t,Y,YP,data)integratedDAEmodel(t,Y,YP,data,model,FBAmodel); 
%Initialize IDA Solver
IDAInit(DAEmodel,t0,y0,yp0,options);

%Integrate DAE
starttime = Options.tout;
endtime = Options.tmax;
maxdata = Options.MaxDataPoints;

%Calculate consistent initial conditions for algebraic variables
[status, y0_mod, yp0_mod] = IDACalcIC(starttime, 'FindAlgebraic');

tout = starttime:endtime/(maxdata-1):endtime;
Sol.t = zeros(maxdata+1,1);
Sol.y = zeros(model.nvar,maxdata+1);
Sol.t(1) = t0;
Sol.y(:,1) = y0_mod;
fprintf('\nODE Solver Statistics\n');
fprintf('time\t\t\tnst\t\t\tnfe\t\t\tncfn\t\tqcur\t\thcur\n');
oldtime = 1;

for ktime = 1:length(tout)
    itstart = tic;
    [status,t,dY] = IDASolve(tout(ktime),'OneStep');
    if status < 0
        tfinish = toc(itstart);
        return
    end
    %Print Sover Statistics
    if ktime == 1 || ktime >= oldtime+50
        stats = IDAGetStats;
        fprintf('%2.3g\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%2.3g\n',...
        stats.tcur,stats.nst,stats.nfe,stats.ncfn,stats.qcur,stats.hcur);
        oldtime = ktime;
    end
    %Collect Solution
    Sol.t(ktime+1) = t;%[Sol.t;tout(ktime)];
    Sol.y(:,ktime+1) = dY;%[Sol.y, dY];
%     Sol.ys = [Sol.ys, sY];        
end

                    
                      
                      