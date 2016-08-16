function [sol,finalSS,status] = callNLAsolver(model,pvec,initval,solverP,sol)
if nargin<5
    sol = struct([]);
end
nvar = length(initval);
data = [];

% initialize KINSol
% options = KINSetOptions('FuncNormTol',solverP.FuncNormTol,...
%                         'ScaledStepTol',solverP.ScaledStepTol,...                        
%                         'MaxNumSetups',solverP.MaxNumSetups,...
%                         'Verbose',solverP.Verbose,...
%                         'LinearSolver','Dense',...
%                         'UserData',data);
options = optimset('Display','iter',...
                    'TolFun',1e-13,...
                   'TolX',1e-13,...
                   'MaxFunEvals',500000,...
                   'MaxIter',100000);
                    
callRHSfunction = @(y)ToyNLAmodel(y,model,pvec);
% KINInit(callRHSfunction,nvar,options);

yscale = ones(nvar,1);
fscale = ones(nvar,1);

% solve NLAE
% [status,y] = KINSol(initval,'LineSearch',yscale,fscale);
[x,fval,exitflag,output] = fsolve(callRHSfunction,initval,options);

sol.y = x;
sol.fval = fval;

% si = KINGetStats;
% ls_stats = si.LSInfo;
% 
% if status < 0
%     fprintf('KINSOL failed. status = %d\n',status);
% else
%     fprintf('KINSOL succeded. status = %d\n',status);
%     fprintf('%d   %d   %d   %6.2f   %6.2f\n',si.nfe,si.nni,si.nbops,si.fnorm,si.step);
% end

% KINFree;

                            