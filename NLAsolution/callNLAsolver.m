function [sol,finalSS,status] = callNLAsolver(model,pvec,initval,solverP,sol)
if nargin<5
    sol = struct([]);
end
nvar = length(initval);

% initialize KINSol
options = KINSetOptions('FuncNomTol',solverP.FuncNormTol,...
                        'ScaledStepTol',solverP.ScaledStepTol,...
                        'KrylovMaxDim',solverP.KrylovMaxDim,...
                        'MaxNumRestarts',solverP.MaxNumRestarts,...
                        'MaxNumSteps',solverP.MaxNumSteps,...                        
                        'PrecSetupFn',solverP.PrecSetupFn,...
                        'PrecSolveFn',solverP.PrecSolveFn,...
                        'Verbose',solverP.Verbose,...
                        'LinearSolver','GMRES',...
                        'UserData',data);
                    
callRHSfunction = @(y,data)ToyODEmodel(0,y,data,model,pvec);
KINInit(callRHSfunction,nvar,options);

% solve NLAE
[status,y] = KINSol(initval,'None',yscale,fscale);
si = KINGetStats;
ls_stats = si.LSInfo;

if status < 0
    fprintf('KINSOL failed. status = %d\n',status);
else
    fprintf('KINSOL succeded. status = %d\n',status);
    fprintf('%d   %d   %d   %6.2f   %6.2f',si.nfe,si.nne,si.nbops,si.fnorm,si.step);
end

KINFree;

                            
